#include "stabaal/detail/inverse_tableau.hpp"

#include <algorithm>
#include <array>
#include <cassert>
#include <cstddef>
#include <initializer_list>
#include <iostream>
#include <limits>
#include <memory>
#include <random>
#include <utility>
#include <vector>

#include "stabaal/detail/bit_word.hpp"
#include "stabaal/detail/pauli_string_view.hpp"
#include "stabaal/types.hpp"

namespace stabaal::detail {

namespace {

constexpr std::size_t kBasisCount{2};
constexpr std::size_t kPauliComponents{2};
constexpr std::size_t kSwapFlags{kBasisCount};
constexpr std::size_t kFlagArrays{kBasisCount + 1};
constexpr std::size_t kStackTargetWords{16};
constexpr QubitIdx kGrowthDivisor{4};

// The complete QubitIdx range needs up to twice its bit width for storage sizes.
static_assert(std::numeric_limits<std::size_t>::digits >=
              2 * std::numeric_limits<QubitIdx>::digits);

[[nodiscard]] std::mt19937& RandomEngine() {
    static thread_local std::mt19937 engine{std::random_device{}()};
    return engine;
}

[[nodiscard]] bool Phase(const PhaseExponent exponent, const bool a_sign, const bool b_sign,
                         const PhaseExponent phase_rule) noexcept {
    assert(phase_rule >= 0 && phase_rule < 4);
    const auto phase{((exponent + phase_rule) % 4 + 4) % 4};
    assert(phase % 2 == 0);
    return (phase == 2) != (a_sign != b_sign);
}

}  // namespace

std::size_t InverseTableau::StorageWords(const QubitIdx capacity) noexcept {
    return (kPauliComponents * kBasisCount * capacity + kFlagArrays) * WordCount(capacity);
}

BitWordSpan InverseTableau::RowWords(const std::size_t row) noexcept {
    const auto stride{WordCount(capacity_)};
    return BitWordSpan{words_.get(), StorageWords(capacity_)}.subspan(row * stride, stride);
}

ConstBitWordSpan InverseTableau::RowWords(const std::size_t row) const noexcept {
    const auto stride{WordCount(capacity_)};
    return ConstBitWordSpan{words_.get(), StorageWords(capacity_)}.subspan(row * stride, stride);
}

BitWordSpan InverseTableau::Flags(const std::size_t index) noexcept {
    return RowWords(kPauliComponents * kBasisCount * capacity_ + index);
}

ConstBitWordSpan InverseTableau::Flags(const std::size_t index) const noexcept {
    return RowWords(kPauliComponents * kBasisCount * capacity_ + index);
}

std::size_t InverseTableau::PhysicalBasis(const QubitIdx q, const Basis basis) const noexcept {
    const auto swapped{static_cast<std::size_t>(GetBit(Flags(kSwapFlags), q))};
    return static_cast<std::size_t>(basis) ^ swapped;
}

PauliStringView InverseTableau::Generator(const QubitIdx q, const Basis basis) noexcept {
    return PhysicalGenerator(q, PhysicalBasis(q, basis));
}

PauliStringView InverseTableau::PhysicalGenerator(const QubitIdx q,
                                                  const std::size_t basis) noexcept {
    const auto row{kPauliComponents * (basis * capacity_ + q)};
    return {RowWords(row), RowWords(row + 1)};
}

bool InverseTableau::Sign(const QubitIdx q, const Basis basis) const noexcept {
    return GetBit(Flags(PhysicalBasis(q, basis)), q);
}

void InverseTableau::FlipSign(const QubitIdx q, const Basis basis) noexcept {
    FlipBit(Flags(PhysicalBasis(q, basis)), q);
}

InverseTableau::InverseTableau(const QubitIdx num_qubits)
    : words_{
          // NOLINTNEXTLINE(cppcoreguidelines-avoid-c-arrays)
          num_qubits == 0 ? nullptr : std::make_unique<BitWord[]>(StorageWords(num_qubits)),
      },
      num_qubits_{num_qubits},
      capacity_{num_qubits} {
    for (QubitIdx i{0}; i < num_qubits; ++i) {
        Generator(i, Basis::kComputational).SetZ(i);
        Generator(i, Basis::kHadamard).SetX(i);
    }
}

void InverseTableau::CopyRows(const InverseTableau& other) noexcept {
    const auto used_words{WordCount(other.num_qubits_)};
    if (used_words == 0) return;
    if (capacity_ == other.capacity_) {
        std::copy_n(other.words_.get(), StorageWords(capacity_), words_.get());
        return;
    }
    for (std::size_t basis{0}; basis < kBasisCount; ++basis) {
        for (QubitIdx q{0}; q < other.num_qubits_; ++q) {
            for (std::size_t axis{0}; axis < kPauliComponents; ++axis) {
                const auto source{
                    other.RowWords(kPauliComponents * (basis * other.capacity_ + q) + axis),
                };
                const auto target{RowWords(kPauliComponents * (basis * capacity_ + q) + axis)};
                std::copy_n(source.data(), used_words, target.data());
            }
        }
    }
    for (std::size_t i{0}; i < kFlagArrays; ++i) {
        std::copy_n(other.Flags(i).data(), used_words, Flags(i).data());
    }
}

InverseTableau::InverseTableau(const InverseTableau& other)
    : words_{
          other.num_qubits_ == 0
              ? nullptr
              // NOLINTNEXTLINE(cppcoreguidelines-avoid-c-arrays)
              : std::make_unique_for_overwrite<BitWord[]>(StorageWords(other.num_qubits_)),
      },
      num_qubits_{other.num_qubits_},
      capacity_{other.num_qubits_} {
    CopyRows(other);
}

InverseTableau& InverseTableau::operator=(const InverseTableau& other) {
    if (this != &other) {
        if (words_ && capacity_ == other.num_qubits_) {
            CopyRows(other);
            num_qubits_ = other.num_qubits_;
        } else {
            *this = InverseTableau{other};
        }
    }

    return *this;
}

void InverseTableau::Add() {
    assert(num_qubits_ < std::numeric_limits<QubitIdx>::max());

    if (num_qubits_ == capacity_) {
        const auto growth{std::max(capacity_ / kGrowthDivisor, QubitIdx{1})};
        const auto extra{std::min(growth, std::numeric_limits<QubitIdx>::max() - capacity_)};
        InverseTableau grown{0};
        grown.capacity_ = capacity_ + extra;
        // NOLINTNEXTLINE(cppcoreguidelines-avoid-c-arrays)
        grown.words_ = std::make_unique<BitWord[]>(StorageWords(grown.capacity_));
        grown.CopyRows(*this);
        grown.num_qubits_ = num_qubits_;
        *this = std::move(grown);
    }
    const auto last{num_qubits_++};
    Generator(last, Basis::kComputational).SetZ(last);
    Generator(last, Basis::kHadamard).SetX(last);
}

void InverseTableau::CleanUp(const QubitIdx pivot, const ConstBitWordSpan targets,
                             const bool apply_s_dagger) noexcept {
    for (std::size_t basis{0}; basis < kBasisCount; ++basis) {
        const auto signs{Flags(basis)};
        for (QubitIdx q{0}; q < num_qubits_; ++q) {
            const auto row{kPauliComponents * (basis * capacity_ + q)};
            PauliStringView gen{RowWords(row), RowWords(row + 1)};
            const bool flip{gen.CleanUp(pivot, targets, apply_s_dagger)};
            FlipBit(signs, q, flip);
        }
    }
}

void InverseTableau::Collapse(const QubitIdx pivot, const bool p) noexcept {
    for (std::size_t basis{0}; basis < kBasisCount; ++basis) {
        const auto signs{Flags(basis)};
        for (QubitIdx q{0}; q < num_qubits_; ++q) {
            const auto row{kPauliComponents * (basis * capacity_ + q)};
            PauliStringView gen{RowWords(row), RowWords(row + 1)};
            gen.SwapXZ(pivot);
            const bool flip{gen.GetZ(pivot) && (p != gen.GetX(pivot))};
            FlipBit(signs, q, flip);
        }
    }
}

bool InverseTableau::RandomMeasurement(const QubitIdx q, const Basis b) {
    const auto gen{Generator(q, b)};
    QubitIdx pivot{0};

    for (QubitIdx i{0}; i < num_qubits_; ++i) {
        if (gen.GetX(i)) {
            pivot = i;
            break;
        }
    }

    const auto num_words{WordCount(num_qubits_)};
    std::array<BitWord, kStackTargetWords> local_targets{};
    std::vector<BitWord> heap_targets(num_words > local_targets.size() ? num_words : 0);
    const BitWordSpan targets{heap_targets.empty() ? local_targets.data() : heap_targets.data(),
                              num_words};
    const bool apply_s_dagger{gen.MeasurementTargets(pivot, targets)};
    if (apply_s_dagger ||
        std::any_of(targets.begin(), targets.end(), [](const auto word) { return word != 0; })) {
        CleanUp(pivot, targets, apply_s_dagger);
    }

    const bool p{std::bernoulli_distribution{0.5}(RandomEngine())};
    Collapse(pivot, p);
    return Sign(q, b);
}

bool InverseTableau::Measure(const QubitIdx q, const Basis basis) {
    assert(q < num_qubits_);
    return Generator(q, basis).AnyX(num_qubits_) ? RandomMeasurement(q, basis) : Sign(q, basis);
}

bool InverseTableau::MeasureX(const QubitIdx q) { return Measure(q, Basis::kHadamard); }

bool InverseTableau::MeasureZ(const QubitIdx q) { return Measure(q, Basis::kComputational); }

void InverseTableau::ApplyH(const QubitIdx q) noexcept {
    assert(q < num_qubits_);

    FlipBit(Flags(kSwapFlags), q);
}

void InverseTableau::ApplyS(const QubitIdx q) noexcept {
    assert(q < num_qubits_);

    const auto computational{PhysicalBasis(q, Basis::kComputational)};
    const auto hadamard{computational ^ 1U};
    const auto computational_signs{Flags(computational)};
    const auto hadamard_signs{Flags(hadamard)};
    const bool computational_sign{GetBit(computational_signs, q)};
    const bool hadamard_sign{GetBit(hadamard_signs, q)};
    const auto exponent{
        PhysicalGenerator(q, hadamard)
            .XorWithPhase(PhysicalGenerator(q, computational), num_qubits_),
    };
    SetBit(hadamard_signs, q, Phase(exponent, hadamard_sign, computational_sign, 3));
}

void InverseTableau::ApplyX(const QubitIdx q) noexcept {
    assert(q < num_qubits_);

    FlipSign(q, Basis::kComputational);
}

void InverseTableau::ApplyZ(const QubitIdx q) noexcept {
    assert(q < num_qubits_);

    FlipSign(q, Basis::kHadamard);
}

void InverseTableau::ApplyCx(const QubitIdx c, const QubitIdx t) noexcept {
    assert(c < num_qubits_);
    assert(t < num_qubits_);

    const auto c_computational{PhysicalBasis(c, Basis::kComputational)};
    const auto t_computational{PhysicalBasis(t, Basis::kComputational)};
    const auto c_hadamard{c_computational ^ 1U};
    const auto t_hadamard{t_computational ^ 1U};
    const auto c_hadamard_signs{Flags(c_hadamard)};
    const auto t_hadamard_signs{Flags(t_hadamard)};
    const bool c_hadamard_sign{GetBit(c_hadamard_signs, c)};
    const bool t_hadamard_sign{GetBit(t_hadamard_signs, t)};
    const auto hadamard_phase{
        PhysicalGenerator(c, c_hadamard)
            .XorWithPhase(PhysicalGenerator(t, t_hadamard), num_qubits_),
    };
    SetBit(c_hadamard_signs, c, Phase(hadamard_phase, c_hadamard_sign, t_hadamard_sign, 0));

    const auto c_computational_signs{Flags(c_computational)};
    const auto t_computational_signs{Flags(t_computational)};
    const bool c_computational_sign{GetBit(c_computational_signs, c)};
    const bool t_computational_sign{GetBit(t_computational_signs, t)};
    const auto computational_phase{
        PhysicalGenerator(t, t_computational)
            .XorWithPhase(PhysicalGenerator(c, c_computational), num_qubits_),
    };
    SetBit(t_computational_signs, t,
           Phase(computational_phase, t_computational_sign, c_computational_sign, 0));
}

std::size_t InverseTableau::MemoryUsage() const noexcept {
    return sizeof(*this) + StorageWords(capacity_) * sizeof(BitWord);
}

void InverseTableau::Print() const {
    for (std::size_t j{0}; j < num_qubits_; ++j) {
        std::cout << (j == 0 ? "\n    | " : "") << "x" << (j == num_qubits_ - 1 ? " | " : " ");
    }

    for (std::size_t j{0}; j < num_qubits_; ++j) {
        std::cout << "z" << (j == num_qubits_ - 1 ? " |" : " ");
    }

    std::cout << " r |";

    for (const auto basis : {Basis::kHadamard, Basis::kComputational}) {
        const char* const label{basis == Basis::kHadamard ? "X" : "Z"};
        for (QubitIdx i{0}; i < num_qubits_; ++i) {
            std::cout << "\n" << label << i << ": | ";
            const auto row{kPauliComponents * (PhysicalBasis(i, basis) * capacity_ + i)};
            for (std::size_t axis{0}; axis < kPauliComponents; ++axis) {
                const auto bits{RowWords(row + axis)};
                for (QubitIdx j{0}; j < num_qubits_; ++j) {
                    std::cout << ((bits[j / kBitsPerWord] >> (j % kBitsPerWord)) & 1U)
                              << (j == num_qubits_ - 1 ? " | " : " ");
                }
            }

            std::cout << Sign(i, basis) << " |";
        }
    }

    std::cout << '\n';
}

}  // namespace stabaal::detail
