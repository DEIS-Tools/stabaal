#include "stabaal/detail/inverse_tableau.hpp"

#include <algorithm>
#include <array>
#include <cassert>
#include <cstddef>
#include <cstdint>
#include <initializer_list>
#include <iostream>
#include <limits>
#include <memory>
#include <random>
#include <span>
#include <stdexcept>
#include <utility>
#include <vector>

#include "stabaal/detail/bit_word.hpp"
#include "stabaal/detail/pauli_string.hpp"
#include "stabaal/types.hpp"

namespace stabaal::detail {

namespace {

constexpr std::size_t kBasisCount{2};
constexpr std::size_t kPauliComponents{2};
constexpr std::size_t kSwapFlags{kBasisCount};
constexpr std::size_t kFlagArrays{kBasisCount + 1};
constexpr std::size_t kStackTargetWords{16};
constexpr std::uint64_t kGrowthRatePercent{25};

[[nodiscard]] std::mt19937& RandomEngine() {
    static thread_local std::mt19937 engine{std::random_device{}()};
    return engine;
}

[[nodiscard]] bool Phase(const std::int64_t exponent, const bool a_sign, const bool b_sign,
                         const std::uint8_t phase_rule) noexcept {
    assert(phase_rule < 4);
    const auto phase{static_cast<std::uint64_t>(exponent + phase_rule) & 3U};
    assert(phase % 2 == 0);
    return (phase == 2) != (a_sign != b_sign);
}

}  // namespace

std::size_t InverseTableau::StorageWords(const QubitIdx capacity) noexcept {
    return (kPauliComponents * kBasisCount * capacity + kFlagArrays) * WordCount(capacity);
}

std::span<std::uint64_t> InverseTableau::RowWords(const std::size_t row) noexcept {
    const auto stride{WordCount(capacity_)};
    return std::span{words_.get(), StorageWords(capacity_)}.subspan(row * stride, stride);
}

std::span<const std::uint64_t> InverseTableau::RowWords(const std::size_t row) const noexcept {
    const auto stride{WordCount(capacity_)};
    return std::span<const std::uint64_t>{words_.get(), StorageWords(capacity_)}.subspan(
        row * stride, stride);
}

std::span<std::uint64_t> InverseTableau::Flags(const std::size_t index) noexcept {
    return RowWords(kPauliComponents * kBasisCount * capacity_ + index);
}

std::span<const std::uint64_t> InverseTableau::Flags(const std::size_t index) const noexcept {
    return RowWords(kPauliComponents * kBasisCount * capacity_ + index);
}

std::size_t InverseTableau::PhysicalBasis(const QubitIdx q, const Basis basis) const noexcept {
    const auto swapped{static_cast<std::size_t>(GetBit(Flags(kSwapFlags), q))};
    return static_cast<std::size_t>(basis) ^ swapped;
}

PauliStringView InverseTableau::Generator(const QubitIdx q, const Basis basis) noexcept {
    const auto row{kPauliComponents * (PhysicalBasis(q, basis) * capacity_ + q)};
    return {RowWords(row), RowWords(row + 1)};
}

bool InverseTableau::Sign(const QubitIdx q, const Basis basis) const noexcept {
    return GetBit(Flags(PhysicalBasis(q, basis)), q);
}

void InverseTableau::FlipSign(const QubitIdx q, const Basis basis) noexcept {
    FlipBit(Flags(PhysicalBasis(q, basis)), q);
}

void InverseTableau::SetSign(const QubitIdx q, const Basis basis, const bool value) noexcept {
    SetBit(Flags(PhysicalBasis(q, basis)), q, value);
}

InverseTableau::InverseTableau(const QubitIdx num_qubits)
    : words_{
          // NOLINTNEXTLINE(cppcoreguidelines-avoid-c-arrays)
          num_qubits == 0 ? nullptr : std::make_unique<std::uint64_t[]>(StorageWords(num_qubits)),
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
              : std::make_unique_for_overwrite<std::uint64_t[]>(StorageWords(other.num_qubits_)),
      },
      num_qubits_{other.num_qubits_},
      capacity_{other.num_qubits_} {
    CopyRows(other);
}

InverseTableau& InverseTableau::operator=(const InverseTableau& other) {
    if (this != &other) *this = InverseTableau{other};
    return *this;
}

void InverseTableau::Add() {
    if (num_qubits_ == std::numeric_limits<QubitIdx>::max()) {
        throw std::length_error{"Too many qubits"};
    }
    if (num_qubits_ == capacity_) {
        const auto growth{std::max(capacity_ * kGrowthRatePercent / 100, std::uint64_t{1})};
        const auto extra{
            static_cast<QubitIdx>(
                std::min<std::uint64_t>(growth, std::numeric_limits<QubitIdx>::max() - capacity_)),
        };
        InverseTableau grown{0};
        grown.capacity_ = capacity_ + extra;
        // NOLINTNEXTLINE(cppcoreguidelines-avoid-c-arrays)
        grown.words_ = std::make_unique<std::uint64_t[]>(StorageWords(grown.capacity_));
        grown.CopyRows(*this);
        grown.num_qubits_ = num_qubits_;
        *this = std::move(grown);
    }
    const auto last{num_qubits_++};
    Generator(last, Basis::kComputational).SetZ(last);
    Generator(last, Basis::kHadamard).SetX(last);
}

void InverseTableau::CleanUp(const QubitIdx pivot, const std::span<const std::uint64_t> targets,
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
    std::array<std::uint64_t, kStackTargetWords> local_targets{};
    std::vector<std::uint64_t> heap_targets(num_words > local_targets.size() ? num_words : 0);
    const std::span<std::uint64_t> targets{
        heap_targets.empty() ? local_targets.data() : heap_targets.data(), num_words};
    const bool apply_s_dagger{gen.MeasurementTargets(pivot, targets)};
    if (apply_s_dagger ||
        std::any_of(targets.begin(), targets.end(), [](const auto word) { return word != 0; })) {
        CleanUp(pivot, targets, apply_s_dagger);
    }
    const bool p{std::bernoulli_distribution{0.5}(RandomEngine())};
    Collapse(pivot, p);
    return Sign(q, b);
}

void InverseTableau::CheckQubit(const QubitIdx q) const {
    if (q >= num_qubits_) [[unlikely]] {
        throw std::out_of_range{"Qubit index out of range"};
    }
}

bool InverseTableau::Measure(const QubitIdx q, const Basis basis) {
    CheckQubit(q);
    return Generator(q, basis).AnyX(num_qubits_) ? RandomMeasurement(q, basis) : Sign(q, basis);
}

bool InverseTableau::MeasureX(const QubitIdx q) { return Measure(q, Basis::kHadamard); }

bool InverseTableau::MeasureZ(const QubitIdx q) { return Measure(q, Basis::kComputational); }

void InverseTableau::ApplyH(const QubitIdx q) {
    CheckQubit(q);

    FlipBit(Flags(kSwapFlags), q);
}

void InverseTableau::ApplyS(const QubitIdx q) {
    CheckQubit(q);

    const auto exponent{
        Generator(q, Basis::kHadamard)
            .XorWithPhase(Generator(q, Basis::kComputational), num_qubits_),
    };
    SetSign(q, Basis::kHadamard,
            Phase(exponent, Sign(q, Basis::kHadamard), Sign(q, Basis::kComputational), 3));
}

void InverseTableau::ApplyX(const QubitIdx q) {
    CheckQubit(q);

    FlipSign(q, Basis::kComputational);
}

void InverseTableau::ApplyZ(const QubitIdx q) {
    CheckQubit(q);

    FlipSign(q, Basis::kHadamard);
}

void InverseTableau::ApplyCx(const QubitIdx c, const QubitIdx t) {
    CheckQubit(c);
    CheckQubit(t);

    const auto hadamard_phase{
        Generator(c, Basis::kHadamard).XorWithPhase(Generator(t, Basis::kHadamard), num_qubits_),
    };
    SetSign(c, Basis::kHadamard,
            Phase(hadamard_phase, Sign(c, Basis::kHadamard), Sign(t, Basis::kHadamard), 0));

    const auto computational_phase{
        Generator(t, Basis::kComputational)
            .XorWithPhase(Generator(c, Basis::kComputational), num_qubits_),
    };
    SetSign(t, Basis::kComputational,
            Phase(computational_phase, Sign(t, Basis::kComputational),
                  Sign(c, Basis::kComputational), 0));
}

std::size_t InverseTableau::MemoryUsage() const noexcept {
    return sizeof(*this) + StorageWords(capacity_) * sizeof(std::uint64_t);
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
