#include "stabaal/backends/inverse_tableau.hpp"

#include <algorithm>
#include <array>
#include <bit>
#include <cassert>
#include <cstddef>
#include <cstdint>
#include <initializer_list>
#include <iostream>
#include <limits>
#include <memory>
#include <new>
#include <optional>
#include <span>
#include <stdexcept>
#include <utility>
#include <vector>

#include "stabaal/detail/bit_ops.hpp"
#include "stabaal/detail/pauli_string_view.hpp"
#include "stabaal/random.hpp"
#include "stabaal/types.hpp"

namespace stabaal::detail {

namespace {

constexpr std::size_t kBasisCount{2};
constexpr std::size_t kPauliComponents{2};
constexpr std::size_t kSwapFlags{kBasisCount};
constexpr std::size_t kTombstoneFlags{3};
constexpr std::size_t kPurificationFlags{4};
constexpr std::size_t kFlagArrays{kBasisCount + 3};
constexpr std::size_t kStackTargetWords{16};
constexpr QubitSize kGrowthDivisor{4};

// The complete QubitSize range needs up to twice its bit width for storage sizes.
static_assert(std::numeric_limits<std::size_t>::digits >=
              2 * std::numeric_limits<QubitSize>::digits);

[[nodiscard]] bool Phase(const PhaseExponent exponent, const bool a_sign, const bool b_sign,
                         const PhaseExponent phase_rule) noexcept {
    assert(phase_rule >= 0 && phase_rule < 4);
    const auto phase{((exponent + phase_rule) % 4 + 4) % 4};
    assert(phase % 2 == 0);
    return (phase == 2) != (a_sign != b_sign);
}

}  // namespace

std::size_t InverseTableau::NumWords(const QubitSize capacity) noexcept {
    return (std::size_t{4} * capacity + kFlagArrays) * WordCount(capacity);
}

BitWordSpan InverseTableau::RowWords(const std::size_t row) noexcept {
    const auto stride{WordCount(capacity_)};
    return BitWordSpan{words_.get(), NumWords(capacity_)}.subspan(row * stride, stride);
}

ConstBitWordSpan InverseTableau::RowWords(const std::size_t row) const noexcept {
    const auto stride{WordCount(capacity_)};
    return ConstBitWordSpan{words_.get(), NumWords(capacity_)}.subspan(row * stride, stride);
}

BitWordSpan InverseTableau::Flags(const std::size_t index) noexcept {
    return RowWords(std::size_t{4} * capacity_ + index);
}

ConstBitWordSpan InverseTableau::Flags(const std::size_t index) const noexcept {
    return RowWords(std::size_t{4} * capacity_ + index);
}

std::size_t InverseTableau::PhysicalBasis(const QubitSize generator_index,
                                          const Basis basis) const noexcept {
    const auto swapped{static_cast<std::size_t>(GetBit(Flags(kSwapFlags), generator_index))};
    return static_cast<std::size_t>(basis) ^ swapped;
}

PauliStringView InverseTableau::Generator(const QubitSize generator_index,
                                          const Basis basis) noexcept {
    return PhysicalGenerator(generator_index, PhysicalBasis(generator_index, basis));
}

PauliStringView InverseTableau::PhysicalGenerator(const QubitSize generator_index,
                                                  const std::size_t basis) noexcept {
    const auto row{kPauliComponents * (basis * capacity_ + generator_index)};
    return {RowWords(row), RowWords(row + 1)};
}

bool InverseTableau::Sign(const QubitSize generator_index, const Basis basis) const noexcept {
    return GetBit(Flags(PhysicalBasis(generator_index, basis)), generator_index);
}

void InverseTableau::FlipSign(const QubitSize generator_index, const Basis basis) noexcept {
    FlipBit(Flags(PhysicalBasis(generator_index, basis)), generator_index);
}

InverseTableau::InverseTableau(const QubitSize num_qubits)
    : InverseTableau(num_qubits, std::nullopt) {}

InverseTableau::InverseTableau(const QubitSize num_qubits, const std::optional<std::uint64_t> seed)
    : words_{
          num_qubits == 0 ? nullptr : std::make_unique<BitWord[]>(NumWords(num_qubits)),
      },
      num_qubits_{num_qubits},
      num_slots_{num_qubits},
      capacity_{num_qubits},
      measurement_engine_{seed ? Sfc64{*seed} : Sfc64{}} {
    for (QubitSize i{0}; i < num_qubits; ++i) {
        Generator(i, Basis::kComputational).SetZ(i);
        Generator(i, Basis::kHadamard).SetX(i);
    }
}

void InverseTableau::CopyRows(const InverseTableau& other) noexcept {
    const auto used_words{WordCount(other.num_slots_)};
    if (used_words == 0) {
        return;
    }

    if (capacity_ == other.capacity_) {
        std::copy_n(other.words_.get(), NumWords(capacity_), words_.get());
        return;
    }

    for (std::size_t basis{0}; basis < kBasisCount; ++basis) {
        for (QubitSize generator_index{0}; generator_index < other.num_slots_; ++generator_index) {
            for (std::size_t axis{0}; axis < kPauliComponents; ++axis) {
                const auto source{
                    other.RowWords(kPauliComponents * (basis * other.capacity_ + generator_index) +
                                   axis),
                };
                const auto target{
                    RowWords(kPauliComponents * (basis * capacity_ + generator_index) + axis),
                };
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
          other.num_slots_ == 0
              ? nullptr
              : std::make_unique_for_overwrite<BitWord[]>(NumWords(other.num_slots_)),
      },
      num_qubits_{other.num_qubits_},
      num_slots_{other.num_slots_},
      capacity_{other.num_slots_},
      first_free_slot_{other.first_free_slot_},
      num_free_slots_{other.num_free_slots_},
      num_ids_{other.num_ids_},
      id_capacity_{other.num_ids_},
      first_free_id_{other.first_free_id_},
      generation_epoch_{other.generation_epoch_},
      measurement_engine_{other.measurement_engine_},
      measurement_random_{other.measurement_random_} {
    CopyRows(other);
    if (other.generations_) {
        const auto count{id_capacity_ != 0 ? std::size_t{2} * id_capacity_ : capacity_};
        generations_ = std::make_unique_for_overwrite<std::uint32_t[]>(count);
        std::copy_n(other.generations_.get(), count, generations_.get());
    }
}

InverseTableau& InverseTableau::operator=(const InverseTableau& other) {
    if (this == &other) {
        return *this;
    }

    if (words_ && capacity_ == other.num_slots_) {
        const auto count{other.id_capacity_ != 0 ? std::size_t{2} * other.num_ids_ : capacity_};
        auto generations{
            other.generations_ ? std::make_unique_for_overwrite<std::uint32_t[]>(count) : nullptr,
        };
        if (generations) {
            std::copy_n(other.generations_.get(), count, generations.get());
        }

        CopyRows(other);
        num_qubits_ = other.num_qubits_;
        num_slots_ = other.num_slots_;
        first_free_slot_ = other.first_free_slot_;
        num_free_slots_ = other.num_free_slots_;
        num_ids_ = other.num_ids_;
        id_capacity_ = other.num_ids_;
        first_free_id_ = other.first_free_id_;
        generations_ = std::move(generations);
        generation_epoch_ = other.generation_epoch_;
        measurement_engine_ = other.measurement_engine_;
        measurement_random_ = other.measurement_random_;
    } else {
        *this = InverseTableau{other};
    }

    return *this;
}

QubitSize InverseTableau::NumQubits() const noexcept { return num_qubits_; }

std::uint32_t InverseTableau::Generation(const QubitSize id_slot) const noexcept {
    if (!generations_) {
        return static_cast<std::uint32_t>(generation_epoch_);
    }

    const auto count{id_capacity_ != 0 ? std::size_t{2} * id_capacity_ : capacity_};
    const auto offset{id_capacity_ != 0 ? std::size_t{2} * id_slot + 1 : id_slot};
    return std::span<const std::uint32_t>{generations_.get(), count}[offset];
}

std::uint32_t& InverseTableau::Generation(const QubitSize id_slot) noexcept {
    assert(generations_);
    const auto count{id_capacity_ != 0 ? std::size_t{2} * id_capacity_ : capacity_};
    const auto offset{id_capacity_ != 0 ? std::size_t{2} * id_slot + 1 : id_slot};
    return std::span{generations_.get(), count}[offset];
}

QubitSize InverseTableau::IdIndex(const QubitSize id_slot) const noexcept {
    assert(id_capacity_ != 0);
    return std::span<const std::uint32_t>{generations_.get(),
                                          std::size_t{2} * id_capacity_}[std::size_t{2} * id_slot];
}

QubitSize& InverseTableau::IdIndex(const QubitSize id_slot) noexcept {
    assert(id_capacity_ != 0);
    return std::span{generations_.get(), std::size_t{2} * id_capacity_}[std::size_t{2} * id_slot];
}

QubitSize InverseTableau::ResolveQubit(const QubitId id) const {
    if (!generations_ && generation_epoch_ == 0) {
        if (id < num_slots_ &&
            (num_qubits_ == num_slots_ ||
             (!GetBit(Flags(kTombstoneFlags), id) && !GetBit(Flags(kPurificationFlags), id)))) {
            return static_cast<QubitSize>(id);
        }

        throw std::out_of_range{"Qubit id is not active"};
    }

    const auto slot{static_cast<QubitSize>(id)};
    const auto generation{static_cast<std::uint32_t>(id >> 32U)};
    const auto count{id_capacity_ != 0 ? num_ids_ : num_slots_};
    if (slot < count && (generation & 1U) == 0 && Generation(slot) == generation) {
        if (id_capacity_ != 0) {
            return IdIndex(slot);
        }

        if (num_qubits_ == num_slots_ ||
            (!GetBit(Flags(kTombstoneFlags), slot) && !GetBit(Flags(kPurificationFlags), slot))) {
            return slot;
        }
    }

    throw std::out_of_range{"Qubit id is not active"};
}

void InverseTableau::MapIds() {
    if (id_capacity_ != 0) {
        return;
    }

    assert(num_slots_ != 0);
    auto mapped{std::make_unique_for_overwrite<std::uint32_t[]>(std::size_t{2} * num_slots_)};
    const std::span entries{mapped.get(), std::size_t{2} * num_slots_};
    QubitSize free_id{kNoSlot};
    for (QubitSize slot{0}; slot < num_slots_; ++slot) {
        auto generation{std::as_const(*this).Generation(slot)};
        if (GetBit(Flags(kTombstoneFlags), slot) || GetBit(Flags(kPurificationFlags), slot)) {
            generation |= 1U;
            entries[std::size_t{2} * slot] = free_id;
            if (generation != std::numeric_limits<std::uint32_t>::max()) {
                free_id = slot;
            }
        } else {
            entries[std::size_t{2} * slot] = slot;
        }

        entries[std::size_t{2} * slot + 1] = generation;
    }

    generations_ = std::move(mapped);
    num_ids_ = num_slots_;
    id_capacity_ = num_slots_;
    first_free_id_ = free_id;
}

void InverseTableau::ReserveId() {
    if (num_ids_ == kNoSlot) {
        throw std::length_error{"Qubit id slots exhausted"};
    }

    if (num_ids_ < id_capacity_) {
        return;
    }

    const auto growth{
        std::min(std::max(id_capacity_ / kGrowthDivisor, QubitSize{1}), kNoSlot - id_capacity_),
    };
    const auto capacity{id_capacity_ + growth};
    auto entries{std::make_unique_for_overwrite<std::uint32_t[]>(std::size_t{2} * capacity)};
    std::copy_n(generations_.get(), std::size_t{2} * num_ids_, entries.get());
    generations_ = std::move(entries);
    id_capacity_ = capacity;
}

void InverseTableau::AdvanceFreeSlot() noexcept {
    assert(num_free_slots_ != 0);
    const auto num_words{WordCount(num_slots_)};
    const auto tombstones{Flags(kTombstoneFlags)};
    for (auto word{first_free_slot_ / kBitsPerWord}; word < num_words; ++word) {
        const auto free{tombstones[word]};
        if (free != 0) {
            first_free_slot_ = static_cast<QubitSize>(word * kBitsPerWord + std::countr_zero(free));
            return;
        }
    }

    assert(false && "Reusable slot count disagrees with tombstones");
}

QubitSize InverseTableau::AppendSlot() {
    if (num_slots_ == std::numeric_limits<QubitSize>::max()) {
        throw std::length_error{"Too many stored qubit slots"};
    }

    if (num_slots_ == capacity_) {
        const auto growth{std::max(capacity_ / kGrowthDivisor, QubitSize{1})};
        const auto extra{std::min(growth, std::numeric_limits<QubitSize>::max() - capacity_)};
        InverseTableau grown{0};
        grown.capacity_ = capacity_ + extra;
        grown.words_ = std::make_unique<BitWord[]>(NumWords(grown.capacity_));
        grown.CopyRows(*this);
        if (generations_ && id_capacity_ == 0) {
            grown.generations_ = std::make_unique<std::uint32_t[]>(grown.capacity_);
            std::copy_n(generations_.get(), num_slots_, grown.generations_.get());
        }

        words_ = std::move(grown.words_);
        if (id_capacity_ == 0) {
            generations_ = std::move(grown.generations_);
        }

        capacity_ = grown.capacity_;
    }

    const auto index{num_slots_++};
    Generator(index, Basis::kComputational).SetZ(index);
    Generator(index, Basis::kHadamard).SetX(index);
    return index;
}

QubitId InverseTableau::Add() {
    if (num_qubits_ == std::numeric_limits<QubitSize>::max()) {
        throw std::length_error{"Too many qubits"};
    }

    if (id_capacity_ != 0 && first_free_id_ == kNoSlot) {
        ReserveId();
    }

    if (generation_epoch_ > std::numeric_limits<std::uint32_t>::max()) {
        throw std::length_error{"Qubit generations exhausted"};
    }

    const auto slot{first_free_slot_};
    if (slot != kNoSlot && !generations_) {
        generations_ = std::make_unique_for_overwrite<std::uint32_t[]>(capacity_);
        const std::span generations{generations_.get(), capacity_};
        const auto tombstones{Flags(kTombstoneFlags)};
        const auto purification{Flags(kPurificationFlags)};
        for (QubitSize q{0}; q < capacity_; ++q) {
            generations[q] =
                static_cast<std::uint32_t>(generation_epoch_) +
                static_cast<std::uint32_t>(GetBit(tombstones, q) || GetBit(purification, q));
        }
    }

    const auto index{slot != kNoSlot ? slot : AppendSlot()};
    QubitSize id_slot{index};
    auto generation{static_cast<std::uint32_t>(generation_epoch_)};
    if (id_capacity_ != 0) {
        if (first_free_id_ != kNoSlot) {
            id_slot = first_free_id_;
            first_free_id_ = IdIndex(id_slot);
            generation = Generation(id_slot) + 1;
        } else {
            id_slot = num_ids_++;
        }

        IdIndex(id_slot) = index;
        Generation(id_slot) = generation;
    } else if (generations_) {
        if (slot != kNoSlot) {
            generation = Generation(slot) + 1;
        }

        Generation(index) = generation;
    }
    generation_epoch_ = std::max(generation_epoch_, std::uint64_t{generation});
    if (slot != kNoSlot) {
        SetBit(Flags(kTombstoneFlags), slot, false);
        if (--num_free_slots_ == 0) {
            first_free_slot_ = kNoSlot;
        } else {
            AdvanceFreeSlot();
        }
    }

    ++num_qubits_;
    return (QubitId{generation} << 32U) | id_slot;
}

void InverseTableau::ResetProductSlot(const QubitSize q) noexcept {
    if (!IsDeterministicMeasure(q, Basis::kComputational)) {
        if (!IsDeterministicMeasure(q, Basis::kHadamard)) {
            ApplyPhase(q, 1);
        }
        FlipBit(Flags(kSwapFlags), q);
    }

    if (Sign(q, Basis::kComputational)) {
        FlipSign(q, Basis::kComputational);
    }
}

void InverseTableau::Discard(const QubitId id) {
    const auto q{ResolveQubit(id)};
    const auto id_slot{static_cast<QubitSize>(id)};
    const bool pure{!RequiresPurification(q)};
    if (id_capacity_ == 0 && generations_ &&
        Generation(id_slot) == std::numeric_limits<std::uint32_t>::max() - 1) {
        MapIds();
    }

    if (pure) {
        ResetProductSlot(q);
    }

    if (generations_) {
        ++Generation(id_slot);
    }

    if (id_capacity_ != 0) {
        IdIndex(id_slot) = first_free_id_;
        if (Generation(id_slot) != std::numeric_limits<std::uint32_t>::max()) {
            first_free_id_ = id_slot;
        }
    }

    SetBit(Flags(kTombstoneFlags), q, pure);
    SetBit(Flags(kPurificationFlags), q, !pure);
    if (pure) {
        ++num_free_slots_;
        first_free_slot_ = std::min(first_free_slot_, q);
    }

    --num_qubits_;
    try {
        Contract();
    } catch (const std::bad_alloc&) {
        // Discard succeeded. Automatic contraction is best-effort when allocation fails.
        return;
    }
}

void InverseTableau::RemoveProductSlot(const QubitSize q, const BitWordSpan targets,
                                       const std::span<QubitSize> origins) noexcept {
    const auto num_words{WordCount(num_slots_)};
    const auto used{targets.first(num_words)};
    const auto z_row{kPauliComponents * (PhysicalBasis(q, Basis::kComputational) * capacity_ + q)};
    std::copy_n(RowWords(z_row + 1).begin(), num_words, used.begin());
    QubitSize pivot{0};
    while (!GetBit(used, pivot)) {
        ++pivot;
    }
    SetBit(used, pivot, false);
    if (std::any_of(used.begin(), used.end(), [](const auto word) { return word != 0; })) {
        for (std::size_t basis{0}; basis < kBasisCount; ++basis) {
            const auto signs{Flags(basis)};
            for (QubitSize gen{0}; gen < num_slots_; ++gen) {
                const auto row{kPauliComponents * (basis * capacity_ + gen)};
                PauliStringView reversed{RowWords(row + 1), RowWords(row)};
                FlipBit(signs, gen, reversed.CleanUp(pivot, used, false));
            }
        }
    }

    const auto x{Generator(q, Basis::kHadamard)};
    const bool phase{x.MeasurementTargets(pivot, used)};
    if (phase || std::any_of(used.begin(), used.end(), [](const auto word) { return word != 0; })) {
        CleanUp(pivot, used, phase);
    }

    const auto last{num_slots_ - 1};
    for (std::size_t basis{0}; basis < kBasisCount; ++basis) {
        for (QubitSize gen{0}; gen < num_slots_; ++gen) {
            const auto row{kPauliComponents * (basis * capacity_ + gen)};
            for (std::size_t axis{0}; axis < kPauliComponents; ++axis) {
                const auto bits{RowWords(row + axis)};
                SetBit(bits, pivot, GetBit(bits, last));
                SetBit(bits, last, false);
            }
        }

        if (q != last) {
            for (std::size_t axis{0}; axis < kPauliComponents; ++axis) {
                const auto source{RowWords(kPauliComponents * (basis * capacity_ + last) + axis)};
                const auto target{RowWords(kPauliComponents * (basis * capacity_ + q) + axis)};
                std::copy(source.begin(), source.end(), target.begin());
            }
        }
    }

    for (std::size_t flag{0}; flag < kFlagArrays; ++flag) {
        SetBit(Flags(flag), q, GetBit(Flags(flag), last));
        SetBit(Flags(flag), last, false);
    }

    origins[q] = origins[last];
    --num_slots_;
}

void InverseTableau::Contract() {
    if (num_qubits_ == 0 && generation_epoch_ < std::numeric_limits<std::uint32_t>::max() - 1) {
        words_.reset();
        generations_.reset();
        num_slots_ = 0;
        capacity_ = 0;
        num_free_slots_ = 0;
        first_free_slot_ = kNoSlot;
        num_ids_ = 0;
        id_capacity_ = 0;
        first_free_id_ = kNoSlot;
        generation_epoch_ += 2;
        return;
    }

    const auto retained{num_slots_ - num_free_slots_};
    if (num_free_slots_ == 0 || retained > capacity_ / 2) {
        return;
    }

    InverseTableau compacted{0};
    compacted.capacity_ = retained;

    compacted.words_ = std::make_unique<BitWord[]>(NumWords(retained));
    std::vector<BitWord> targets(WordCount(num_slots_));
    std::vector<QubitSize> origins(num_slots_);
    std::vector<QubitSize> relocation(num_slots_, kNoSlot);
    for (QubitSize q{0}; q < num_slots_; ++q) {
        origins[q] = q;
    }

    MapIds();
    QubitSize q{0};
    while (q < num_slots_) {
        if (GetBit(Flags(kTombstoneFlags), q)) {
            RemoveProductSlot(q, targets, origins);
        } else {
            ++q;
        }
    }

    for (QubitSize index{0}; index < num_slots_; ++index) {
        relocation[origins[index]] = index;
    }

    for (QubitSize id_slot{0}; id_slot < num_ids_; ++id_slot) {
        if ((Generation(id_slot) & 1U) == 0) {
            IdIndex(id_slot) = relocation[IdIndex(id_slot)];
        }
    }

    compacted.CopyRows(*this);
    words_ = std::move(compacted.words_);
    capacity_ = retained;
    num_free_slots_ = 0;
    first_free_slot_ = kNoSlot;
}

void InverseTableau::CleanUp(const QubitSize pivot, const ConstBitWordSpan targets,
                             const bool apply_s_dagger) noexcept {
    for (std::size_t basis{0}; basis < kBasisCount; ++basis) {
        const auto signs{Flags(basis)};
        for (QubitSize generator_index{0}; generator_index < num_slots_; ++generator_index) {
            const auto row{kPauliComponents * (basis * capacity_ + generator_index)};
            PauliStringView gen{RowWords(row), RowWords(row + 1)};
            const bool flip{gen.CleanUp(pivot, targets, apply_s_dagger)};
            FlipBit(signs, generator_index, flip);
        }
    }
}

void InverseTableau::Collapse(const QubitSize pivot, const bool p) noexcept {
    for (std::size_t basis{0}; basis < kBasisCount; ++basis) {
        const auto signs{Flags(basis)};
        for (QubitSize generator_index{0}; generator_index < num_slots_; ++generator_index) {
            const auto row{kPauliComponents * (basis * capacity_ + generator_index)};
            PauliStringView gen{RowWords(row), RowWords(row + 1)};
            gen.SwapXZ(pivot);
            const bool flip{gen.GetZ(pivot) && (p != gen.GetX(pivot))};
            FlipBit(signs, generator_index, flip);
        }
    }
}

bool InverseTableau::RandomMeasurement(const QubitSize q, const Basis b,
                                       std::optional<bool> measurement_opt) {
    const auto gen{Generator(q, b)};
    QubitSize pivot{0};

    for (QubitSize i{0}; i < num_slots_; ++i) {
        if (gen.GetX(i)) {
            pivot = i;
            break;
        }
    }

    const auto num_words{WordCount(num_slots_)};
    std::array<BitWord, kStackTargetWords> local_targets{};
    std::vector<BitWord> heap_targets(num_words > local_targets.size() ? num_words : 0);
    const BitWordSpan targets{heap_targets.empty() ? local_targets.data() : heap_targets.data(),
                              num_words};
    const bool apply_s_dagger{gen.MeasurementTargets(pivot, targets)};
    if (apply_s_dagger ||
        std::any_of(targets.begin(), targets.end(), [](const auto word) { return word != 0; })) {
        CleanUp(pivot, targets, apply_s_dagger);
    }

    bool p{false};
    if (measurement_opt) {
        p = measurement_opt.value();
    } else {
        p = measurement_random_(measurement_engine_);
    }

    Collapse(pivot, p);
    return Sign(q, b);
}

bool InverseTableau::Measure(const QubitSize q, const Basis basis) {
    assert(q < num_slots_);
    return IsDeterministicMeasure(q, basis) ? Sign(q, basis) : RandomMeasurement(q, basis);
}

bool InverseTableau::IsDeterministicMeasure(const QubitSize q, const Basis basis) const noexcept {
    assert(q < num_slots_);
    const auto row{kPauliComponents * (PhysicalBasis(q, basis) * capacity_ + q)};
    const auto words{RowWords(row).first(WordCount(num_slots_))};
    return std::none_of(words.begin(), words.end(), [](const auto word) { return word != 0; });
}

bool InverseTableau::IsDeterministicMeasureZ(const QubitId id) const {
    return IsDeterministicMeasure(ResolveQubit(id), Basis::kComputational);
}

bool InverseTableau::IsDeterministicMeasureX(const QubitId id) const {
    return IsDeterministicMeasure(ResolveQubit(id), Basis::kHadamard);
}

bool InverseTableau::RequiresPurification(const QubitSize q) const noexcept {
    const auto first_x{RowWords(kPauliComponents * q)};
    const auto second_x{RowWords(kPauliComponents * (std::size_t{capacity_} + q))};
    BitWord any_first_x{0};
    BitWord any_second_x{0};
    BitWord any_y_x{0};

    const auto num_words{WordCount(num_slots_)};
    for (std::size_t i{0}; i < num_words; ++i) {
        any_first_x |= first_x[i];
        any_second_x |= second_x[i];
        any_y_x |= first_x[i] ^ second_x[i];
        if (any_first_x != 0 && any_second_x != 0 && any_y_x != 0) {
            return true;
        }
    }

    return false;
}

bool InverseTableau::MeasureX(const QubitId id) {
    return Measure(ResolveQubit(id), Basis::kHadamard);
}

bool InverseTableau::MeasureZ(const QubitId id) {
    return Measure(ResolveQubit(id), Basis::kComputational);
}

bool InverseTableau::MeasureAndDiscardX(const QubitId id) {
    const auto outcome = MeasureX(id);
    Discard(id);
    return outcome;
}

bool InverseTableau::MeasureAndDiscardZ(const QubitId id) {
    const auto outcome = MeasureZ(id);
    Discard(id);
    return outcome;
}

bool InverseTableau::ForceMeasure(const QubitSize q, bool measurement, Basis basis) {
    return RandomMeasurement(q, basis, measurement);
}

bool InverseTableau::ForceMeasureZ(const QubitId id, const bool measurement) {
    return ForceMeasure(ResolveQubit(id), measurement, Basis::kComputational);
}

bool InverseTableau::ForceMeasureX(const QubitId id, const bool measurement) {
    return ForceMeasure(ResolveQubit(id), measurement, Basis::kHadamard);
}

void InverseTableau::ApplyH(const QubitId id) {
    const auto q{ResolveQubit(id)};
    assert(q < num_slots_);

    FlipBit(Flags(kSwapFlags), q);
}

void InverseTableau::ApplyPhase(const QubitSize q, const PhaseExponent phase_rule) noexcept {
    const auto computational{PhysicalBasis(q, Basis::kComputational)};
    const auto hadamard{computational ^ 1U};
    const auto computational_signs{Flags(computational)};
    const auto hadamard_signs{Flags(hadamard)};
    const bool computational_sign{GetBit(computational_signs, q)};
    const bool hadamard_sign{GetBit(hadamard_signs, q)};
    const auto exponent{
        PhysicalGenerator(q, hadamard)
            .XorWithPhase(PhysicalGenerator(q, computational), num_slots_),
    };
    SetBit(hadamard_signs, q, Phase(exponent, hadamard_sign, computational_sign, phase_rule));
}

void InverseTableau::ApplyS(const QubitId id) { ApplyPhase(ResolveQubit(id), 3); }

void InverseTableau::ApplySdg(const QubitId id) { ApplyPhase(ResolveQubit(id), 1); }

void InverseTableau::ApplySx(const QubitId id) {
    ApplyH(id);
    ApplyS(id);
    ApplyH(id);
}

void InverseTableau::ApplyI(const QubitId id) { static_cast<void>(ResolveQubit(id)); }

void InverseTableau::ApplyX(const QubitId id) {
    const auto q{ResolveQubit(id)};
    assert(q < num_slots_);

    FlipSign(q, Basis::kComputational);
}

void InverseTableau::ApplyZ(const QubitId id) {
    const auto q{ResolveQubit(id)};
    assert(q < num_slots_);

    FlipSign(q, Basis::kHadamard);
}

void InverseTableau::ApplyY(const QubitId id) {
    const auto q{ResolveQubit(id)};
    // We dont need to check swap flags so no need to look those up.
    // Equivalent to ApplyX(q); ApplyZ(q);
    FlipBit(Flags(0), q);
    FlipBit(Flags(1), q);
}

void InverseTableau::ApplyCx(const QubitId control, const QubitId target) {
    const auto control_slot{ResolveQubit(control)};
    const auto target_slot{ResolveQubit(target)};
    assert(control_slot < num_slots_);
    assert(target_slot < num_slots_);

    const auto control_computational{PhysicalBasis(control_slot, Basis::kComputational)};
    const auto target_computational{PhysicalBasis(target_slot, Basis::kComputational)};
    const auto control_hadamard{control_computational ^ 1U};
    const auto target_hadamard{target_computational ^ 1U};
    const auto control_hadamard_signs{Flags(control_hadamard)};
    const auto target_hadamard_signs{Flags(target_hadamard)};
    const bool control_hadamard_sign{GetBit(control_hadamard_signs, control_slot)};
    const bool target_hadamard_sign{GetBit(target_hadamard_signs, target_slot)};
    const auto hadamard_phase{
        PhysicalGenerator(control_slot, control_hadamard)
            .XorWithPhase(PhysicalGenerator(target_slot, target_hadamard), num_slots_),
    };
    SetBit(control_hadamard_signs, control_slot,
           Phase(hadamard_phase, control_hadamard_sign, target_hadamard_sign, 0));

    const auto control_computational_signs{Flags(control_computational)};
    const auto target_computational_signs{Flags(target_computational)};
    const bool control_computational_sign{GetBit(control_computational_signs, control_slot)};
    const bool target_computational_sign{GetBit(target_computational_signs, target_slot)};
    const auto computational_phase{
        PhysicalGenerator(target_slot, target_computational)
            .XorWithPhase(PhysicalGenerator(control_slot, control_computational), num_slots_),
    };
    SetBit(target_computational_signs, target_slot,
           Phase(computational_phase, target_computational_sign, control_computational_sign, 0));
}

void InverseTableau::ApplyCz(const QubitId control, const QubitId target) {
    static_cast<void>(ResolveQubit(control));
    static_cast<void>(ResolveQubit(target));
    ApplyH(target);
    ApplyCx(control, target);
    ApplyH(target);
}

void InverseTableau::ApplyCy(const QubitId control, const QubitId target) {
    static_cast<void>(ResolveQubit(control));
    static_cast<void>(ResolveQubit(target));
    ApplySdg(target);
    ApplyCx(control, target);
    ApplyS(target);
}

std::size_t InverseTableau::MemoryUsage() const noexcept {
    return sizeof(*this) + NumWords(capacity_) * sizeof(BitWord) +
           (generations_ ? (id_capacity_ != 0 ? std::size_t{2} * id_capacity_ : capacity_) *
                               sizeof(std::uint32_t)
                         : 0);
}

void InverseTableau::Print() const {
    for (std::size_t j{0}; j < num_slots_; ++j) {
        std::cout << (j == 0 ? "\n    | " : "") << "x" << (j == num_slots_ - 1 ? " | " : " ");
    }

    for (std::size_t j{0}; j < num_slots_; ++j) {
        std::cout << "z" << (j == num_slots_ - 1 ? " |" : " ");
    }

    std::cout << " r |";

    for (const auto basis : {Basis::kHadamard, Basis::kComputational}) {
        const char* const label{basis == Basis::kHadamard ? "X" : "Z"};
        for (QubitSize generator_index{0}; generator_index < num_slots_; ++generator_index) {
            std::cout << "\n" << label << generator_index;
            if (GetBit(Flags(kPurificationFlags), generator_index)) {
                std::cout << " [p]";
            } else if (GetBit(Flags(kTombstoneFlags), generator_index)) {
                std::cout << " [†]";
            } else {
                std::cout << " [a]";
            }

            std::cout << ": | ";
            const auto row{
                kPauliComponents *
                    (PhysicalBasis(generator_index, basis) * capacity_ + generator_index),
            };
            for (std::size_t axis{0}; axis < kPauliComponents; ++axis) {
                const auto bits{RowWords(row + axis)};
                for (QubitSize j{0}; j < num_slots_; ++j) {
                    std::cout << ((bits[j / kBitsPerWord] >> (j % kBitsPerWord)) & 1U)
                              << (j == num_slots_ - 1 ? " | " : " ");
                }
            }

            std::cout << Sign(generator_index, basis) << " |";
        }
    }

    std::cout << '\n';
}

InverseTableau::GateInformation InverseTableau::GetGateInformation(UnaryGateFunction gate) {
    if (gate == &InverseTableau::ApplyX) {
        return {"X gate", "Applies the X gate", 1};
    }
    if (gate == &InverseTableau::ApplyY) {
        return {"Y gate", "Applies the Y gate", 1};
    }
    if (gate == &InverseTableau::ApplyZ) {
        return {"Z gate", "Applies the Z gate", 1};
    }
    if (gate == &InverseTableau::ApplyH) {
        return {"Hadamard gate", "Applies the Hadamard gate", 1};
    }
    if (gate == &InverseTableau::ApplySx) {
        return {"Sx gate", "Applies the Square Root X gate", 1};
    }
    if (gate == &InverseTableau::ApplyS) {
        return {"S gate", "Applies the S gate", 1};
    }
    if (gate == &InverseTableau::ApplySdg) {
        return {"S-dagger gate", "Applies the S-Dagger gate", 1};
    }
    if (gate == &InverseTableau::ApplyI) {
        return {"Identity gate", "Applies the Identity gate", 1};
    };
    throw std::range_error("Function is not a supported gate function");
}

InverseTableau::GateInformation InverseTableau::GetGateInformation(BinaryGateFunction gate) {
    if (gate == &InverseTableau::ApplyCx) {
        return {"Control-X gate", "Applies the Control-X gate", 2};
    }
    if (gate == &InverseTableau::ApplyCy) {
        return {"Control-Y gate", "Applies the Control-Y gate", 2};
    }
    if (gate == &InverseTableau::ApplyCz) {
        return {"Control-Z gate", "Applies the Control-Z gate", 2};
    };
    throw std::range_error("Function is not a supported gate function");
}

std::array<InverseTableau::UnaryGateFunction, 8> InverseTableau::SupportedUnaryGates() {
    return std::array<InverseTableau::UnaryGateFunction, 8>{
        &InverseTableau::ApplyX,   &InverseTableau::ApplyY,  &InverseTableau::ApplyZ,
        &InverseTableau::ApplyH,   &InverseTableau::ApplySx, &InverseTableau::ApplyS,
        &InverseTableau::ApplySdg, &InverseTableau::ApplyI,
    };
}

std::array<InverseTableau::BinaryGateFunction, 3> InverseTableau::SupportedBinaryGates() {
    return std::array<InverseTableau::BinaryGateFunction, 3>{
        &InverseTableau::ApplyCx,
        &InverseTableau::ApplyCy,
        &InverseTableau::ApplyCz,
    };
}

}  // namespace stabaal::detail
