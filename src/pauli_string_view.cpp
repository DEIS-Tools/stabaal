#include "stabaal/detail/pauli_string_view.hpp"

#include <algorithm>
#include <bit>
#include <cstddef>
#include <cstdint>

#include "stabaal/detail/bit_word.hpp"
#include "stabaal/types.hpp"

namespace stabaal::detail {

namespace {

__attribute__((target_clones("popcnt", "default"))) std::int64_t XorWordsWithPhase(
    const BitWordSpan xs, const BitWordSpan zs, const ConstBitWordSpan other_xs,
    const ConstBitWordSpan other_zs) noexcept {
    std::int64_t exponent{0};
    for (std::size_t i{0}; i < xs.size(); ++i) {
        const auto x{xs[i]};
        const auto z{zs[i]};
        const auto other_x{other_xs[i]};
        const auto other_z{other_zs[i]};
        const auto new_x{x ^ other_x};
        const auto new_z{z ^ other_z};
        exponent += std::popcount(z & x) + std::popcount(other_z & other_x) +
                    2 * std::popcount(z & other_x) - std::popcount(new_x & new_z);
        xs[i] = new_x;
        zs[i] = new_z;
    }

    return exponent;
}

__attribute__((target_clones("popcnt", "default"))) bool MeasurementTargetWords(
    const ConstBitWordSpan xs, const ConstBitWordSpan zs, const QubitIdx pivot,
    const BitWordSpan targets) noexcept {
    BitWord parity{0};
    for (std::size_t i{0}; i < targets.size(); ++i) {
        targets[i] = xs[i];
        if (i == pivot / kBitsPerWord) targets[i] &= ~(BitWord{1} << (pivot % kBitsPerWord));
        parity ^= zs[i] & targets[i];
    }

    return GetBit(zs, pivot) != (std::popcount(parity) % 2 != 0);
}

struct CleanupCounts {
    std::size_t z_count;
    bool overlap_parity;
};

__attribute__((target_clones("popcnt", "default"))) CleanupCounts
CleanUpWords(const BitWordSpan xs, const ConstBitWordSpan zs, const ConstBitWordSpan targets,
             const bool x) noexcept {
    CleanupCounts result{0, false};
    BitWord overlap{0};
    for (std::size_t i{0}; i < targets.size(); ++i) {
        const auto selected_z{zs[i] & targets[i]};
        result.z_count += std::popcount(selected_z);
        if (x) {
            overlap ^= xs[i] & selected_z;
            xs[i] ^= targets[i];
        }
    }

    result.overlap_parity = std::popcount(overlap) % 2 != 0;
    return result;
}

}  // namespace

PauliStringView::PauliStringView(const BitWordSpan xs, const BitWordSpan zs) noexcept
    : xs_{xs}, zs_{zs} {}

bool PauliStringView::AnyX(const std::size_t num_qubits) const noexcept {
    const auto words{xs_.first(WordCount(num_qubits))};
    return std::any_of(words.begin(), words.end(), [](const auto word) { return word != 0; });
}

bool PauliStringView::GetX(const QubitIdx idx) const noexcept { return GetBit(xs_, idx); }

bool PauliStringView::GetZ(const QubitIdx idx) const noexcept { return GetBit(zs_, idx); }

void PauliStringView::SetX(const QubitIdx idx, const bool val) noexcept { SetBit(xs_, idx, val); }

void PauliStringView::SetZ(const QubitIdx idx, const bool val) noexcept { SetBit(zs_, idx, val); }

void PauliStringView::SwapXZ(const QubitIdx idx) noexcept {
    const bool x{GetX(idx)};
    SetX(idx, GetZ(idx));
    SetZ(idx, x);
}

bool PauliStringView::MeasurementTargets(const QubitIdx pivot,
                                         const BitWordSpan targets) const noexcept {
    return MeasurementTargetWords(xs_, zs_, pivot, targets);
}

bool PauliStringView::CleanUp(const QubitIdx pivot, const ConstBitWordSpan targets,
                              const bool apply_s_dagger) noexcept {
    const bool x{GetX(pivot)};
    bool z{GetZ(pivot)};
    const auto [z_count, overlap_parity]{CleanUpWords(xs_, zs_, targets, x)};
    const bool parity{z_count % 2 != 0};
    bool flip{x && ((parity && !z) != ((z_count & 2U) != 0))};
    flip = flip != (x && overlap_parity);
    z = z != parity;
    if (apply_s_dagger) {
        flip = flip != (x && !z);
        z = z != x;
    }

    SetZ(pivot, z);
    return flip;
}

std::int64_t PauliStringView::XorWithPhase(const PauliStringView& other,
                                           const std::size_t num_qubits) noexcept {
    const auto num_words{WordCount(num_qubits)};
    return XorWordsWithPhase(xs_.first(num_words), zs_.first(num_words), other.xs_.first(num_words),
                             other.zs_.first(num_words));
}
}  // namespace stabaal::detail
