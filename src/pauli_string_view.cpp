#include "stabaal/detail/pauli_string_view.hpp"

#include <algorithm>
#include <bit>
#include <cstddef>

#include "stabaal/detail/bit_ops.hpp"
#include "stabaal/types.hpp"

namespace stabaal::detail {

namespace {

/**
 * Multiply packed Pauli strings, storing the product's components in xs and zs.
 * Uses POPCNT when supported, with a default CPU fallback.
 * @param xs Writable X components of the left-hand string.
 * @param zs Writable Z components of the left-hand string.
 * @param other_xs X components of the right-hand string.
 * @param other_zs Z components of the right-hand string.
 * @return Exponent p in [0, 3] for the product's phase factor i^p, excluding stored signs.
 */
__attribute__((target_clones("popcnt", "default"))) PhaseExponent
XorWordsWithPhase(const BitWordSpan xs, const BitWordSpan zs, const ConstBitWordSpan other_xs,
                  const ConstBitWordSpan other_zs) noexcept {
    PhaseExponent exponent{0};
    BitWord negative_phases{0};
    for (std::size_t i{0}; i < xs.size(); ++i) {
        const auto x{xs[i]};
        const auto z{zs[i]};
        const auto other_x{other_xs[i]};
        const auto other_z{other_zs[i]};
        const auto new_x{x ^ other_x};
        const auto new_z{z ^ other_z};
        const auto cross_xz{x & other_z};
        const auto cross_zx{z & other_x};
        const auto anticommuting{cross_xz ^ cross_zx};
        // Anticommuting Paulis contribute i, with another factor of -1 for XZ, ZY, YX.
        exponent += std::popcount(anticommuting);
        negative_phases ^= anticommuting & ((new_x & other_z) ^ (x & z));
        xs[i] = new_x;
        zs[i] = new_z;
    }

    return (exponent + PhaseExponent{2} * (std::popcount(negative_phases) % 2)) % 4;
}

/**
 * Fill measurement targets and determine whether the pivot needs S dagger.
 * Uses POPCNT when supported, with a default CPU fallback.
 * @param xs X components of the measured generator.
 * @param zs Z components of the measured generator.
 * @param pivot Bit position with an X component, excluded from the targets.
 * @param targets Output words containing the X-component bits with the pivot bit cleared.
 * @return Whether cleanup needs S dagger on the pivot.
 */
__attribute__((target_clones("popcnt", "default"))) bool MeasurementTargetWords(
    const ConstBitWordSpan xs, const ConstBitWordSpan zs, const QubitSize pivot,
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

/**
 * Count target Z components and toggle target X components if the pivot has an X component.
 * Uses POPCNT when supported, with a default CPU fallback.
 * @param xs Writable X components of the generator being updated.
 * @param zs Z components of the generator being updated.
 * @param targets Packed bits marking cleanup target positions.
 * @param x X component at the pivot.
 * @return Target Z-component count and, if x is true, the parity of target X/Z overlaps
 * before the update. The overlap parity is false when x is false.
 */
__attribute__((target_clones("popcnt", "default"))) CleanupCounts
CleanUpWords(const BitWordSpan xs, const ConstBitWordSpan zs, const ConstBitWordSpan targets,
             const bool x) noexcept {
    CleanupCounts result{0, false};
    BitWord overlap{0};
    for (std::size_t i{0}; i < targets.size(); ++i) {
        const auto selected_z{zs[i] & targets[i]};
        result.z_count += static_cast<std::size_t>(std::popcount(selected_z));
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

bool PauliStringView::AnyX(const QubitSize num_qubits) const noexcept {
    const auto words{xs_.first(WordCount(num_qubits))};
    return std::any_of(words.begin(), words.end(), [](const auto word) { return word != 0; });
}

bool PauliStringView::GetX(const QubitSize index) const noexcept { return GetBit(xs_, index); }

bool PauliStringView::GetZ(const QubitSize index) const noexcept { return GetBit(zs_, index); }

void PauliStringView::SetX(const QubitSize index, const bool value) noexcept {
    SetBit(xs_, index, value);
}

void PauliStringView::SetZ(const QubitSize index, const bool value) noexcept {
    SetBit(zs_, index, value);
}

void PauliStringView::SwapXZ(const QubitSize index) noexcept {
    const bool x{GetX(index)};
    SetX(index, GetZ(index));
    SetZ(index, x);
}

bool PauliStringView::MeasurementTargets(const QubitSize pivot,
                                         const BitWordSpan targets) const noexcept {
    return MeasurementTargetWords(xs_, zs_, pivot, targets);
}

bool PauliStringView::CleanUp(const QubitSize pivot, const ConstBitWordSpan targets,
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

PhaseExponent PauliStringView::XorWithPhase(const PauliStringView& other,
                                            const QubitSize num_qubits) noexcept {
    const auto num_words{WordCount(num_qubits)};
    return XorWordsWithPhase(xs_.first(num_words), zs_.first(num_words), other.xs_.first(num_words),
                             other.zs_.first(num_words));
}

}  // namespace stabaal::detail
