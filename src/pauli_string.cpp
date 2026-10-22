#include "stabaal/detail/pauli_string.hpp"

#include <bit>
#include <cstddef>
#include <cstdint>
#include <span>

#include "stabaal/detail/dynamic_bitset.hpp"
#include "stabaal/detail/types.hpp"

namespace stabaal::detail {

namespace {

constexpr std::size_t kQubitsPerWord{64};

constexpr std::size_t getXIdx(const QubitIdx idx) noexcept {
    return (idx / kQubitsPerWord) * (2 * kQubitsPerWord) + idx % kQubitsPerWord;
}

constexpr std::size_t getZIdx(const QubitIdx idx) noexcept { return getXIdx(idx) + kQubitsPerWord; }

constexpr std::size_t bitsRequired(const std::size_t num_qubits) noexcept {
    return (num_qubits / kQubitsPerWord + (num_qubits % kQubitsPerWord != 0 ? 1 : 0)) * 2 *
           kQubitsPerWord;
}

// Use hardware POPCNT instruction if available
__attribute__((target_clones("popcnt", "default"))) std::int64_t xorWordsWithPhase(
    const std::span<std::uint64_t> words,
    const std::span<const std::uint64_t> other_words) noexcept {
    std::int64_t exponent{0};
    for (std::size_t i{0}; i < words.size(); i += 2) {
        const auto x{words[i]};
        const auto z{words[i + 1]};
        const auto other_x{other_words[i]};
        const auto other_z{other_words[i + 1]};
        const auto new_x{x ^ other_x};
        const auto new_z{z ^ other_z};
        exponent += std::popcount(z & x) + std::popcount(other_z & other_x) +
                    2 * std::popcount(z & other_x) - std::popcount(new_x & new_z);
        words[i] = new_x;
        words[i + 1] = new_z;
    }

    return exponent;
}

}  // namespace

PauliString::PauliString() noexcept = default;
PauliString::PauliString(const std::size_t num_qubits) : xz_{bitsRequired(num_qubits)} {}

PauliString::PauliString(const PauliString& other, const std::size_t num_qubits)
    : xz_{other.xz_, bitsRequired(num_qubits)} {}

void PauliString::addI(const std::size_t old_num_qubits) {
    if (old_num_qubits % kQubitsPerWord == 0) {
        xz_.resize(bitsRequired(old_num_qubits), bitsRequired(old_num_qubits + 1));
    }
}

std::size_t PauliString::allocatedBytes(const std::size_t num_qubits) noexcept {
    return DynamicBitset::allocatedBytes(bitsRequired(num_qubits));
}

bool PauliString::anyX(const std::size_t num_qubits) const noexcept {
    const auto words{xz_.words(bitsRequired(num_qubits))};
    for (std::size_t i{0}; i < words.size(); i += 2) {
        if (words[i] != 0) return true;
    }

    return false;
}

bool PauliString::getX(const QubitIdx idx) const noexcept { return xz_.get(getXIdx(idx)); }

bool PauliString::getZ(const QubitIdx idx) const noexcept { return xz_.get(getZIdx(idx)); }

void PauliString::setX(const QubitIdx idx, const bool val) noexcept { xz_.set(getXIdx(idx), val); }

void PauliString::setZ(const QubitIdx idx, const bool val) noexcept { xz_.set(getZIdx(idx), val); }

void PauliString::swapXZ(const QubitIdx idx) noexcept {
    const bool x{getX(idx)};
    setX(idx, getZ(idx));
    setZ(idx, x);
}

void PauliString::copyFrom(const PauliString& other, const std::size_t num_qubits) noexcept {
    xz_.copyFrom(other.xz_, bitsRequired(num_qubits));
}

std::int64_t PauliString::xorWithPhase(const PauliString& other,
                                       const std::size_t num_qubits) noexcept {
    return xorWordsWithPhase(xz_.words(bitsRequired(num_qubits)),
                             other.xz_.words(bitsRequired(num_qubits)));
}

}  // namespace stabaal::detail
