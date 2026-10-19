#include "pauli_string.hpp"

#include <bit>
#include <cstddef>
#include <cstdint>

#include "types.hpp"

using namespace stabaal;

PauliString::PauliString() noexcept = default;
PauliString::PauliString(const std::size_t num_qubits) : xz_{BitsRequired(num_qubits)} {}

PauliString::PauliString(const PauliString& other, const std::size_t num_qubits)
    : xz_{other.xz_, BitsRequired(num_qubits)} {}

void PauliString::AddI(const std::size_t old_num_qubits) {
    if (old_num_qubits % kQubitsPerWord == 0) {
        xz_.Resize(BitsRequired(old_num_qubits), BitsRequired(old_num_qubits + 1));
    }
}

std::size_t PauliString::AllocatedBytes(const std::size_t num_qubits) noexcept {
    return DynamicBitset::AllocatedBytes(BitsRequired(num_qubits));
}

bool PauliString::AnyX(const std::size_t num_qubits) const noexcept {
    const auto words{xz_.Words(BitsRequired(num_qubits))};
    for (std::size_t i{0}; i < words.size(); i += 2) {
        if (words[i] != 0) return true;
    }

    return false;
}

bool PauliString::GetX(const QubitIdx idx) const noexcept { return xz_.Get(GetXIdx(idx)); }

bool PauliString::GetZ(const QubitIdx idx) const noexcept { return xz_.Get(GetZIdx(idx)); }

void PauliString::SetX(const QubitIdx idx, const bool val) noexcept { xz_.Set(GetXIdx(idx), val); }

void PauliString::SetZ(const QubitIdx idx, const bool val) noexcept { xz_.Set(GetZIdx(idx), val); }

void PauliString::SwapXZ(const QubitIdx idx) noexcept {
    const bool x{GetX(idx)};
    SetX(idx, GetZ(idx));
    SetZ(idx, x);
}

__attribute__((target_clones("popcnt", "default"))) std::int64_t PauliString::XorWithPhase(
    const PauliString& other, const std::size_t num_qubits) noexcept {
    const auto words{xz_.Words(BitsRequired(num_qubits))};
    const auto other_words{other.xz_.Words(BitsRequired(num_qubits))};
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

std::size_t PauliString::GetXIdx(const QubitIdx idx) noexcept {
    return (idx / kQubitsPerWord) * (2 * kQubitsPerWord) + idx % kQubitsPerWord;
}

std::size_t PauliString::GetZIdx(const QubitIdx idx) noexcept {
    return GetXIdx(idx) + kQubitsPerWord;
}

std::size_t PauliString::BitsRequired(const std::size_t num_qubits) noexcept {
    return (num_qubits / kQubitsPerWord + (num_qubits % kQubitsPerWord != 0 ? 1 : 0)) * 2 *
           kQubitsPerWord;
}
