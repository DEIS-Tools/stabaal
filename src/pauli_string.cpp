#include "pauli_string.hpp"

#include <bit>
#include <cassert>
#include <cstddef>
#include <cstdint>

#include "types.hpp"

using namespace stabaal;

PauliString::PauliString() noexcept = default;
PauliString::PauliString(const QubitIdx num_qubits)
    : xz_{(num_qubits / kQubitsPerWord + (num_qubits % kQubitsPerWord != 0 ? 1 : 0)) * 2 *
          kQubitsPerWord},
      num_qubits_{num_qubits} {}

void PauliString::AddI() {
    if (num_qubits_ % kQubitsPerWord == 0) {
        xz_.Resize(2 * (num_qubits_ + kQubitsPerWord));
    }

    ++num_qubits_;
}

std::size_t PauliString::Size() const noexcept { return num_qubits_; }

std::size_t PauliString::AllocatedBytes() const noexcept { return xz_.AllocatedBytes(); }

bool PauliString::AnyX() const noexcept {
    const auto words{xz_.Words()};
    for (std::size_t i{0}; i < words.size(); i += 2) {
        if (words[i] != 0) return true;
    }

    return false;
}

bool PauliString::GetX(const QubitIdx idx) const noexcept {
    assert(idx < num_qubits_);
    return xz_.Get(GetXIdx(idx));
}

bool PauliString::GetZ(const QubitIdx idx) const noexcept {
    assert(idx < num_qubits_);
    return xz_.Get(GetZIdx(idx));
}

void PauliString::SetX(const QubitIdx idx, const bool val) noexcept {
    assert(idx < num_qubits_);
    xz_.Set(GetXIdx(idx), val);
}

void PauliString::SetZ(const QubitIdx idx, const bool val) noexcept {
    assert(idx < num_qubits_);
    xz_.Set(GetZIdx(idx), val);
}

void PauliString::SwapXZ(const QubitIdx idx) noexcept {
    const bool x{GetX(idx)};
    SetX(idx, GetZ(idx));
    SetZ(idx, x);
}

__attribute__((target_clones("popcnt", "default"))) std::int64_t PauliString::XorWithPhase(
    const PauliString& other) noexcept {
    assert(Size() == other.Size());
    const auto words{xz_.Words()};
    const auto other_words{other.xz_.Words()};
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
