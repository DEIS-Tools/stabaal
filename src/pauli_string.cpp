#include "pauli_string.hpp"

#include <bit>
#include <cassert>
#include <cstddef>
#include <cstdint>

#include "types.hpp"

using namespace stabaal;

PauliString::PauliString() noexcept = default;
PauliString::PauliString(const QubitIdx num_qubits) : x_{num_qubits}, z_{num_qubits} {}

void PauliString::AddI() {
    x_.PushBack();
    z_.PushBack();
}

std::size_t PauliString::Size() const noexcept { return x_.Size(); }

bool PauliString::AnyX() const noexcept { return x_.Any(); }

bool PauliString::GetX(const QubitIdx idx) const { return x_.Get(idx); }

bool PauliString::GetZ(const QubitIdx idx) const { return z_.Get(idx); }

void PauliString::SetX(const QubitIdx idx, const bool val) { x_.Set(idx, val); }

void PauliString::SetZ(const QubitIdx idx, const bool val) { z_.Set(idx, val); }

void PauliString::SwapXZ(const QubitIdx idx) {
    const bool x{x_.Get(idx)};
    x_.Set(idx, z_.Get(idx));
    z_.Set(idx, x);
}

PauliString& PauliString::operator^=(const PauliString& other) noexcept {
    x_ ^= other.x_;
    z_ ^= other.z_;
    return *this;
}

std::int64_t PauliString::PhaseExponent(const PauliString& other) const noexcept {
    assert(Size() == other.Size());
    const auto x{x_.Words()};
    const auto z{z_.Words()};
    const auto other_x{other.x_.Words()};
    const auto other_z{other.z_.Words()};
    std::int64_t exponent{0};
    for (std::size_t i{0}; i < x.size(); ++i) {
        exponent += std::popcount(z[i] & x[i]) + std::popcount(other_z[i] & other_x[i]) +
                    2 * std::popcount(z[i] & other_x[i]) -
                    std::popcount((x[i] ^ other_x[i]) & (z[i] ^ other_z[i]));
    }
    
    return exponent;
}
