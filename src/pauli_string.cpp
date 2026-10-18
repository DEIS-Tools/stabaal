#include "pauli_string.hpp"

#include <cstddef>
#include <cstdint>

#include "types.hpp"

using namespace stabaal;

PauliString::PauliString() = default;
PauliString::PauliString(QubitIdx num_qubits) : x_{num_qubits}, z_{num_qubits} {}

void PauliString::AddI() noexcept {
    x_.push_back(false);
    z_.push_back(false);
}

std::size_t PauliString::Size() const noexcept { return x_.size(); }

bool PauliString::AnyX() const noexcept { return x_.any(); }

bool PauliString::GetX(QubitIdx idx) const noexcept { return x_[idx]; }

bool PauliString::GetZ(QubitIdx idx) const noexcept { return z_[idx]; }

void PauliString::SetX(QubitIdx idx, bool val) noexcept { x_[idx] = val; }

void PauliString::SetZ(QubitIdx idx, bool val) noexcept { z_[idx] = val; }

void PauliString::SwapXZ(QubitIdx idx) noexcept {
    const bool x{x_[idx]};
    x_[idx] = z_[idx];
    z_[idx] = x;
}

void PauliString::XorWith(const PauliString& other) noexcept {
    x_ ^= other.x_;
    z_ ^= other.z_;
}

std::int64_t PauliString::PhaseExponent(const PauliString& other) const {
    return static_cast<std::int64_t>((z_ & x_).count()) +
           static_cast<std::int64_t>((other.z_ & other.x_).count()) +
           2 * static_cast<std::int64_t>((z_ & other.x_).count()) -
           static_cast<std::int64_t>(((x_ ^ other.x_) & (z_ ^ other.z_)).count());
}
