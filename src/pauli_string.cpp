#include "pauli_string.hpp"

#include <bit>
#include <cassert>
#include <cstddef>
#include <cstdint>

#include "types.hpp"

using namespace stabaal;

PauliString::PauliString() noexcept = default;
PauliString::PauliString(const QubitIdx num_qubits) : x_{num_qubits}, z_{num_qubits} {}

void PauliString::addI() {
    x_.pushBack();
    z_.pushBack();
}

std::size_t PauliString::size() const noexcept { return x_.size(); }

std::size_t PauliString::allocatedBytes() const noexcept {
    return x_.allocatedBytes() + z_.allocatedBytes();
}

bool PauliString::anyX() const noexcept { return x_.any(); }

bool PauliString::getX(const QubitIdx idx) const noexcept { return x_.get(idx); }

bool PauliString::getZ(const QubitIdx idx) const noexcept { return z_.get(idx); }

void PauliString::setX(const QubitIdx idx, const bool val) noexcept { x_.set(idx, val); }

void PauliString::setZ(const QubitIdx idx, const bool val) noexcept { z_.set(idx, val); }

void PauliString::swapXZ(const QubitIdx idx) noexcept {
    const bool x{x_.get(idx)};
    x_.set(idx, z_.get(idx));
    z_.set(idx, x);
}

__attribute__((target_clones("popcnt", "default"))) std::int64_t PauliString::xorWithPhase(
    const PauliString &other) noexcept {
    assert(size() == other.size());
    const auto x{x_.words()};
    const auto z{z_.words()};
    const auto other_x{other.x_.words()};
    const auto other_z{other.z_.words()};
    std::int64_t exponent{0};
    for (std::size_t i{0}; i < x.size(); ++i) {
        const auto new_x{x[i] ^ other_x[i]};
        const auto new_z{z[i] ^ other_z[i]};
        exponent += std::popcount(z[i] & x[i]) + std::popcount(other_z[i] & other_x[i]) +
                    2 * std::popcount(z[i] & other_x[i]) - std::popcount(new_x & new_z);
        x[i] = new_x;
        z[i] = new_z;
    }

    return exponent;
}
