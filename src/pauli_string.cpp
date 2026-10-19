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

std::size_t PauliString::AllocatedBytes() const noexcept {
    return x_.AllocatedBytes() + z_.AllocatedBytes();
}

bool PauliString::AnyX() const noexcept { return x_.Any(); }

bool PauliString::GetX(const QubitIdx idx) const noexcept { return x_.Get(idx); }

bool PauliString::GetZ(const QubitIdx idx) const noexcept { return z_.Get(idx); }

void PauliString::SetX(const QubitIdx idx, const bool val) noexcept { x_.Set(idx, val); }

void PauliString::SetZ(const QubitIdx idx, const bool val) noexcept { z_.Set(idx, val); }

void PauliString::SwapXZ(const QubitIdx idx) noexcept {
    const bool x{x_.Get(idx)};
    x_.Set(idx, z_.Get(idx));
    z_.Set(idx, x);
}

std::int64_t PauliString::XorWithPhase(const PauliString& other) noexcept {
    assert(Size() == other.Size());
    const auto x{x_.Words()};
    const auto z{z_.Words()};
    const auto other_x{other.x_.Words()};
    const auto other_z{other.z_.Words()};
    std::int64_t exponent{0};
    for (std::size_t i{0}; i < x.size(); ++i) {
        const auto new_x{x[i] ^ other_x[i]};
        const auto new_z{z[i] ^ other_z[i]};
        exponent += FastPopCount(z[i] & x[i]) + FastPopCount(other_z[i] & other_x[i]) +
                    2 * FastPopCount(z[i] & other_x[i]) - FastPopCount(new_x & new_z);
        x[i] = new_x;
        z[i] = new_z;
    }
    return exponent;
}

std::uint64_t PauliString::PopCount(const std::uint64_t word) const noexcept {
    __builtin_cpu_init();
    if (__builtin_cpu_supports("popcnt")) {
        return FastPopCount(word);
    }

    return NormalPopCount(word);
}

__attribute__((target("popcnt")))
std::uint64_t PauliString::FastPopCount(const std::uint64_t word) const noexcept {
    return std::popcount(word);
}

__attribute__((target("default")))
std::uint64_t PauliString::NormalPopCount(const std::uint64_t word) const noexcept {
    return std::popcount(word);
}
