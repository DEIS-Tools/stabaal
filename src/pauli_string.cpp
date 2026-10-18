#include "pauli_string.hpp"

#include "types.hpp"

using namespace stabaal;

PauliString::PauliString() = default;
PauliString::PauliString(QubitIdx num_qubits) : x_(num_qubits), z_(num_qubits) {}

void PauliString::AddI() noexcept {
    x_.push_back(false);
    z_.push_back(false);
}
