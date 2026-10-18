#include "pauli_string.hpp"

using namespace stabaal;

PauliString::PauliString() = default;
PauliString::PauliString(int num_qubits) : x_(num_qubits), z_(num_qubits) {}
