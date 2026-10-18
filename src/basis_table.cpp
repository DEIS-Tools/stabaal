#include "basis_table.hpp"

#include "pauli_string.hpp"
#include "types.hpp"

using namespace stabaal;

BasisTable::BasisTable(const QubitIdx num_qubits)
    : gens{num_qubits, PauliString{num_qubits}}, signs{num_qubits} {}

void BasisTable::Add() {
    for (auto& row : gens) {
        row.AddI();
    }

    gens.emplace_back(static_cast<QubitIdx>(gens.size() + 1));
    signs.PushBack();
}
