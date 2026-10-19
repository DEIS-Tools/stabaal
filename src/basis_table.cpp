#include "basis_table.hpp"

#include <cstddef>

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
    signs.Resize(gens.size());
}

std::size_t BasisTable::AllocatedBytes() const noexcept {
    auto bytes{gens.capacity() * sizeof(PauliString) + signs.AllocatedBytes()};
    for (const auto& gen : gens) {
        bytes += gen.AllocatedBytes();
    }

    return bytes;
}
