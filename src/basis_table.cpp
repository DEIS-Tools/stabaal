#include "basis_table.hpp"

#include <cstddef>

#include "pauli_string.hpp"
#include "types.hpp"

using namespace stabaal;

BasisTable::BasisTable(const QubitIdx num_qubits) : signs{num_qubits} {
    gens.reserve(num_qubits);
    for (QubitIdx i{0}; i < num_qubits; ++i) {
        gens.emplace_back(num_qubits);
    }
}

BasisTable::BasisTable(const BasisTable& other) : signs{other.signs, other.gens.size()} {
    gens.reserve(other.gens.size());
    for (const auto& gen : other.gens) {
        gens.emplace_back(gen, other.gens.size());
    }
}

BasisTable& BasisTable::operator=(const BasisTable& other) {
    if (this != &other) *this = BasisTable{other};
    return *this;
}

void BasisTable::Add() {
    for (auto& row : gens) {
        row.AddI(gens.size());
    }

    gens.emplace_back(gens.size() + 1);
    signs.Resize(gens.size() - 1, gens.size());
}

std::size_t BasisTable::AllocatedBytes() const noexcept {
    auto bytes{gens.capacity() * sizeof(PauliString) + DynamicBitset::AllocatedBytes(gens.size())};
    bytes += gens.size() * PauliString::AllocatedBytes(gens.size());

    return bytes;
}
