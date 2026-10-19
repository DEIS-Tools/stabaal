#include "basis_table.hpp"

#include <cstddef>
#include <utility>

#include "pauli_string.hpp"
#include "types.hpp"

using namespace stabaal;

BasisTable::BasisTable(const QubitIdx num_qubits)
    : signs{num_qubits}, gens{num_qubits, PauliString{num_qubits}} {}

BasisTable::BasisTable(const BasisTable& other)
    : signs{other.signs, other.gens.size()}, gens{other.gens} {}

BasisTable& BasisTable::operator=(const BasisTable& other) {
    if (this != &other) *this = BasisTable{other};
    return *this;
}

void BasisTable::Add() {
    for (auto& row : gens) {
        row.AddI();
    }

    gens.emplace_back(static_cast<QubitIdx>(gens.size() + 1));
    signs.Resize(gens.size() - 1, gens.size());
}

std::size_t BasisTable::AllocatedBytes() const noexcept {
    auto bytes{gens.capacity() * sizeof(PauliString) + signs.AllocatedBytes(gens.size())};
    for (const auto& gen : gens) {
        bytes += gen.AllocatedBytes();
    }

    return bytes;
}
