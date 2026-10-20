#include "basis_table.hpp"

#include <cstddef>

#include "pauli_string.hpp"
#include "types.hpp"

using namespace stabaal;

BasisTable::BasisTable(const QubitIdx num_qubits)
    : gens{num_qubits, PauliString{num_qubits}}, signs{num_qubits} {}

void BasisTable::add() {
    for (auto &row : gens) {
        row.addI();
    }

    gens.emplace_back(static_cast<QubitIdx>(gens.size() + 1));
    signs.pushBack();
}

std::size_t BasisTable::allocatedBytes() const noexcept {
    auto bytes{gens.capacity() * sizeof(PauliString) + signs.allocatedBytes()};
    for (const auto &gen : gens) {
        bytes += gen.allocatedBytes();
    }

    return bytes;
}
