#include "basis_table.hpp"

#include <bit>
#include <cstddef>
#include <utility>

#include "pauli_string.hpp"
#include "types.hpp"

using namespace stabaal;

namespace {

std::size_t CapacityRows(const std::size_t num_qubits) noexcept {
    return num_qubits == 0 ? 0 : std::bit_ceil(num_qubits);
}

}  // namespace

BasisTable::BasisTable(const std::size_t num_qubits)
    : signs{num_qubits}, gens{MakeDynamicArray<PauliString>(CapacityRows(num_qubits))} {
    for (std::size_t i{0}; i < num_qubits; ++i) {
        gens[i] = PauliString{num_qubits};
    }
}

BasisTable::BasisTable(const BasisTable& other, const std::size_t num_qubits)
    : signs{other.signs, num_qubits},
      gens{MakeDynamicArray<PauliString>(CapacityRows(num_qubits))} {
    for (std::size_t i{0}; i < num_qubits; ++i) {
        gens[i] = PauliString{other.gens[i], num_qubits};
    }
}

void BasisTable::Add(const std::size_t old_num_qubits) {
    const auto num_qubits{old_num_qubits + 1};
    if (CapacityRows(old_num_qubits) != CapacityRows(num_qubits)) {
        auto resized{MakeDynamicArray<PauliString>(CapacityRows(num_qubits))};
        for (std::size_t i{0}; i < old_num_qubits; ++i) {
            resized[i] = std::move(gens[i]);
        }
        
        gens = std::move(resized);
    }

    for (std::size_t i{0}; i < old_num_qubits; ++i) {
        gens[i].AddI(old_num_qubits);
    }

    gens[old_num_qubits] = PauliString{num_qubits};
    signs.Resize(old_num_qubits, num_qubits);
}

std::size_t BasisTable::AllocatedBytes(const std::size_t num_qubits) noexcept {
    return CapacityRows(num_qubits) * sizeof(PauliString) +
           DynamicBitset::AllocatedBytes(num_qubits) +
           num_qubits * PauliString::AllocatedBytes(num_qubits);
}
