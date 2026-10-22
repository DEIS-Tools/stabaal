#include "stabaal/detail/basis_table.hpp"

#include <bit>
#include <cstddef>
#include <utility>

#include "stabaal/detail/dynamic_array.hpp"
#include "stabaal/detail/dynamic_bitset.hpp"
#include "stabaal/detail/pauli_string.hpp"

namespace stabaal::detail {

namespace {

std::size_t capacityRows(const std::size_t num_qubits) noexcept {
    return num_qubits == 0 ? 0 : std::bit_ceil(num_qubits);
}

}  // namespace

BasisTable::BasisTable(const std::size_t num_qubits)
    : signs{num_qubits}, gens{capacityRows(num_qubits)} {
    for (std::size_t i{0}; i < num_qubits; ++i) {
        gens[i] = PauliString{num_qubits};
    }
}

BasisTable::BasisTable(const BasisTable& other, const std::size_t num_qubits)
    : signs{other.signs, num_qubits}, gens{capacityRows(num_qubits)} {
    for (std::size_t i{0}; i < num_qubits; ++i) {
        gens[i] = PauliString{other.gens[i], num_qubits};
    }
}

void BasisTable::add(const std::size_t old_num_qubits) {
    const auto num_qubits{old_num_qubits + 1};
    if (capacityRows(old_num_qubits) < capacityRows(num_qubits)) {
        DynamicArray<PauliString> resized{capacityRows(num_qubits)};
        for (std::size_t i{0}; i < old_num_qubits; ++i) {
            resized[i] = std::move(gens[i]);
        }

        gens = std::move(resized);
    }

    for (std::size_t i{0}; i < old_num_qubits; ++i) {
        gens[i].addI(old_num_qubits);
    }

    gens[old_num_qubits] = PauliString{num_qubits};
    signs.resize(old_num_qubits, num_qubits);
}

std::size_t BasisTable::allocatedBytes(const std::size_t num_qubits) noexcept {
    return capacityRows(num_qubits) * sizeof(PauliString) +
           DynamicBitset::allocatedBytes(num_qubits) +
           num_qubits * PauliString::allocatedBytes(num_qubits);
}

}  // namespace stabaal::detail
