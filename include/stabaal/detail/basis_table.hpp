#ifndef STABAAL_DETAIL_BASIS_TABLE_HPP
#define STABAAL_DETAIL_BASIS_TABLE_HPP

#include <cstddef>

#include "dynamic_array.hpp"
#include "dynamic_bitset.hpp"
#include "pauli_string.hpp"

namespace stabaal::detail {

struct BasisTable {
    explicit BasisTable(std::size_t num_qubits);

    BasisTable(const BasisTable& other, std::size_t num_qubits);
    BasisTable(const BasisTable& other) = delete;
    BasisTable& operator=(const BasisTable& other) = delete;
    BasisTable(BasisTable&& other) noexcept = default;
    BasisTable& operator=(BasisTable&& other) noexcept = default;

    ~BasisTable() noexcept = default;

    void add(std::size_t old_num_qubits);
    static std::size_t allocatedBytes(std::size_t num_qubits) noexcept;

    DynamicBitset signs;
    DynamicArray<PauliString> gens;
};

}  // namespace stabaal::detail

#endif
