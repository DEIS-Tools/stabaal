#ifndef STABAAL_BASIS_TABLE_HPP
#define STABAAL_BASIS_TABLE_HPP

#include <cstddef>

#include "dynamic_bitset.hpp"
#include "pauli_string.hpp"
#include "types.hpp"

namespace stabaal {

struct BasisTable {
    explicit BasisTable(std::size_t num_qubits);

    BasisTable(const BasisTable& other, std::size_t num_qubits);
    BasisTable(const BasisTable& other) = delete;
    BasisTable& operator=(const BasisTable& other) = delete;
    BasisTable(BasisTable&& other) noexcept = default;
    BasisTable& operator=(BasisTable&& other) noexcept = default;

    ~BasisTable() noexcept = default;

    void Add(std::size_t old_num_qubits);
    static std::size_t AllocatedBytes(std::size_t num_qubits) noexcept;

    DynamicBitset signs;
    DynamicArray<PauliString> gens;
};

}  // namespace stabaal

#endif
