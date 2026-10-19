#ifndef STABAAL_BASIS_TABLE_HPP
#define STABAAL_BASIS_TABLE_HPP

#include <cstddef>
#include <vector>

#include "dynamic_bitset.hpp"
#include "pauli_string.hpp"
#include "types.hpp"

namespace stabaal {

struct BasisTable {
    explicit BasisTable(QubitIdx num_qubits);

    BasisTable(const BasisTable& other);
    BasisTable& operator=(const BasisTable& other);
    BasisTable(BasisTable&& other) noexcept = default;
    BasisTable& operator=(BasisTable&& other) noexcept = default;

    ~BasisTable() noexcept = default;

    void Add();
    std::size_t AllocatedBytes() const noexcept;

    DynamicBitset signs;
    std::vector<PauliString> gens;
};

}  // namespace stabaal

#endif
