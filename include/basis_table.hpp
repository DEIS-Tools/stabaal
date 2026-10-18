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

    void Add();
    std::size_t AllocatedBytes() const noexcept;

    std::vector<PauliString> gens;
    DynamicBitset signs;
};

}  // namespace stabaal

#endif
