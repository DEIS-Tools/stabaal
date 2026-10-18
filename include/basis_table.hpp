#ifndef STABAAL_BASIS_TABLE_HPP
#define STABAAL_BASIS_TABLE_HPP

#include <boost/dynamic_bitset.hpp>
#include <vector>

#include "pauli_string.hpp"
#include "types.hpp"

namespace stabaal {

struct BasisTable {
    explicit BasisTable(QubitIdx num_qubits);

    void Add();

    std::vector<PauliString> gens;
    boost::dynamic_bitset<> signs;
};

}  // namespace stabaal

#endif
