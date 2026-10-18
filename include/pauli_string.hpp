#ifndef STABAAL_PAULI_STRING_HPP
#define STABAAL_PAULI_STRING_HPP

#include <boost/dynamic_bitset.hpp>

#include "types.hpp"

namespace stabaal {

class PauliString {
   public:
    PauliString();
    explicit PauliString(QubitIdx num_qubits);

    void AddI() noexcept;

    boost::dynamic_bitset<> x_;
    boost::dynamic_bitset<> z_;
    bool sign_{false};
};

}  // namespace stabaal

#endif
