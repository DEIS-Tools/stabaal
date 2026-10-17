#ifndef STABAAL_GENERATOR_HPP
#define STABAAL_GENERATOR_HPP

#include <boost/dynamic_bitset.hpp>

namespace stabaal {

class Generator {
   public:
    Generator();
    explicit Generator(int qubits);

    boost::dynamic_bitset<> x;
    boost::dynamic_bitset<> z;
    bool r = false;
};

}  // namespace stabaal

#endif
