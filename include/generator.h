#pragma once
#ifndef GENERATOR_H

#include <boost/dynamic_bitset.hpp>

class Generator {
   public:
    Generator() {}
    Generator(int qubits) : x(qubits), z(qubits) {}

    boost::dynamic_bitset<> x;
    boost::dynamic_bitset<> z;
    bool r = false;
};

#endif