#ifndef STABAAL_TABLEAU_HPP
#define STABAAL_TABLEAU_HPP

#include <cstdint>
#include <iomanip>
#include <iostream>
#include <random>
#include <sstream>
#include <vector>

#include "generator.hpp"

namespace stabaal {

using std::int64_t;
using std::uint64_t;
using std::uint8_t;

using BitRef = boost::dynamic_bitset<>::reference;
using BasisTable = std::vector<Generator>;

class Tableau {
   public:
    enum class Basis : uint8_t {
        COMPUTATIONAL,
        HADAMARD,
        IMAGINARY,
    };

    explicit Tableau(int qubitCount);

    void add();

    bool measureZ(int q);
    bool measureX(int q);

    void applyH(int q);
    void applyCx(int c, int t);
    void applyX(int q);
    void applyS(int q);
    void applyZ(int q);

    void print() const;

   private:
    friend struct TableauTestAccess;

    static std::mt19937 &_randomEngine();

    BasisTable _computational;
    BasisTable _hadamard;
    uint64_t _qubits;

    void _cleanUp(int pivot, const Generator &gen);
    void _collapse(int pivot, bool p);
    bool _randomMeasurement(int q, Basis b);

    static int64_t _e(const Generator &a, const Generator &b);
    static bool _phase(const Generator &a, const Generator &b, int phaseRule);

    static void _swapBit(BitRef a, BitRef b);
};

}  // namespace stabaal

#endif
