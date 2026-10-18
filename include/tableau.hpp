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
using BitRef = boost::dynamic_bitset<>::reference;
using BasisTable = std::vector<Generator>;

class Tableau {
   public:
    enum class Basis : std::uint8_t {
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
    void applyX(int q) noexcept;
    void applyS(int q);
    void applyZ(int q) noexcept;

    void print() const;

   private:
    friend struct TableauTestAccess;

    static std::mt19937& randomEngine_();

    BasisTable computational_;
    BasisTable hadamard_;
    std::uint64_t qubits_;

    void cleanUp_(int pivot, const Generator& gen) noexcept;
    void collapse_(int pivot, bool p) noexcept;
    bool randomMeasurement_(int q, Basis b);

    static std::int64_t e_(const Generator& a, const Generator& b);
    static bool phase_(const Generator& a, const Generator& b, int phaseRule);

    static void swapBit_(BitRef a, BitRef b) noexcept;
};

}  // namespace stabaal

#endif
