#ifndef STABAAL_TABLEAU_HPP
#define STABAAL_TABLEAU_HPP

#include <cstdint>
#include <iomanip>
#include <iostream>
#include <random>
#include <sstream>
#include <vector>

#include "pauli_string.hpp"

namespace stabaal {
using BitRef = boost::dynamic_bitset<>::reference;
using BasisTable = std::vector<PauliString>;

class Tableau {
   public:
    enum class Basis : std::uint8_t {
        kComputational,
        kHadamard,
        kImaginary,
    };

    explicit Tableau(int num_qubits);

    void Add();

    bool MeasureZ(int q);
    bool MeasureX(int q);

    void ApplyH(int q);
    void ApplyCx(int c, int t);
    void ApplyX(int q) noexcept;
    void ApplyS(int q);
    void ApplyZ(int q) noexcept;

    void Print() const;

   private:
    friend struct TableauTestAccess;

    static std::mt19937& RandomEngine();

    BasisTable computational_;
    BasisTable hadamard_;
    std::uint64_t num_qubits_;

    void CleanUp(int pivot, const PauliString& gen) noexcept;
    void Collapse(int pivot, bool p) noexcept;
    bool RandomMeasurement(int q, Basis b);

    static std::int64_t E(const PauliString& a, const PauliString& b);
    static bool Phase(const PauliString& a, const PauliString& b, int phase_rule);

    static void SwapBit(BitRef a, BitRef b) noexcept;
};

}  // namespace stabaal

#endif
