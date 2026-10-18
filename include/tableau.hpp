#ifndef STABAAL_TABLEAU_HPP
#define STABAAL_TABLEAU_HPP

#include <cstdint>
#include <iomanip>
#include <iostream>
#include <random>
#include <sstream>
#include <vector>

#include "pauli_string.hpp"
#include "types.hpp"

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

    explicit Tableau(QubitIdx num_qubits);

    void Add();

    bool MeasureZ(QubitIdx q);
    bool MeasureX(QubitIdx q);

    void ApplyH(QubitIdx q);
    void ApplyCx(QubitIdx c, QubitIdx t);
    void ApplyX(QubitIdx q) noexcept;
    void ApplyS(QubitIdx q);
    void ApplyZ(QubitIdx q) noexcept;

    void Print() const;

   private:
    friend struct TableauTestAccess;

    static std::mt19937& RandomEngine();

    BasisTable computational_;
    BasisTable hadamard_;
    QubitIdx num_qubits_;

    void CleanUp(QubitIdx pivot, const PauliString& gen) noexcept;
    void Collapse(QubitIdx pivot, bool p) noexcept;
    bool RandomMeasurement(QubitIdx q, Basis b);

    static std::int64_t E(const PauliString& a, const PauliString& b);
    static bool Phase(const PauliString& a, const PauliString& b, std::uint8_t phase_rule);

    static void SwapBit(BitRef a, BitRef b) noexcept;
};

}  // namespace stabaal

#endif
