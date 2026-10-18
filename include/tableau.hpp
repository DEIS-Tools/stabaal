#ifndef STABAAL_TABLEAU_HPP
#define STABAAL_TABLEAU_HPP

#include <cstdint>
#include <iomanip>
#include <iostream>
#include <random>
#include <sstream>

#include "basis_table.hpp"
#include "pauli_string.hpp"
#include "types.hpp"

namespace stabaal {

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
    void ApplyX(QubitIdx q);
    void ApplyS(QubitIdx q);
    void ApplyZ(QubitIdx q);

    void Print() const;

   private:
    static std::mt19937& RandomEngine();

    BasisTable computational_;
    BasisTable hadamard_;
    QubitIdx num_qubits_;

    void CleanUp(QubitIdx pivot, const PauliString& gen) noexcept;
    void Collapse(QubitIdx pivot, bool p) noexcept;
    bool RandomMeasurement(QubitIdx q, Basis b);

    static bool Phase(const PauliString& a, bool a_sign, const PauliString& b, bool b_sign,
                      std::uint8_t phase_rule) noexcept;
};

}  // namespace stabaal

#endif
