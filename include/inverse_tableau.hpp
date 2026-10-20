#ifndef STABAAL_INVERSE_TABLEAU_HPP
#define STABAAL_INVERSE_TABLEAU_HPP

#include <cstddef>
#include <cstdint>

#include "basis_table.hpp"
#include "pauli_string.hpp"
#include "types.hpp"

namespace stabaal {

class InverseTableau {
   public:
    enum class Basis : std::uint8_t {
        kComputational,
        kHadamard,
        kImaginary,
    };

    explicit InverseTableau(QubitIdx num_qubits);

    InverseTableau(const InverseTableau& other);
    InverseTableau& operator=(const InverseTableau& other);
    InverseTableau(InverseTableau&& other) noexcept = default;
    InverseTableau& operator=(InverseTableau&& other) noexcept = default;
    ~InverseTableau() noexcept = default;

    void Add();

    bool MeasureX(QubitIdx q);
    bool MeasureZ(QubitIdx q);

    void ApplyH(QubitIdx q);
    void ApplyS(QubitIdx q);
    void ApplyX(QubitIdx q);
    void ApplyZ(QubitIdx q);
    void ApplyCx(QubitIdx c, QubitIdx t);

    void Print() const;

    std::size_t MemoryUsage() const noexcept;

   private:
    BasisTable computational_;
    BasisTable hadamard_;
    QubitIdx num_qubits_;

    void CleanUp(QubitIdx pivot, const PauliString& gen) noexcept;
    void Collapse(QubitIdx pivot, bool p) noexcept;
    bool RandomMeasurement(QubitIdx q, Basis b);
};

}  // namespace stabaal

#endif
