#ifndef STABAAL_DETAIL_INVERSE_TABLEAU_HPP
#define STABAAL_DETAIL_INVERSE_TABLEAU_HPP

#include <cstddef>
#include <cstdint>
#include <memory>

#include "pauli_string_view.hpp"
#include "stabaal/types.hpp"

namespace stabaal::detail {

class InverseTableau {
   public:
    enum class Basis : std::uint8_t {
        kComputational,
        kHadamard,
    };

    explicit InverseTableau(QubitIdx num_qubits);

    InverseTableau(const InverseTableau& other);
    InverseTableau& operator=(const InverseTableau& other);
    InverseTableau(InverseTableau&& other) noexcept = default;
    InverseTableau& operator=(InverseTableau&& other) noexcept = default;
    ~InverseTableau() noexcept = default;

    [[nodiscard]] QubitIdx NumQubits() const noexcept { return num_qubits_; }

    void Add();

    [[nodiscard]] bool MeasureX(QubitIdx q);
    [[nodiscard]] bool MeasureZ(QubitIdx q);

    void ApplyH(QubitIdx q) noexcept;
    void ApplyS(QubitIdx q) noexcept;
    void ApplyX(QubitIdx q) noexcept;
    void ApplyZ(QubitIdx q) noexcept;
    void ApplyCx(QubitIdx c, QubitIdx t) noexcept;

    void Print() const;

    [[nodiscard]] std::size_t MemoryUsage() const noexcept;

   private:
    std::unique_ptr<BitWord[]> words_;  // NOLINT(cppcoreguidelines-avoid-c-arrays)
    QubitIdx num_qubits_;
    QubitIdx capacity_;

    [[nodiscard]] static std::size_t StorageWords(QubitIdx capacity) noexcept;
    [[nodiscard]] BitWordSpan RowWords(std::size_t row) noexcept;
    [[nodiscard]] ConstBitWordSpan RowWords(std::size_t row) const noexcept;
    [[nodiscard]] BitWordSpan Flags(std::size_t index) noexcept;
    [[nodiscard]] ConstBitWordSpan Flags(std::size_t index) const noexcept;
    [[nodiscard]] std::size_t PhysicalBasis(QubitIdx q, Basis basis) const noexcept;
    [[nodiscard]] PauliStringView Generator(QubitIdx q, Basis basis) noexcept;
    [[nodiscard]] PauliStringView PhysicalGenerator(QubitIdx q, std::size_t basis) noexcept;
    [[nodiscard]] bool Sign(QubitIdx q, Basis basis) const noexcept;
    void FlipSign(QubitIdx q, Basis basis) noexcept;
    void CopyRows(const InverseTableau& other) noexcept;

    [[nodiscard]] bool Measure(QubitIdx q, Basis basis);

    void CleanUp(QubitIdx pivot, ConstBitWordSpan targets, bool apply_s_dagger) noexcept;
    void Collapse(QubitIdx pivot, bool p) noexcept;
    [[nodiscard]] bool RandomMeasurement(QubitIdx q, Basis b);
};

}  // namespace stabaal::detail

#endif
