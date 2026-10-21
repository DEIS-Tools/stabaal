#ifndef STABAAL_DETAIL_INVERSE_TABLEAU_HPP
#define STABAAL_DETAIL_INVERSE_TABLEAU_HPP

#include <cstddef>
#include <cstdint>
#include <memory>
#include <span>

#include "pauli_string.hpp"
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

    void Add();

    [[nodiscard]] bool MeasureX(QubitIdx q);
    [[nodiscard]] bool MeasureZ(QubitIdx q);

    void ApplyH(QubitIdx q);
    void ApplyS(QubitIdx q);
    void ApplyX(QubitIdx q);
    void ApplyZ(QubitIdx q);
    void ApplyCx(QubitIdx c, QubitIdx t);

    void Print() const;

    [[nodiscard]] std::size_t MemoryUsage() const noexcept;

   private:
    std::unique_ptr<std::uint64_t[]> words_;  // NOLINT(cppcoreguidelines-avoid-c-arrays)
    QubitIdx num_qubits_;
    QubitIdx capacity_;

    [[nodiscard]] static std::size_t StorageWords(QubitIdx capacity) noexcept;
    [[nodiscard]] std::span<std::uint64_t> RowWords(std::size_t row) noexcept;
    [[nodiscard]] std::span<const std::uint64_t> RowWords(std::size_t row) const noexcept;
    [[nodiscard]] std::span<std::uint64_t> Flags(std::size_t index) noexcept;
    [[nodiscard]] std::span<const std::uint64_t> Flags(std::size_t index) const noexcept;
    [[nodiscard]] std::size_t PhysicalBasis(QubitIdx q, Basis basis) const noexcept;
    [[nodiscard]] PauliStringView Generator(QubitIdx q, Basis basis) noexcept;
    [[nodiscard]] bool Sign(QubitIdx q, Basis basis) const noexcept;
    void FlipSign(QubitIdx q, Basis basis) noexcept;
    void SetSign(QubitIdx q, Basis basis, bool value) noexcept;
    void CopyRows(const InverseTableau& other) noexcept;

    void CheckQubit(QubitIdx q) const;
    [[nodiscard]] bool Measure(QubitIdx q, Basis basis);

    void CleanUp(QubitIdx pivot, std::span<const std::uint64_t> targets,
                 bool apply_s_dagger) noexcept;
    void Collapse(QubitIdx pivot, bool p) noexcept;
    [[nodiscard]] bool RandomMeasurement(QubitIdx q, Basis b);
};

}  // namespace stabaal::detail

#endif
