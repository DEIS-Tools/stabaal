#ifndef STABAAL_DETAIL_PAULI_STRING_VIEW_HPP
#define STABAAL_DETAIL_PAULI_STRING_VIEW_HPP

#include "stabaal/types.hpp"

namespace stabaal::detail {
/**
 * A non-owning, writable view of a Pauli string's packed X and Z components.
 */
class PauliStringView {
   public:
    /**
     * Create a non-owning view of packed X and Z components.
     * Changes through the view modify the supplied storage, which must outlive the view.
     * @param xs Words storing the X components.
     * @param zs Words storing the Z components.
     */
    PauliStringView(BitWordSpan xs, BitWordSpan zs) noexcept;

    /**
     * Returns true if any X component in any of the Paulis in the string.
     * @param num_qubits Number of qubits used.
     * @return Whether any active position contains X or Y.
     */
    [[nodiscard]] bool AnyX(QubitSize num_qubits) const noexcept;

    /**
     * Get X component of pauli at index.
     * @param index Pauli position in the string.
     * @return The X component bit.
     */
    [[nodiscard]] bool GetX(QubitSize index) const noexcept;

    /**
     * Get Z component of pauli at index.
     * @param index Pauli position in the string.
     * @return The Z component bit.
     */
    [[nodiscard]] bool GetZ(QubitSize index) const noexcept;

    /**
     * Set the X component of the Pauli at index to value.
     * @param index Pauli position in the string.
     * @param value Value to write.
     */
    void SetX(QubitSize index, bool value = true) noexcept;

    /**
     * Set the Z component of the Pauli at index to value.
     * @param index Pauli position in the string.
     * @param value Value to write.
     */
    void SetZ(QubitSize index, bool value = true) noexcept;

    /**
     * Swap the X and Z components of the Pauli at index.
     * @param index Pauli position in the string.
     */
    void SwapXZ(QubitSize index) noexcept;

    /**
     * Prepare cleanup for a random measurement.
     * @param pivot Pivot position.
     * @param targets Output parameter. Contains 1 everywhere where Pauli has an X component
     * (excluding pivot), and 0 everywhere else.
     * @return Whether cleanup needs an S dagger gate on the pivot.
     */
    [[nodiscard]] bool MeasurementTargets(QubitSize pivot, BitWordSpan targets) const noexcept;

    /**
     * Update this Pauli string by applying Cx from the pivot to each target,
     * followed by S dagger on the pivot if needed.
     * @param pivot Pivot position.
     * @param targets Contains 1 everywhere where Pauli has an X component (excluding pivot), and 0
     * everywhere else.
     * @param apply_s_dagger Whether to apply S dagger after the CX gates.
     * @return Whether this Pauli string's sign must flip.
     */
    [[nodiscard]] bool CleanUp(QubitSize pivot, ConstBitWordSpan targets,
                               bool apply_s_dagger) noexcept;

    /**
     * Multiply this Pauli string by other, storing the resulting X and Z components.
     * @param other Pauli string to multiply with.
     * @param num_qubits Number of active qubits.
     * @return Exponent p in [0, 3] for the phase factor i^p of the product.
     */
    [[nodiscard]] PhaseExponent XorWithPhase(const PauliStringView& other,
                                             QubitSize num_qubits) noexcept;

   private:
    BitWordSpan xs_;
    BitWordSpan zs_;
};

}  // namespace stabaal::detail

#endif
