#ifndef STABAAL_DETAIL_PAULI_STRING_VIEW_HPP
#define STABAAL_DETAIL_PAULI_STRING_VIEW_HPP

#include <cstddef>
#include <cstdint>

#include "stabaal/types.hpp"

namespace stabaal::detail {
class PauliStringView {
   public:
    PauliStringView(BitWordSpan xs, BitWordSpan zs) noexcept;

    [[nodiscard]] bool AnyX(std::size_t num_qubits) const noexcept;

    [[nodiscard]] bool GetX(QubitIdx idx) const noexcept;
    [[nodiscard]] bool GetZ(QubitIdx idx) const noexcept;

    void SetX(QubitIdx idx, bool val = true) noexcept;
    void SetZ(QubitIdx idx, bool val = true) noexcept;
    void SwapXZ(QubitIdx idx) noexcept;

    [[nodiscard]] bool MeasurementTargets(QubitIdx pivot, BitWordSpan targets) const noexcept;
    [[nodiscard]] bool CleanUp(QubitIdx pivot, ConstBitWordSpan targets,
                               bool apply_s_dagger) noexcept;

    [[nodiscard]] std::int64_t XorWithPhase(const PauliStringView& other,
                                            std::size_t num_qubits) noexcept;

   private:
    BitWordSpan xs_;
    BitWordSpan zs_;
};

}  // namespace stabaal::detail

#endif
