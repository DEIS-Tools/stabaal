#ifndef STABAAL_DETAIL_PAULI_STRING_HPP
#define STABAAL_DETAIL_PAULI_STRING_HPP

#include <cstddef>
#include <cstdint>
#include <span>

#include "stabaal/types.hpp"

namespace stabaal::detail {
class PauliStringView {
   public:
    PauliStringView(std::span<std::uint64_t> xs, std::span<std::uint64_t> zs) noexcept;

    [[nodiscard]] bool AnyX(std::size_t num_qubits) const noexcept;

    [[nodiscard]] bool GetX(QubitIdx idx) const noexcept;
    [[nodiscard]] bool GetZ(QubitIdx idx) const noexcept;

    void SetX(QubitIdx idx, bool val = true) noexcept;
    void SetZ(QubitIdx idx, bool val = true) noexcept;
    void SwapXZ(QubitIdx idx) noexcept;

    [[nodiscard]] bool MeasurementTargets(QubitIdx pivot,
                                          std::span<std::uint64_t> targets) const noexcept;
    [[nodiscard]] bool CleanUp(QubitIdx pivot, std::span<const std::uint64_t> targets,
                               bool apply_s_dagger) noexcept;

    [[nodiscard]] std::int64_t XorWithPhase(const PauliStringView& other,
                                            std::size_t num_qubits) noexcept;

   private:
    std::span<std::uint64_t> xs_;
    std::span<std::uint64_t> zs_;
};

}  // namespace stabaal::detail

#endif
