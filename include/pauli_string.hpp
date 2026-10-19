#ifndef STABAAL_PAULI_STRING_HPP
#define STABAAL_PAULI_STRING_HPP

#include <cstddef>
#include <cstdint>

#include "dynamic_bitset.hpp"
#include "types.hpp"

namespace stabaal {
class PauliString {
   public:
    PauliString() noexcept;
    explicit PauliString(QubitIdx num_qubits);

    void AddI();

    std::size_t Size() const noexcept;
    std::size_t AllocatedBytes() const noexcept;
    bool AnyX() const noexcept;

    bool GetX(QubitIdx idx) const noexcept;
    bool GetZ(QubitIdx idx) const noexcept;

    void SetX(QubitIdx idx, bool val = true) noexcept;
    void SetZ(QubitIdx idx, bool val = true) noexcept;
    void SwapXZ(QubitIdx idx) noexcept;

    std::int64_t XorWithPhase(const PauliString& other) noexcept;

   private:
    DynamicBitset x_;
    DynamicBitset z_;
};

}  // namespace stabaal

#endif
