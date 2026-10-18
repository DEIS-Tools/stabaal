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
    bool AnyX() const noexcept;

    bool GetX(QubitIdx idx) const;
    bool GetZ(QubitIdx idx) const;

    void SetX(QubitIdx idx, bool val = true);
    void SetZ(QubitIdx idx, bool val = true);
    void SwapXZ(QubitIdx idx);

    PauliString& operator^=(const PauliString& other) noexcept;
    std::int64_t PhaseExponent(const PauliString& other) const noexcept;

   private:
    DynamicBitset x_;
    DynamicBitset z_;
};

}  // namespace stabaal

#endif
