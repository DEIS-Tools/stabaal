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

    void addI();

    std::size_t size() const noexcept;
    std::size_t allocatedBytes() const noexcept;
    bool anyX() const noexcept;

    bool getX(QubitIdx idx) const noexcept;
    bool getZ(QubitIdx idx) const noexcept;

    void setX(QubitIdx idx, bool val = true) noexcept;
    void setZ(QubitIdx idx, bool val = true) noexcept;
    void swapXZ(QubitIdx idx) noexcept;

    std::int64_t xorWithPhase(const PauliString &other) noexcept;

   private:
    DynamicBitset x_;
    DynamicBitset z_;
};

}  // namespace stabaal

#endif
