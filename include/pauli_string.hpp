#ifndef STABAAL_PAULI_STRING_HPP
#define STABAAL_PAULI_STRING_HPP

#include <boost/dynamic_bitset.hpp>
#include <cstddef>
#include <cstdint>

#include "types.hpp"

namespace stabaal {

class PauliString {
   public:
    PauliString();
    explicit PauliString(QubitIdx num_qubits);

    void AddI() noexcept;

    std::size_t Size() const noexcept;
    bool AnyX() const noexcept;

    bool GetX(QubitIdx idx) const noexcept;
    bool GetZ(QubitIdx idx) const noexcept;

    void SetX(QubitIdx idx, bool val = true) noexcept;
    void SetZ(QubitIdx idx, bool val = true) noexcept;
    void SwapXZ(QubitIdx idx) noexcept;

    void XorWith(const PauliString& other) noexcept;
    std::int64_t PhaseExponent(const PauliString& other) const;

   private:
    boost::dynamic_bitset<> x_;
    boost::dynamic_bitset<> z_;
};

}  // namespace stabaal

#endif
