#ifndef STABAAL_DETAIL_PAULI_STRING_HPP
#define STABAAL_DETAIL_PAULI_STRING_HPP

#include <cstddef>
#include <cstdint>

#include "dynamic_bitset.hpp"
#include "types.hpp"

namespace stabaal::detail {
class PauliString {
   public:
    PauliString() noexcept;
    explicit PauliString(std::size_t num_qubits);

    PauliString(const PauliString& other, std::size_t num_qubits);
    PauliString(const PauliString& other) = delete;
    PauliString& operator=(const PauliString& other) = delete;
    PauliString(PauliString&& other) noexcept = default;
    PauliString& operator=(PauliString&& other) noexcept = default;

    ~PauliString() noexcept = default;

    void addI(std::size_t old_num_qubits);

    static std::size_t allocatedBytes(std::size_t num_qubits) noexcept;
    bool anyX(std::size_t num_qubits) const noexcept;

    bool getX(QubitIdx idx) const noexcept;
    bool getZ(QubitIdx idx) const noexcept;

    void setX(QubitIdx idx, bool val = true) noexcept;
    void setZ(QubitIdx idx, bool val = true) noexcept;
    void swapXZ(QubitIdx idx) noexcept;

    void copyFrom(const PauliString& other, std::size_t num_qubits) noexcept;

    std::int64_t xorWithPhase(const PauliString& other, std::size_t num_qubits) noexcept;

   private:
    DynamicBitset xz_;
};

}  // namespace stabaal::detail

#endif
