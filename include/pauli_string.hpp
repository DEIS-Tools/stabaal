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
    explicit PauliString(std::size_t num_qubits);

    PauliString(const PauliString& other, std::size_t num_qubits);
    PauliString(const PauliString& other) = delete;
    PauliString& operator=(const PauliString& other) = delete;
    PauliString(PauliString&& other) noexcept = default;
    PauliString& operator=(PauliString&& other) noexcept = default;

    ~PauliString() noexcept = default;

    void AddI(std::size_t old_num_qubits);

    static std::size_t AllocatedBytes(std::size_t num_qubits) noexcept;
    bool AnyX(std::size_t num_qubits) const noexcept;

    bool GetX(QubitIdx idx) const noexcept;
    bool GetZ(QubitIdx idx) const noexcept;

    void SetX(QubitIdx idx, bool val = true) noexcept;
    void SetZ(QubitIdx idx, bool val = true) noexcept;
    void SwapXZ(QubitIdx idx) noexcept;

    std::int64_t XorWithPhase(const PauliString& other, std::size_t num_qubits) noexcept;

   private:
    DynamicBitset xz_;
};

}  // namespace stabaal

#endif
