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

    PauliString(const PauliString& other);
    PauliString& operator=(const PauliString& other);
    PauliString(PauliString&& other) noexcept;
    PauliString& operator=(PauliString&& other) noexcept;

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
    static std::size_t BitsRequired(std::size_t num_qubits) noexcept;
    static std::size_t GetXIdx(QubitIdx idx) noexcept;
    static std::size_t GetZIdx(QubitIdx idx) noexcept;

    DynamicBitset xz_;
    std::size_t num_qubits_{0};

    static constexpr std::size_t kQubitsPerWord{64};
};

}  // namespace stabaal

#endif
