#ifndef STABAAL_DETAIL_INVERSE_TABLEAU_HPP
#define STABAAL_DETAIL_INVERSE_TABLEAU_HPP

#include <array>
#include <bit>
#include <cassert>
#include <concepts>
#include <cstddef>
#include <cstdint>
#include <stdexcept>
#include <vector>

#include "basis_table.hpp"
#include "gate_effect.hpp"
#include "pauli_string.hpp"
#include "types.hpp"
#include "unitary_matrix.hpp"

namespace stabaal::detail {

class InverseTableau {
   public:
    enum class Basis : std::uint8_t {
        kComputational,
        kHadamard,
        kImaginary,
    };

    explicit InverseTableau(QubitIdx num_qubits);

    InverseTableau(const InverseTableau& other);
    InverseTableau& operator=(const InverseTableau& other);
    InverseTableau(InverseTableau&& other) noexcept = default;
    InverseTableau& operator=(InverseTableau&& other) noexcept = default;
    ~InverseTableau() noexcept = default;

    void add();

    bool measureX(QubitIdx q);
    bool measureZ(QubitIdx q);

    void applyH(QubitIdx q);
    void applyS(QubitIdx q);
    void applyX(QubitIdx q);
    void applyZ(QubitIdx q);
    void applyCx(QubitIdx c, QubitIdx t);

    template <UMat M, std::convertible_to<QubitIdx>... Qubits>
    void apply(const Qubits... qubits) {
        constexpr auto k_len{sizeof...(Qubits)};
        static_assert(k_len == static_cast<std::size_t>(std::countr_zero(M.kOrder)),
                      "number of qubits must match the gate size");
        static constexpr auto k_effect{findGateEffect<k_len>(M)};

        const std::array<QubitIdx, k_len> qs{static_cast<QubitIdx>(qubits)...};
        for (const auto q : qs) {
            if (q >= num_qubits_) [[unlikely]] {
                throw std::out_of_range{"Qubit index out of range"};
            }
        }

        static constexpr auto k_copy{findRowsToCopy(k_effect)};
        std::array<PauliString, 2 * k_len> scratch{};
        GateRows<k_len> old{};
        for (std::size_t i{0}; i < k_len; ++i) {
            old.gens[2 * i] = &hadamard_.gens[qs[i]];
            old.signs[2 * i] = hadamard_.signs.get(qs[i]);
            old.gens[2 * i + 1] = &computational_.gens[qs[i]];
            old.signs[2 * i + 1] = computational_.signs.get(qs[i]);
        }

        // copy only the rows which are overwritten before all its uses are done
        for (std::size_t r{0}; r < 2 * k_len; ++r) {
            if (!k_copy[r]) continue;
            scratch[r] = PauliString{*old.gens[r], num_qubits_};
            old.gens[r] = &scratch[r];
        }

        // Written in the same order as the rows are numbered, which findRowsToCopy relies on.
        for (std::size_t i{0}; i < k_len; ++i) {
            hadamard_.signs.set(
                qs[i], combineRows(k_effect.hadamard[i], old, hadamard_.gens[qs[i]], num_qubits_));

            computational_.signs.set(qs[i], combineRows(k_effect.computational[i], old,
                                                        computational_.gens[qs[i]], num_qubits_));
        }
    }

// NOLINTNEXTLINE(cppcoreguidelines-macro-usage)
#define STABAAL_GATE_METHOD(name, mat) \
    void name(const auto... qubits) { apply<mat>(qubits...); }
    STABAAL_GATE_LIST(STABAAL_GATE_METHOD)
#undef STABAAL_GATE_METHOD

    void print() const;

    std::size_t memoryUsage() const noexcept;

   private:
    template <std::size_t N>
    struct GateRows {
        std::array<const PauliString*, 2 * N> gens;
        std::array<bool, 2 * N> signs;
    };

    BasisTable computational_;
    BasisTable hadamard_;
    QubitIdx num_qubits_;

    void cleanUp(QubitIdx pivot, const PauliString& gen) noexcept;
    void collapse(QubitIdx pivot, bool p) noexcept;
    bool randomMeasurement(QubitIdx q, Basis b);

    // Writes the product of old rows selected by effect into out and returns its sign.
    template <std::size_t N>
    static bool combineRows(const GenEffect<N>& effect, const GateRows<N>& old, PauliString& out,
                            std::size_t num_qubits) noexcept;
};

template <std::size_t N>
bool InverseTableau::combineRows(const GenEffect<N>& effect, const GateRows<N>& old,
                                 PauliString& out, const std::size_t num_qubits) noexcept {
    std::int64_t exponent{effect.r ? 2 : 0};
    bool first{true};

    const auto multiply{[&](const std::size_t r) {
        if (first) {
            out.copyFrom(*old.gens[r], num_qubits);
            first = false;
        } else {
            exponent += out.xorWithPhase(*old.gens[r], num_qubits);
        }
        exponent += old.signs[r] ? 2 : 0;
    }};

    for (std::size_t i{0}; i < N; ++i) {
        if (effect.x[i]) multiply(2 * i);
        if (effect.z[i]) multiply(2 * i + 1);
        if (effect.x[i] && effect.z[i]) ++exponent;
    }
    assert(!first);

    const auto phase{static_cast<std::uint64_t>(exponent) & 3U};
    assert(phase % 2 == 0);
    return phase == 2;
}

}  // namespace stabaal::detail

#endif
