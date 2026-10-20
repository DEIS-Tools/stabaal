#ifndef STABAAL_INVERSE_TABLEAU_HPP
#define STABAAL_INVERSE_TABLEAU_HPP

#include <array>
#include <bit>
#include <cassert>
#include <concepts>
#include <cstddef>
#include <cstdint>
#include <iomanip>
#include <iostream>
#include <random>
#include <sstream>
#include <stdexcept>
#include <vector>

#include "basis_table.hpp"
#include "gate_effect.hpp"
#include "pauli_string.hpp"
#include "types.hpp"
#include "unitary_matrix.hpp"

namespace stabaal {

class InverseTableau {
   public:
    enum class Basis : std::uint8_t {
        kComputational,
        kHadamard,
        kImaginary,
    };

    explicit InverseTableau(QubitIdx num_qubits);

    void add();

    bool measureZ(QubitIdx q);
    bool measureX(QubitIdx q);

    void applyH(QubitIdx q);
    void applyCx(QubitIdx c, QubitIdx t);
    void applyX(QubitIdx q);
    void applyS(QubitIdx q);
    void applyZ(QubitIdx q);

    // Applies the Clifford gate M to the given qubits, its effect on the tableau rows is derived
    // from the matrix at compile time.
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

        // The new rows are built from the old ones, so rows that are overwritten before their
        // last read are copied to scratch_ first, the rest are read in place.
        static constexpr auto k_copy{findRowsToCopy(k_effect)};
        if (scratch_.size() < 2 * k_len) scratch_.resize(2 * k_len);
        GateRows<k_len> old{};
        for (std::size_t j{0}; j < k_len; ++j) {
            old.gens[2 * j] = &computational_.gens[qs[j]];
            old.gens[2 * j + 1] = &hadamard_.gens[qs[j]];
            old.signs[2 * j] = computational_.signs.get(qs[j]);
            old.signs[2 * j + 1] = hadamard_.signs.get(qs[j]);
        }
        for (std::size_t r{0}; r < 2 * k_len; ++r) {
            if (!k_copy[r]) continue;
            scratch_[r] = *old.gens[r];
            old.gens[r] = &scratch_[r];
        }

        for (std::size_t j{0}; j < k_len; ++j) {
            computational_.signs.set(
                qs[j], combineRows(k_effect.computational[j], old, computational_.gens[qs[j]]));
            hadamard_.signs.set(qs[j],
                                combineRows(k_effect.hadamard[j], old, hadamard_.gens[qs[j]]));
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
    // The old rows a gate reads, numbered as in findRowsToCopy.
    template <std::size_t N>
    struct GateRows {
        std::array<const PauliString *, 2 * N> gens;
        std::array<bool, 2 * N> signs;
    };

    static std::mt19937 &randomEngine();

    BasisTable computational_;
    BasisTable hadamard_;
    QubitIdx num_qubits_;
    // Copies of old rows for apply, kept between gates to reuse their allocations.
    std::vector<PauliString> scratch_;

    void cleanUp(QubitIdx pivot, const PauliString &gen) noexcept;
    void collapse(QubitIdx pivot, bool p) noexcept;
    bool randomMeasurement(QubitIdx q, Basis b);

    // Writes the product of old rows selected by effect into out and returns its sign.
    template <std::size_t N>
    static bool combineRows(const GenEffect<N> &effect, const GateRows<N> &old, PauliString &out);

    static bool phase(std::int64_t exponent, bool a_sign, bool b_sign,
                      std::uint8_t phase_rule) noexcept;
};

template <std::size_t N>
bool InverseTableau::combineRows(const GenEffect<N> &effect, const GateRows<N> &old,
                                 PauliString &out) {
    std::int64_t exponent{effect.r ? 2 : 0};
    bool first{true};
    // The first row is copied into out, which reuses its allocation and equals multiplying it
    // onto the identity.
    const auto multiply{[&](const std::size_t r) {
        if (first) {
            out = *old.gens[r];
            first = false;
        } else {
            exponent += out.xorWithPhase(*old.gens[r]);
        }
        exponent += old.signs[r] ? 2 : 0;
    }};

    for (std::size_t k{0}; k < N; ++k) {
        if (effect.x[k]) multiply(2 * k + 1);
        if (effect.z[k]) multiply(2 * k);
        // Y = iXZ
        if (effect.x[k] && effect.z[k]) ++exponent;
    }
    // A unitary never maps Z_k or X_k to the identity.
    assert(!first);

    const auto phase{static_cast<std::uint64_t>(exponent) & 3U};
    assert(phase % 2 == 0);
    return phase == 2;
}

}  // namespace stabaal

#endif
