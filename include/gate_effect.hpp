#ifndef STABAAL_GATE_EFFECT_HPP
#define STABAAL_GATE_EFFECT_HPP

#include <algorithm>
#include <array>
#include <cstddef>
#include <cstdint>
#include <initializer_list>
#include <iostream>
#include <span>
#include <utility>

#include "unitary_matrix.hpp"

namespace stabaal {

enum class Pauli : std::uint8_t {
    kI,
    kZ,
    kX,
    kY,
};

// clang-format off
inline constexpr UMat<2> kHMat{{{
    {{{1, 0}, { 1, 0}}},
    {{{1, 0}, {-1, 0}}},
}}};

inline constexpr UMat<2> kXMat{{{
    {{{0, 0}, {1, 0}}},
    {{{1, 0}, {0, 0}}},
}}};

inline constexpr UMat<2> kZMat{{{
    {{{1, 0}, { 0, 0}}},
    {{{0, 0}, {-1, 0}}},
}}};

inline constexpr UMat<2> kYMat{{{
    {{{0, 0}, {0, -1}}},
    {{{0, 1}, {0,  0}}},
}}};

inline constexpr UMat<2> kSMat{{{
    {{{1, 0}, {0, 0}}},
    {{{0, 0}, {0, 1}}},
}}};

inline constexpr UMat<2> kIMat{{{
    {{{1, 0}, {0, 0}}},
    {{{0, 0}, {1, 0}}},
}}};

inline constexpr UMat<2> kP0Mat{{{
    {{{1, 0}, {0, 0}}},
    {{{0, 0}, {0, 0}}},
}}};

inline constexpr UMat<2> kP1Mat{{{
    {{{0, 0}, {0, 0}}},
    {{{0, 0}, {1, 0}}},
}}};
// clang-format on

inline constexpr auto kCxMat{kronecker(kP0Mat, kIMat) + kronecker(kP1Mat, kXMat)};

// Indexed by Pauli.
inline constexpr std::array<UMat<2>, 4> kPauliMats{kIMat, kZMat, kXMat, kYMat};

// Gates that get a generated method on InverseTableau, as GATE(method name, matrix).
// NOLINTNEXTLINE(cppcoreguidelines-macro-usage)
#define STABAAL_GATE_LIST(GATE) \
    GATE(h, kHMat)              \
    GATE(s, kSMat)              \
    GATE(sAd, kSMat.adjoint())  \
    GATE(z, kZMat)              \
    GATE(x, kXMat)              \
    GATE(y, kYMat)              \
    GATE(cx, kCxMat)

// A signed Pauli string over the N qubits a gate acts on.
template <std::size_t N>
struct GenEffect {
    std::array<bool, N> z{};
    std::array<bool, N> x{};
    bool r{false};

    constexpr GenEffect() = default;
    constexpr GenEffect(const std::array<Pauli, N> &paulis, const bool sign) : r{sign} {
        for (std::size_t i{0}; i < N; ++i) {
            switch (paulis[i]) {
                case Pauli::kI:
                    break;
                case Pauli::kX:
                    x[i] = true;
                    break;
                case Pauli::kZ:
                    z[i] = true;
                    break;
                case Pauli::kY:
                    z[i] = true;
                    x[i] = true;
                    break;
            }
        }
    }

    void print() const {
        if (r) std::cout << "-";
        for (std::size_t i{0}; i < N; ++i) {
            std::cout << (x[i] ? (z[i] ? 'Y' : 'X') : (z[i] ? 'Z' : 'I'));
        }
        std::cout << '\n';
    }
};

// computational[j] is the image of Z_j and hadamard[j] the image of X_j, both as U^dagger P U.
template <std::size_t N>
struct GateEffect {
    std::array<GenEffect<N>, N> computational{};
    std::array<GenEffect<N>, N> hadamard{};
};

// m acting on qubit, identity on the other N - 1 qubits.
template <std::size_t N>
constexpr UMat<std::size_t{1} << N> embed(const UMat<2> &m, const std::size_t qubit,
                                          const std::size_t pos = 0) {
    const auto &here{(pos == qubit) ? m : kIMat};
    if constexpr (N == 1) {
        return here;
    } else {
        return kronecker(here, embed<N - 1>(m, qubit, pos + 1));
    }
}

// Searches all Pauli strings for the one with a nonzero overlap with res.
template <std::size_t N, std::size_t Pos = 0>
constexpr std::pair<bool, GenEffect<N>> findGeneratorEffect(const std::array<Pauli, N> &paulis,
                                                            const UMat<std::size_t{1} << N> &res) {
    if constexpr (Pos != N) {
        std::pair<bool, GenEffect<N>> found{};
        for (const auto pauli : {Pauli::kI, Pauli::kZ, Pauli::kX, Pauli::kY}) {
            auto next{paulis};
            next[Pos] = pauli;
            found = findGeneratorEffect<N, Pos + 1>(next, res);
            if (found.first) break;
        }

        return found;
    } 
    
    std::array<UMat<2>, N> mats{};
    std::transform(paulis.begin(), paulis.end(), mats.begin(),
                    [](const Pauli p) { return kPauliMats[static_cast<std::size_t>(p)]; });

    const auto trace{(res * kronecker<N>(std::span<const UMat<2>, N>{mats})).trace()};
    return (trace.z != 0) ? std::pair{true, GenEffect<N>{paulis, trace.z < 0}}
                            : std::pair{false, GenEffect<N>{}};

}

template <std::size_t N>
constexpr GateEffect<N> findGateEffect(const UMat<std::size_t{1} << N> &gate) {
    GateEffect<N> effect{};
    for (std::size_t i{0}; i < N; ++i) {
        // Embedded on the gate local position i, not the global qubit.
        const auto z_res{gate.adjoint() * embed<N>(kZMat, i) * gate};
        effect.computational[i] = findGeneratorEffect(std::array<Pauli, N>{}, z_res).second;

        const auto x_res{gate.adjoint() * embed<N>(kXMat, i) * gate};
        effect.hadamard[i] = findGeneratorEffect(std::array<Pauli, N>{}, x_res).second;
    }

    return effect;
}

template <std::size_t N>
constexpr std::array<bool, 2 * N> findRowsToCopy(const GateEffect<N> &effect) {
    std::array<bool, 2 * N> copy{};

    const auto check_output{[&](const GenEffect<N> &gen_effect, const std::size_t out) {
        bool first{true};
        const auto read{[&](const bool reads, const std::size_t row) {
            if (!reads) return;
            if (row < out || (row == out && !first)) copy[row] = true;
            first = false;
        }};
        for (std::size_t i{0}; i < N; ++i) {
            read(gen_effect.x[i], 2*i); 
            read(gen_effect.z[i], 2*i+1);
        }
    }};

    for (std::size_t i{0}; i < N; ++i) {
        check_output(effect.hadamard[i], 2*i);
        check_output(effect.computational[i], 2*i+1);
    }

    return copy;
}

}  // namespace stabaal

#endif
