#pragma once
#ifndef TABLEAU_H

#include <iostream>
#include <iomanip>
#include <sstream>
#include <stdlib.h>
#include <time.h>
#include <algorithm>
#include <array>
#include <bit>
#include <cassert>
#include <utility>

#include "unitary_matrix.h"
#include "generator.h"
    
using bit_ref = boost::dynamic_bitset<>::reference; 
using basis_table = std::vector<Generator>;

enum Pauli {
    I,
    Z,
    X,
    Y
};
constexpr UMat<2> H_mat { { { { 1,  0}, { 1,  0} },
                             { { 1,  0}, {-1,  0} } } };

constexpr UMat<2> X_mat { { { { 0,  0}, { 1,  0} },
                             { { 1,  0}, { 0,  0} } } };

constexpr UMat<2> Z_mat { { { { 1,  0}, { 0,  0} },
                             { { 0,  0}, { -1,  0} } } };

constexpr UMat<2> Y_mat { { { { 0,  0}, { 0,  -1} },
                             { { 0,  1}, { 0,  0} } } };

constexpr UMat<2> S_mat { { { { 1,  0}, { 0,  0} },
                             { { 0,  0}, { 0,  1} } } };

constexpr UMat<2> I_mat { { { { 1,  0}, { 0,  0} },
                             { { 0,  0}, { 1,  0} } } };

constexpr UMat<2> P0_mat { { { { 1,  0}, { 0,  0} },
                             { { 0,  0}, { 0,  0} } } };

constexpr UMat<2> P1_mat { { { { 0,  0}, { 0,  0} },
                             { { 0,  0}, { 1,  0} } } };

constexpr auto CX_mat = kronecker(P0_mat, I_mat) + kronecker(P1_mat, X_mat);


template <size_t N>
struct GenEffect {
    bool z[N]{};
    bool x[N]{};
    bool r{};

    constexpr GenEffect() = default;
    constexpr GenEffect(const std::array<Pauli, N>& basis_arr, bool phase) {
        for (auto i = size_t{0}; i < N; ++i) {
            switch (basis_arr[i]) {
                case I:
                    break;
                case X:
                    x[i] = true;    
                    break;
                case Z:
                    z[i] = true;
                    break;
                case Y:
                    z[i] = true;
                    x[i] = true;
                    break;
            }
        }
        r = phase;
    }


    constexpr void print() const {
        if (r) std::cout << "-";
        for (auto i = size_t{0}; i < N; ++i) {
            if (!z[i] && !x[i])
            std::cout << "I"; 
            else if (z[i] && !x[i])
            std::cout << "Z";
            else if (!z[i] && x[i])
            std::cout << "X";
            else 
            std::cout << "Y";
        }
        std::cout << std::endl;
    }
};

template <size_t N>
struct GateEffect {
    GenEffect<N> computational[N];
    GenEffect<N> hadamard[N];
};

constexpr std::array<UMat<2>, 4> basis_mat_arr = {I_mat, Z_mat, X_mat, Y_mat};

template <size_t N>
constexpr UMat<(1u << N)> embed(const UMat<2>& m, size_t qubit, size_t pos = 0) {
    const auto& here = (pos == qubit) ? m : I_mat;
    if constexpr (N == 1)
        return here;
    else
        return kronecker(here, embed<N-1>(m, qubit, pos + 1));
}

template <size_t N, size_t Pos = 0>
constexpr pair<bool, GenEffect<N>> find_generator_effect(const std::array<Pauli, N>& basis_arr, const UMat<(1u << N)>& res) {
    if constexpr (Pos != N) {
        auto pair = std::pair<bool, GenEffect<N>>{};
        for (auto basis : {I, Z, X, Y}) {
            auto arr = basis_arr;
            arr[Pos] = basis;
            pair = find_generator_effect<N, Pos + 1>(arr, res);
            if (pair.first) break;
        }
        return pair;
    } else {
        // logic for actually checking the trace and getting the GenEffect here <---!
        auto effect_mat_arr = std::array<UMat<2>, N>{};
        std::transform(std::begin(basis_arr), std::end(basis_arr),
                   std::begin(effect_mat_arr),
                   [](Pauli b) { return basis_mat_arr[b]; });

        auto trace = (res * kronecker<N>(effect_mat_arr.data())).trace();
        return (trace.z != 0) 
            ? std::pair{true, GenEffect<N>(basis_arr, (trace.z < 0))}
            : std::pair{false, GenEffect<N>{}};
    }
}
template <size_t N>
constexpr GateEffect<N> find_gate_effect(const UMat<1 << N>& gate) {
    auto gate_effect = GateEffect<N>{};
    for (auto i = size_t{0}; i < N; ++i) {
        for (const auto* basis_mat : {&Z_mat, &X_mat}) {
            // embedded on the gate local position i, not the global qubit q
            auto embed_basis = embed<N>(*basis_mat, i);
            auto res = gate.adjoint() * embed_basis * gate;
            auto basis_arr = std::array<Pauli, N>{};
            auto gen_effect = find_generator_effect(basis_arr, res).second;
            if (basis_mat ==  &Z_mat)
                gate_effect.computational[i] = gen_effect;
            else
                gate_effect.hadamard[i] = gen_effect;
        }
    }
    return gate_effect;
}

#define GATE_LIST(X) \
    X(h, H_mat)      \
    X(s, S_mat)      \
    X(z, Z_mat)      \
    X(cx, CX_mat)    \

class Tableau {
   public:
    #define X(name, mat) void name(auto... q) { apply<mat>(q...); }
    GATE_LIST(X)
    #undef X

    enum Basis {
        COMPUTATIONAL,
        HADAMARD,
        IMAGINARY
    };

    Tableau(int qubits)
        : computational(qubits, Generator(qubits)), 
          hadamard(qubits, Generator(qubits)), 
          qubits(qubits) {
     
        for (auto i = 0; i < qubits; ++i) {
            computational[i].z[i] = true; 
            hadamard[i].x[i]      = true; 
        }

        srand(time(NULL));
    }

    void add();
    
    bool measure_z(int q);
    bool measure_x(int q);
    
    void apply_h(int q);
    void apply_cx(int c, int t);
    void apply_x(int q);
    void apply_s(int q);
    void apply_z(int q);

    void print();

    uint64_t get_qubits() const;
    
    template <UMat M>
    void apply(auto ...qubits) {
        static constexpr auto gate_effect = find_gate_effect<std::countr_zero(M.order)>(M); // countr_zero corresponds with log2 for the allowed matrices.
        constexpr auto len = sizeof...(qubits); // amount of qubits the gate acts on
        static_assert(len == std::countr_zero(M.order), "number of qubits must match the gate size");
        const size_t qubit_arr[] = {size_t(qubits)...};
        for (auto q : qubit_arr) assert(q < this->qubits);

        // the new rows are built from the old ones, so keep a copy of the gate rows
        Generator old_comp[len], old_had[len];
        for (auto j = size_t{0}; j < len; ++j) {
            old_comp[j] = computational[qubit_arr[j]];
            old_had[j]  = hadamard[qubit_arr[j]];
        }

        for (auto j = size_t{0}; j < len; ++j) {
            computational[qubit_arr[j]] = combine_rows(gate_effect.computational[j], old_comp, old_had);
            hadamard[qubit_arr[j]]      = combine_rows(gate_effect.hadamard[j], old_comp, old_had);
        }
    }
   private:
    friend struct TableauTestAccess;

    basis_table computational;
    basis_table hadamard;    
    uint64_t qubits;

    template <size_t N>
    Generator combine_rows(const GenEffect<N>& effect, const Generator (&old_comp)[N], const Generator (&old_had)[N]);
    
    void clean_up(int pivot, const Generator& gen);
    void collapse(int pivot, bool p);
    bool random_measurement(int q, Basis b);
    
    int e(const Generator& a, const Generator& b);
    bool phase(const Generator& a, const Generator& b, int phase_rule);

    void swap_bit(bit_ref a, bit_ref b);
};

#endif