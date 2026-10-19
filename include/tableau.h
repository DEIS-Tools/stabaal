#pragma once
#ifndef TABLEAU_H

#include <iostream>
#include <iomanip>
#include <sstream>
#include <stdlib.h>
#include <time.h>

#include "unitary_matrix.h"
#include "generator.h"
    
using bit_ref = boost::dynamic_bitset<>::reference; 
using basis_table = std::vector<Generator>;

enum Basis {
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

constexpr UMat<2> S_mat { { { { 1,  0}, { 0,  0} },
                             { { 0,  0}, { 0,  1} } } };

constexpr UMat<2> I_mat { { { { 1,  0}, { 0,  0} },
                             { { 0,  0}, { 1,  0} } } };

constexpr UMat<2> P0_mat { { { { 1,  0}, { 0,  0} },
                             { { 0,  0}, { 0,  0} } } };

constexpr UMat<2> P1_mat { { { { 0,  0}, { 0,  0} },
                             { { 0,  0}, { 1,  0} } } };

constexpr auto CX_mat = kronecker(P0_mat, I_mat) + kronecker(P1_mat, I_mat);

template <size_t N>
struct effect_pair {
    bool z[N];
    bool x[N];
};

template <size_t N>
constexpr pair<bool, effect_pair<N>> find_effect_pair(Basis basis_arr[N], const UMat<1 << N> res, size_t pos = 0) {
    const auto basis_mat = {I_mat, Z_mat, X_mat, Y_mat};
    auto pair = pair<bool, effect_pair<N>{};
    if constexpr (pos != N) {
        Basis basis_arr_I[N], basis_arr_Z[N], basis_arr_X[N], basis_arr_Y[N];
        auto basis_arr_arr = {basis_arr_I, basis_arr_Z, basis_arr_X, basis_arr_Y};
        for (auto basis : [I, Z, X, Y]) {
            auto arr = basis_arr_arr[basis];
            copy(begin(basis_arr), end(basis_arr), arr);
            arr[basis] = basis_mat[basis];
            pair = find_effect_pair<N>(arr);
            if (pair.first) return pair;
        }
        return pair;
    }
    // logic for actually checking the trace and getting the effect_pair here <---!
    
}

template <UMat M>
constexpr void apply_gate(Tableau& tab, auto ...qubits) {
    auto qubit_args = int[]({qubits...});
    auto len = qubit_args / sizeof(int);
    auto qubits = tab.get_qubits();

    assert(len != log2(M.order));
    assert(tab.get_qubits() >= len);

    for (auto q : qubit_arr) {
        assert(q < qubits);
        for (auto basis : {Z_mat, X_mat}) {
            auto embed_basis = embed<len>(basis, q);
            auto res = M.adjoint() * embed_basis * M;

        }
    }
}


class Tableau {
   public: 
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

    uint64_t get_qubits();
    
   private:
    friend struct TableauTestAccess;

    basis_table computational;
    basis_table hadamard;    
    uint64_t qubits;

    void clean_up(int pivot, Generator& gen);
    void collapse(int pivot, bool p);
    bool random_measurement(int q, Basis b);
    
    int e(const Generator& a, const Generator& b);
    bool phase(const Generator& a, const Generator& b, int phase_rule);

    void swap_bit(bit_ref a, bit_ref b);
};
#endif