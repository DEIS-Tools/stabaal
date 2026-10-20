#pragma once

#include <cstddef>

#include "complex_integer.h"

using namespace std;

template <size_t O>
struct UMat {
   public:
    CInt mat[O][O];
    const static size_t order = O;

    template<class Op>
    constexpr UMat<O> elementwise(const UMat<O>& o, Op op) const {
        auto out = UMat<O>{};
        for (auto i = 0; i < O; ++i) 
            for (auto j = 0; j < O; ++j) 
                out[i][j] = op(mat[i][j], o[i][j]);

        return out;
    }

    constexpr UMat<O> operator+(const UMat<O>& o) const { return elementwise(o, std::plus{}); }
    constexpr UMat<O>& operator+=(const UMat<O>& o) { return *this = *this + o; }

    constexpr UMat<O> operator-() const { return *this * CInt{-1, 0}; }
    constexpr UMat<O> operator-(const UMat<O>& o) const { return elementwise(o, std::minus{}); }
    constexpr UMat<O>& operator-=(const UMat<O>& o) { return *this = *this - o; }

    constexpr UMat<O> operator*(const UMat<O>& o) const {
        auto out = UMat<O>{};
        for (auto i = 0; i < O; ++i) 
            for (auto j = 0; j < O; ++j) 
                for (auto k = 0; k < O; ++k)
                    out[i][j] += mat[i][k] * o[k][j];
                      
        return out;
    }
    constexpr UMat<O>& operator*=(const UMat<O>& o) { return *this = *this * o; }

    constexpr UMat<O> operator*(const CInt& o) const {
        auto out = UMat<O>{};
        for (auto i = 0; i < O; ++i) 
            for (auto j = 0; j < O; ++j) 
                out[i][j] = mat[i][j] * o;
        
        return out;
    }
    constexpr UMat<O>& operator*=(const CInt& o) {
        return *this = *this * o;
    }

    using Row = CInt[O];
    constexpr       Row& operator[](size_t i)       { return mat[i]; } 
    constexpr const Row& operator[](size_t i) const { return mat[i]; } 


    constexpr UMat<O> adjoint() const {
        auto out = UMat<O>{};
        for (auto i = 0; i < O; ++i)
            for (auto j = 0; j < O; ++j)
                out[i][j] = mat[j][i].conjugate();

        return out;
    }

    constexpr CInt trace() const {
        auto sum = CInt{0, 0};
        for (int i = 0; i < O; ++i) sum += mat[i][i];
        return sum;
    }

    constexpr bool operator==(const UMat<O>& o) const = default;

};

template <size_t N, size_t M>
constexpr UMat<N * M> kronecker(const UMat<N>& a, const UMat<M>& b) {
    auto out = UMat<N * M>{};
    for (auto i = 0; i < N; ++i)
        for (auto j = 0; j < N; ++j)
            for (auto k = 0; k < M; ++k)
                for (auto l = 0; l < M; ++l)
                    out[i * M + k][j * M + l] = a[i][j] * b[k][l];
    return out;
}

template <size_t N>
constexpr UMat<1u << N> kronecker(const UMat<2> mat_list[N]) {
    constexpr auto mid = N / 2;
    if constexpr (N == 1)
        return mat_list[0];
    else
        return kronecker(
                kronecker<mid>(mat_list),
                kronecker<N - mid>(mat_list + mid)
            );
}
