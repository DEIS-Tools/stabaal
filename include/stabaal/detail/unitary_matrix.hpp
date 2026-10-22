#ifndef STABAAL_DETAIL_UNITARY_MATRIX_HPP
#define STABAAL_DETAIL_UNITARY_MATRIX_HPP

#include <array>
#include <cstddef>
#include <functional>
#include <span>

#include "complex_integer.hpp"

namespace stabaal::detail {

template <std::size_t O>
struct UMat {
    using Row = std::array<CInt, O>;

    static constexpr std::size_t kOrder{O};

    std::array<Row, O> mat{};

    template <typename Op>
    constexpr UMat elementwise(const UMat &o, Op op) const {
        UMat out{};
        for (std::size_t i{0}; i < O; ++i) {
            for (std::size_t j{0}; j < O; ++j) out[i][j] = op(mat[i][j], o[i][j]);
        }

        return out;
    }

    constexpr UMat operator+(const UMat &o) const { return elementwise(o, std::plus{}); }
    constexpr UMat &operator+=(const UMat &o) { return *this = *this + o; }

    constexpr UMat operator-() const { return *this * CInt{-1, 0}; }
    constexpr UMat operator-(const UMat &o) const { return elementwise(o, std::minus{}); }
    constexpr UMat &operator-=(const UMat &o) { return *this = *this - o; }

    constexpr UMat operator*(const UMat &o) const {
        UMat out{};
        for (std::size_t i{0}; i < O; ++i) {
            for (std::size_t j{0}; j < O; ++j) {
                for (std::size_t k{0}; k < O; ++k) out[i][j] += mat[i][k] * o[k][j];
            }
        }

        return out;
    }
    constexpr UMat &operator*=(const UMat &o) { return *this = *this * o; }

    constexpr UMat operator*(const CInt &o) const {
        UMat out{};
        for (std::size_t i{0}; i < O; ++i) {
            for (std::size_t j{0}; j < O; ++j) out[i][j] = mat[i][j] * o;
        }

        return out;
    }
    constexpr UMat &operator*=(const CInt &o) { return *this = *this * o; }

    constexpr Row &operator[](const std::size_t i) { return mat[i]; }
    constexpr const Row &operator[](const std::size_t i) const { return mat[i]; }

    constexpr UMat adjoint() const {
        UMat out{};
        for (std::size_t i{0}; i < O; ++i) {
            for (std::size_t j{0}; j < O; ++j) out[i][j] = mat[j][i].conjugate();
        }

        return out;
    }

    constexpr CInt trace() const {
        CInt sum{};
        for (std::size_t i{0}; i < O; ++i) sum += mat[i][i];
        return sum;
    }

    constexpr bool operator==(const UMat &o) const = default;
};

template <std::size_t N, std::size_t M>
constexpr UMat<N * M> kronecker(const UMat<N> &a, const UMat<M> &b) {
    UMat<N * M> out{};
    for (std::size_t i{0}; i < N; ++i) {
        for (std::size_t j{0}; j < N; ++j) {
            for (std::size_t k{0}; k < M; ++k) {
                for (std::size_t l{0}; l < M; ++l) out[i * M + k][j * M + l] = a[i][j] * b[k][l];
            }
        }
    }

    return out;
}

// Kronecker product of N single-qubit matrices, mats[0] being the most significant factor.
template <std::size_t N>
constexpr UMat<std::size_t{1} << N> kronecker(const std::span<const UMat<2>, N> mats) {
    if constexpr (N == 1) {
        return mats[0];
    } else {
        constexpr auto k_mid{N / 2};
        return kronecker(kronecker<k_mid>(mats.template first<k_mid>()),
                         kronecker<N - k_mid>(mats.template last<N - k_mid>()));
    }
}

}  // namespace stabaal::detail

#endif
