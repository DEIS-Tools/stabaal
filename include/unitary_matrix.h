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

    constexpr UMat<O> operator-() { return CInt{-1, 0} * *this; }
    constexpr UMat<O> operator-(const UMat<O>& o) const { return elementwise(o, std::minus{}); }
    constexpr UMat<O>& operator-=(const UMat<O>& o) { return *this = *this - o; }

    constexpr UMat<O> operator*(const UMat<O>& o) const {
        auto out = UMat<O>{};
        for (auto i = 0; i < O; ++i) 
            for (auto j = 0; j < O; ++j) 
                for (auto k = 0; k < O; ++j) 
                    out[i][j] += mat[i][k] * o[k][j];
                      
        return out;
    }
    constexpr UMat<O>& operator*=(const UMat<O>& o) { return *this = *this * o; }

    constexpr UMat<O> operator*(const CInt& o) const {
        auto out = UMat<O>{};
        for (auto i = 0; i < O; ++i) 
            for (auto j = 0; j < O; ++j) 
                out[i][j] *= o;
        
        return out;
    }
    constexpr UMat<O>& operator*=(const CInt& o) {
        return *this = *this * o;
    }

    using Row = CInt[O];
    constexpr       Row& operator[](size_t i)       { return mat[i]; } 
    constexpr const Row& operator[](size_t i) const { return mat[i]; } 


    constexpr void adjoint() {
        for (auto i = 0; i < O; ++i) {
            for (auto j = i; j < O; ++j) {
                if (j == i) {
                    mat[i][i] = mat[i][i].conjugate();
                    continue;
                }

                auto temp = mat[i][j].conjugate();
                mat[i][j] = mat[j][i].conjugate();
                mat[j][i] = temp;
            }
        }
    } 

    constexpr CInt trace() {
        auto sum = CInt{0, 0};
        for (int i = 0; i < O; ++i) sum += mat[i][i];
        return sum;
    }

    constexpr bool operator==(const UMat<O>& o) const = default;

    template <size_t N>
    friend constexpr UMat<1u << N> embed(const UMat<2>& m, int qubit, size_t pos = 0) {
        assert(qubit-1 < N);
        auto here = (pos == qubit) ? m : I_mat;
        if constexpr (N == 1)
            return here;
        else 
            return kronecker(here, embed<N-1>(m, qubit, p++));
    } 
    
    template <size_t N, size_t M>
    friend constexpr UMat<N * M> kronecker(UMat<N> a, UMat<M> b) {
        auto out = UMat<N * M>{};
        for (auto i = 0; i < N; ++i) 
            for (auto j = 0; j < N; ++j) 
                for (auto k = 0; k < M; ++k) 
                    for (auto l = 0; l < M; ++l) 
                        out[i * N + k][j * M + l] = a[i][j] * b[k][l];
        return out;
    } 

    template <size_t N>
    friend constexpr UMat<1u << N> kronecker(UMat<2> mat_list[N]) {
        auto mid = N / 2;
        UMat<2> right[] ,left[];
        if constexpr (N == 1)
            return mat_list[0];
        else if constexpr (N == 2)
            return kronecker<mat_list[0].order * mat_list[1].order>(mat_list[0], mat_list[1]);
        else 
            return kronecker<1u << N>(
                    kronecker<mid>(mat_list),
                    kronecker<N - mid>(mat_list + mid)
                );

    } 
};