#ifndef STABAAL_DETAIL_COMPLEX_INTEGER_HPP
#define STABAAL_DETAIL_COMPLEX_INTEGER_HPP

#include <functional>

namespace stabaal::detail {

// Gaussian integer z + i * i, enough to represent Clifford gates up to normalisation.
struct CInt {
    int z{0};
    int i{0};

    template <typename Op>
    constexpr CInt elementwise(const CInt &o, Op op) const {
        return {op(z, o.z), op(i, o.i)};
    }

    constexpr CInt operator+(const CInt &o) const { return elementwise(o, std::plus{}); }
    constexpr CInt &operator+=(const CInt &o) { return *this = *this + o; }

    constexpr CInt operator-(const CInt &o) const { return elementwise(o, std::minus{}); }
    constexpr CInt &operator-=(const CInt &o) { return *this = *this - o; }
    constexpr CInt operator-() const { return {-z, -i}; }

    constexpr CInt operator*(const CInt &o) const { return {z * o.z - i * o.i, z * o.i + o.z * i}; }
    constexpr CInt &operator*=(const CInt &o) { return *this = *this * o; }

    constexpr bool operator==(const CInt &o) const = default;

    constexpr CInt conjugate() const { return {z, -i}; }
};

}  // namespace stabaal::detail

#endif
