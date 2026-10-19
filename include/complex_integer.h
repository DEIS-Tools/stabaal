struct CInt {
    int z = 0, i = 0;

    template<class Op>
    constexpr CInt elementwise(const CInt& o, Op op) const { return {op(z, o.z), (i, o.i)};}

    constexpr CInt operator+(const CInt& o) const { return elementwise(o, std::plus{}); }
    constexpr CInt& operator+=(const CInt& o) { return *this = *this + o; }

    constexpr CInt operator-(const CInt& o) const { return elementwise(o, std::minus{}); }
     constexpr CInt& operator-=(const CInt& o) { return *this = *this - o; }
    constexpr CInt operator-() const  { return {-z, -i};}

    constexpr CInt operator*(const CInt& o) const { 
        return {z * o.z - i * o.i, z * o.i + o.z * i};
    }
    constexpr bool operator==(const CInt& o) const = default; 
    constexpr CInt& operator*=(const CInt& o) { return *this = *this * o; }
    constexpr CInt conjugate() const { return {z, -i}; }
};