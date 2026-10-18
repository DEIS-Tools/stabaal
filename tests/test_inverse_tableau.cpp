#include <array>
#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators.hpp>
#include <cmath>
#include <complex>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <random>
#include <set>
#include <stdexcept>
#include <string>
#include <vector>

#include "basis_table.hpp"
#include "inverse_tableau.hpp"
#include "types.hpp"

namespace stabaal {

// Builds unsigned operators from a Pauli string, e.g. Pauli("XIZY") = X⊗I⊗Z⊗Y.
// Character i is the Pauli acting on qubit i.
static PauliString Pauli(const std::string& s) {
    PauliString g{static_cast<QubitIdx>(s.size())};
    for (std::size_t i{0}; i < s.size(); ++i) {
        switch (s[i]) {
            case 'I':
                break;
            case 'X':
                g.SetX(i);
                break;
            case 'Z':
                g.SetZ(i);
                break;
            case 'Y':
                g.SetX(i);
                g.SetZ(i);
                break;
            default:
                throw std::invalid_argument{"bad pauli char"};
        }
    }
    return g;
}

// Renders a generator as a signed Pauli string, e.g. "-XIZY".
static std::string Str(const PauliString& g, const bool sign = false) {
    std::string s{sign ? "-" : "+"};
    const auto n{g.Size()};
    for (std::size_t i{0}; i < n; ++i) {
        const bool x{g.GetX(i)};
        const bool z{g.GetZ(i)};
        s += x ? (z ? 'Y' : 'X') : (z ? 'Z' : 'I');
    }
    return s;
}

// Redirects std::cout for the lifetime of the object.
struct CoutCapture {
    std::ostringstream ss{};
    std::streambuf* old;
    CoutCapture() : old{std::cout.rdbuf(ss.rdbuf())} {}
    ~CoutCapture() { std::cout.rdbuf(old); }
};

// 2x2 Pauli matrices used as an independent reference for e().
using Mat = std::array<std::complex<double>, 4>;
static Mat PauliMatrix(const char p) {
    using Complex = std::complex<double>;
    switch (p) {
        case 'X':
            return {Complex{0}, Complex{1}, Complex{1}, Complex{0}};
        case 'Y':
            return {Complex{0}, Complex{0, -1}, Complex{0, 1}, Complex{0}};
        case 'Z':
            return {Complex{1}, Complex{0}, Complex{0}, Complex{-1}};
        default:
            return {Complex{1}, Complex{0}, Complex{0}, Complex{1}};
    }
}
static Mat Multiply(const Mat& a, const Mat& b) {
    return {a[0] * b[0] + a[1] * b[2], a[0] * b[1] + a[1] * b[3], a[2] * b[0] + a[3] * b[2],
            a[2] * b[1] + a[3] * b[3]};
}

struct StateVector {
    std::vector<std::complex<double>> amplitudes;

    explicit StateVector(const QubitIdx n) : amplitudes(std::size_t{1} << n) { amplitudes[0] = 1; }

    void Add() { amplitudes.resize(amplitudes.size() * 2); }

    void Apply(const std::uint8_t gate, const QubitIdx q, const QubitIdx target) {
        const auto mask{std::size_t{1} << q};
        for (std::size_t i{0}; i < amplitudes.size(); ++i) {
            if ((i & mask) != 0) continue;
            const auto j{i | mask};
            switch (gate) {
                case 0: {
                    const auto a{amplitudes[i]};
                    const auto b{amplitudes[j]};
                    amplitudes[i] = (a + b) / std::sqrt(2.0);
                    amplitudes[j] = (a - b) / std::sqrt(2.0);
                    break;
                }
                case 1:
                    std::swap(amplitudes[i], amplitudes[j]);
                    break;
                case 2:
                    amplitudes[j] *= std::complex<double>{0, 1};
                    break;
                case 3:
                    amplitudes[j] = -amplitudes[j];
                    break;
                case 4:
                    if ((j & (std::size_t{1} << target)) == 0) {
                        std::swap(amplitudes[j], amplitudes[j | (std::size_t{1} << target)]);
                    }
                    break;
            }
        }
    }

    void MeasureZ(const QubitIdx q, const bool result) {
        const auto mask{std::size_t{1} << q};
        double probability{0};
        for (std::size_t i{0}; i < amplitudes.size(); ++i) {
            if (((i & mask) != 0) == result) probability += std::norm(amplitudes[i]);
        }
        REQUIRE(probability > 1e-10);
        for (std::size_t i{0}; i < amplitudes.size(); ++i) {
            amplitudes[i] = (((i & mask) != 0) == result) ? amplitudes[i] / std::sqrt(probability)
                                                          : std::complex<double>{};
        }
    }
};

TEST_CASE("InverseTableau gates and measurements agree with state-vector simulation",
          "[inverse_tableau]") {
    const QubitIdx initial_qubits{GENERATE(QubitIdx{1}, QubitIdx{2}, QubitIdx{3}, QubitIdx{4})};
    std::mt19937 choices{1234};
    for (std::uint32_t trial{0}; trial < 16; ++trial) {
        auto n{initial_qubits};
        InverseTableau tableau{n};
        StateVector reference{n};
        for (std::uint32_t step{0}; step < 48; ++step) {
            if (step == 24) {
                tableau.Add();
                reference.Add();
                ++n;
            }
            const auto gate{static_cast<std::uint8_t>(choices() % (n == 1 ? 4 : 5))};
            const QubitIdx q{static_cast<QubitIdx>(choices() % n)};
            const QubitIdx target{
                static_cast<QubitIdx>((q + 1 + choices() % (n == 1 ? 1 : n - 1)) % n)};
            CAPTURE(trial, step, gate, q, target);
            switch (gate) {
                case 0:
                    tableau.ApplyH(q);
                    break;
                case 1:
                    tableau.ApplyX(q);
                    break;
                case 2:
                    tableau.ApplyS(q);
                    break;
                case 3:
                    tableau.ApplyZ(q);
                    break;
                case 4:
                    tableau.ApplyCx(q, target);
                    break;
            }
            reference.Apply(gate, q, target);
            if (step % 8 == 7) {
                const bool x_basis{(choices() & 1U) != 0};
                if (x_basis) reference.Apply(0, q, 0);
                const bool result{x_basis ? tableau.MeasureX(q) : tableau.MeasureZ(q)};
                reference.MeasureZ(q, result);
                CHECK((x_basis ? tableau.MeasureX(q) : tableau.MeasureZ(q)) == result);
                if (x_basis) reference.Apply(0, q, 0);
            }
        }
        for (QubitIdx q{0}; q < n; ++q) {
            reference.MeasureZ(q, tableau.MeasureZ(q));
        }
    }
}

TEST_CASE("New qubits start in zero", "[inverse_tableau][add]") {
    const QubitIdx n{GENERATE(QubitIdx{0}, QubitIdx{1}, QubitIdx{4})};
    InverseTableau t{n};
    for (QubitIdx q{0}; q < n; ++q) CHECK_FALSE(t.MeasureZ(q));
    t.Add();
    CHECK_FALSE(t.MeasureZ(n));
    t.ApplyX(n);
    CHECK(t.MeasureZ(n));
}

TEST_CASE("X and Z measurements sample both outcomes and remain repeatable",
          "[inverse_tableau][measure]") {
    const bool x_basis{GENERATE(false, true)};
    std::set<bool> seen{};
    for (std::uint32_t trial{0}; trial < 64; ++trial) {
        InverseTableau t{1};
        if (!x_basis) t.ApplyH(0);
        const bool result{x_basis ? t.MeasureX(0) : t.MeasureZ(0)};
        seen.insert(result);
        CHECK((x_basis ? t.MeasureX(0) : t.MeasureZ(0)) == result);
    }
    CHECK(seen.size() == 2);
}

TEST_CASE("Bell state measurements remain correlated", "[inverse_tableau][entangle]") {
    const bool x_basis{GENERATE(false, true)};
    const bool reverse_order{GENERATE(false, true)};
    std::set<bool> seen{};
    for (std::uint32_t trial{0}; trial < 64; ++trial) {
        InverseTableau t{2};
        t.ApplyH(0);
        t.ApplyCx(0, 1);
        const QubitIdx first{reverse_order ? 1U : 0U};
        const QubitIdx second{1 - first};
        const bool result{x_basis ? t.MeasureX(first) : t.MeasureZ(first)};
        seen.insert(result);
        CHECK((x_basis ? t.MeasureX(second) : t.MeasureZ(second)) == result);
    }
    CHECK(seen.size() == 2);
}

TEST_CASE("Print", "[inverse_tableau][Print]") {
    SECTION("fresh 1-qubit tableau") {
        const InverseTableau t{1};
        CoutCapture cap{};
        t.Print();
        CHECK(cap.ss.str() ==
              "\n    | x | z | r |"
              "\nX0: | 1 | 0 | 0 |"
              "\nZ0: | 0 | 1 | 0 |"
              "\n");
    }
    SECTION("fresh 2-qubit tableau") {
        const InverseTableau t{2};
        CoutCapture cap{};
        t.Print();
        CHECK(cap.ss.str() ==
              "\n    | x x | z z | r |"
              "\nX0: | 1 0 | 0 0 | 0 |"
              "\nX1: | 0 1 | 0 0 | 0 |"
              "\nZ0: | 0 0 | 1 0 | 0 |"
              "\nZ1: | 0 0 | 0 1 | 0 |"
              "\n");
    }
    SECTION("sign column after X") {
        InverseTableau t{1};
        t.ApplyX(0);
        CoutCapture cap{};
        t.Print();
        CHECK(cap.ss.str() ==
              "\n    | x | z | r |"
              "\nX0: | 1 | 0 | 0 |"
              "\nZ0: | 0 | 1 | 1 |"
              "\n");
    }
    SECTION("rows swapped after H") {
        InverseTableau t{1};
        t.ApplyH(0);
        CoutCapture cap{};
        t.Print();
        CHECK(cap.ss.str() ==
              "\n    | x | z | r |"
              "\nX0: | 0 | 1 | 0 |"
              "\nZ0: | 1 | 0 | 0 |"
              "\n");
    }
}

TEST_CASE("Packed signs survive swaps and growth at block boundaries", "[inverse_tableau][signs]") {
    const QubitIdx n{
        GENERATE(QubitIdx{1}, QubitIdx{63}, QubitIdx{64}, QubitIdx{65}, QubitIdx{129})};
    const auto q{n - 1};
    InverseTableau t{n};
    t.ApplyX(q);
    t.ApplyH(q);
    CHECK(t.MeasureX(q));
    t.Add();
    CHECK(t.MeasureX(q));
    CHECK_FALSE(t.MeasureZ(n));
}

TEST_CASE("BasisTable appends identity rows and positive signs", "[basis][add]") {
    BasisTable table{2};
    REQUIRE(table.gens.size() == 2);
    REQUIRE(table.signs.Size() == 2);
    CHECK_FALSE(table.signs.Any());
    table.gens[0].SetX(0);
    table.signs.Set(0, true);
    table.Add();
    REQUIRE(table.gens.size() == 3);
    REQUIRE(table.signs.Size() == 3);
    CHECK(table.signs.Get(0));
    CHECK_FALSE(table.signs.Get(1));
    CHECK_FALSE(table.signs.Get(2));
    CHECK(Str(table.gens[0], table.signs.Get(0)) == "-XII");
    CHECK(Str(table.gens[1], table.signs.Get(1)) == "+III");
    CHECK(Str(table.gens[2], table.signs.Get(2)) == "+III");
}

TEST_CASE("PauliString swaps X/Z components", "[pauli][swap]") {
    const auto idx{GENERATE(0, 1, 2, 3)};
    auto g{Pauli("IXZY")};
    auto expected{std::string{"+IXZY"}};
    expected[idx + 1] = std::string{"IZXY"}[idx];
    g.SwapXZ(idx);
    CHECK(Str(g) == expected);
    g.SwapXZ(idx);
    CHECK(Str(g) == "+IXZY");
}

TEST_CASE("PauliString multiplication across block boundaries", "[pauli][xor]") {
    const auto n{GENERATE(0, 1, 63, 64, 65, 129)};
    auto a{Pauli(std::string(n, 'X'))};
    const auto b{Pauli(std::string(n, 'Y'))};
    CHECK(a.XorWithPhase(b) == n);
    CHECK(Str(a) == "+" + std::string(n, 'Z'));
    CHECK(Str(b) == "+" + std::string(n, 'Y'));
    CHECK(a.XorWithPhase(b) == 3 * n);
    CHECK(Str(a) == "+" + std::string(n, 'X'));
    CHECK(a.XorWithPhase(a) == 0);
    CHECK(Str(a) == "+" + std::string(n, 'I'));
}

TEST_CASE("Fused XOR preserves phase, result, and source", "[pauli][phase][xor]") {
    constexpr std::array<std::array<std::int64_t, 4>, 4> kPhases{{
        {0, 0, 0, 0},
        {0, 0, 1, -1},
        {0, 3, 4, 1},
        {0, 1, 3, 0},
    }};
    constexpr std::array<const char*, 4> kProducts{"IXYZ", "XIZY", "YZIX", "ZYXI"};
    const auto n{GENERATE(0, 1, 63, 64, 65, 129)};
    std::mt19937 rng{42};
    for (std::uint32_t trial{0}; trial < 32; ++trial) {
        std::string lhs(n, 'I');
        std::string rhs(n, 'I');
        std::string expected(n, 'I');
        std::int64_t exponent{0};
        std::int64_t self_exponent{0};
        for (std::size_t i{0}; i < lhs.size(); ++i) {
            const auto a{rng() % 4};
            const auto b{rng() % 4};
            lhs[i] = "IXYZ"[a];
            rhs[i] = "IXYZ"[b];
            expected[i] = kProducts[a][b];
            exponent += kPhases[a][b];
            if (expected[i] == 'Y') self_exponent += 4;
        }
        auto a{Pauli(lhs)};
        const auto b{Pauli(rhs)};
        CHECK(a.XorWithPhase(b) == exponent);
        CHECK(Str(a) == "+" + expected);
        CHECK(Str(b) == "+" + rhs);

        CHECK(a.XorWithPhase(a) == self_exponent);
        CHECK(Str(a) == "+" + std::string(n, 'I'));
    }
}

TEST_CASE("e: single-qubit table", "[pauli][phase]") {
    struct Case {
        const char* a;
        const char* b;
        std::int64_t e;
    };
    const auto c{GENERATE(values<Case>({
        {"I", "I", 0},
        {"I", "X", 0},
        {"I", "Y", 0},
        {"I", "Z", 0},
        {"X", "I", 0},
        {"X", "X", 0},
        {"X", "Y", 1},
        {"X", "Z", -1},
        {"Y", "I", 0},
        {"Y", "X", 3},
        {"Y", "Y", 4},
        {"Y", "Z", 1},
        {"Z", "I", 0},
        {"Z", "X", 1},
        {"Z", "Y", 3},
        {"Z", "Z", 0},
    }))};
    CAPTURE(c.a, c.b);
    CHECK(Pauli(c.a).XorWithPhase(Pauli(c.b)) == c.e);
}

TEST_CASE("e: matches Pauli matrix multiplication", "[pauli][phase]") {
    const std::complex<double> imaginary_unit{0, 1};
    for (const char a : std::string{"IXYZ"}) {
        for (const char b : std::string{"IXYZ"}) {
            CAPTURE(a, b);
            const PauliString ga{Pauli(std::string(1, a))};
            const PauliString gb{Pauli(std::string(1, b))};
            PauliString prod{ga};
            const auto exponent{prod.XorWithPhase(gb)};
            const char p{Str(prod)[1]};

            const std::int64_t e{((exponent % 4) + 4) % 4};
            const std::complex<double> factor{std::pow(imaginary_unit, e)};
            const Mat lhs{Multiply(PauliMatrix(a), PauliMatrix(b))};
            const Mat rhs{PauliMatrix(p)};
            for (std::size_t k{0}; k < 4; ++k) {
                CHECK(std::abs(lhs[k] - factor * rhs[k]) < 1e-12);
            }
        }
    }
}

TEST_CASE("e: is additive over tensor products", "[pauli][phase]") {
    const auto pair{GENERATE(table<std::string, std::string>({
        {"XY", "ZZ"},
        {"YYYY", "XZXZ"},
        {"IXYZ", "ZYXI"},
        {"ZZZ", "XXX"},
        {"YIY", "IYI"},
        {std::string(65, 'X'), std::string(65, 'Z')},
        {std::string(129, 'Y'), std::string(129, 'Y')},
    }))};
    const auto& [a, b]{pair};
    std::int64_t sum{0};
    for (std::size_t i{0}; i < a.size(); ++i) {
        sum += Pauli(std::string(1, a[i])).XorWithPhase(Pauli(std::string(1, b[i])));
    }
    CAPTURE(a, b);
    CHECK(Pauli(a).XorWithPhase(Pauli(b)) == sum);
}

TEST_CASE("InverseTableau rejects invalid gate and measurement indices",
          "[inverse_tableau][bounds]") {
    const QubitIdx n{GENERATE(QubitIdx{0}, QubitIdx{1}, QubitIdx{65})};
    InverseTableau tableau{n};
    for (const auto q : {n, n + 1, std::numeric_limits<QubitIdx>::max()}) {
        CHECK_THROWS_AS(tableau.ApplyH(q), std::out_of_range);
        CHECK_THROWS_AS(tableau.ApplyX(q), std::out_of_range);
        CHECK_THROWS_AS(tableau.ApplyS(q), std::out_of_range);
        CHECK_THROWS_AS(tableau.ApplyZ(q), std::out_of_range);
        CHECK_THROWS_AS(tableau.ApplyCx(q, 0), std::out_of_range);
        CHECK_THROWS_AS(tableau.ApplyCx(0, q), std::out_of_range);
        CHECK_THROWS_AS(tableau.MeasureX(q), std::out_of_range);
        CHECK_THROWS_AS(tableau.MeasureZ(q), std::out_of_range);
    }
    for (QubitIdx q{0}; q < n; ++q) CHECK_FALSE(tableau.MeasureZ(q));
}

}  // namespace stabaal
