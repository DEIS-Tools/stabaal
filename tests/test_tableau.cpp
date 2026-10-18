#include <array>
#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators.hpp>
#include <cmath>
#include <complex>
#include <cstddef>
#include <cstdint>
#include <random>
#include <set>
#include <stdexcept>
#include <string>
#include <vector>

#include "basis_table.hpp"
#include "tableau.hpp"
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
static std::string Str(const PauliString& g, bool sign = false) {
    std::string s{sign ? "-" : "+"};
    auto n{g.Size()};
    for (std::size_t i{0}; i < n; ++i) {
        bool x{g.GetX(i)};
        bool z{g.GetZ(i)};
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
static Mat PauliMatrix(char p) {
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

    explicit StateVector(QubitIdx n) : amplitudes(std::size_t{1} << n) { amplitudes[0] = 1; }

    void Add() { amplitudes.resize(amplitudes.size() * 2); }

    void Apply(std::uint8_t gate, QubitIdx q, QubitIdx target) {
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

    void MeasureZ(QubitIdx q, bool result) {
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

TEST_CASE("Tableau gates and measurements agree with state-vector simulation", "[tableau]") {
    const QubitIdx initial_qubits{GENERATE(QubitIdx{1}, QubitIdx{2}, QubitIdx{3}, QubitIdx{4})};
    std::mt19937 choices{1234};
    for (std::uint32_t trial{0}; trial < 16; ++trial) {
        auto n{initial_qubits};
        Tableau tableau{n};
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

TEST_CASE("New qubits start in zero", "[tableau][add]") {
    QubitIdx n{GENERATE(QubitIdx{0}, QubitIdx{1}, QubitIdx{4})};
    Tableau t{n};
    for (QubitIdx q{0}; q < n; ++q) CHECK_FALSE(t.MeasureZ(q));
    t.Add();
    CHECK_FALSE(t.MeasureZ(n));
    t.ApplyX(n);
    CHECK(t.MeasureZ(n));
}

TEST_CASE("X and Z measurements sample both outcomes and remain repeatable", "[tableau][measure]") {
    const bool x_basis{GENERATE(false, true)};
    std::set<bool> seen{};
    for (std::uint32_t trial{0}; trial < 64; ++trial) {
        Tableau t{1};
        if (!x_basis) t.ApplyH(0);
        const bool result{x_basis ? t.MeasureX(0) : t.MeasureZ(0)};
        seen.insert(result);
        CHECK((x_basis ? t.MeasureX(0) : t.MeasureZ(0)) == result);
    }
    CHECK(seen.size() == 2);
}

TEST_CASE("Bell state measurements remain correlated", "[tableau][entangle]") {
    const bool x_basis{GENERATE(false, true)};
    const bool reverse_order{GENERATE(false, true)};
    std::set<bool> seen{};
    for (std::uint32_t trial{0}; trial < 64; ++trial) {
        Tableau t{2};
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

TEST_CASE("Print", "[tableau][Print]") {
    SECTION("fresh 1-qubit tableau") {
        const Tableau t{1};
        CoutCapture cap{};
        t.Print();
        CHECK(cap.ss.str() ==
              "\n    | x | z | r |"
              "\nX0: | 1 | 0 | 0 |"
              "\nZ0: | 0 | 1 | 0 |"
              "\n");
    }
    SECTION("fresh 2-qubit tableau") {
        Tableau t{2};
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
        Tableau t{1};
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
        Tableau t{1};
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

TEST_CASE("Packed signs survive swaps and growth at block boundaries", "[tableau][signs]") {
    QubitIdx n{GENERATE(QubitIdx{1}, QubitIdx{63}, QubitIdx{64}, QubitIdx{65}, QubitIdx{129})};
    const auto q{n - 1};
    Tableau t{n};
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
    REQUIRE(table.signs == boost::dynamic_bitset<>{2});
    table.gens[0].SetX(0);
    table.signs[0] = true;
    table.Add();
    REQUIRE(table.gens.size() == 3);
    REQUIRE(table.signs == boost::dynamic_bitset<>{3, 1});
    CHECK(Str(table.gens[0], table.signs[0]) == "-XII");
    CHECK(Str(table.gens[1], table.signs[1]) == "+III");
    CHECK(Str(table.gens[2], table.signs[2]) == "+III");
}

TEST_CASE("PauliString swaps X/Z components", "[pauli][swap]") {
    auto idx{GENERATE(0, 1, 2, 3)};
    auto g{Pauli("IXZY")};
    auto expected{std::string{"+IXZY"}};
    expected[idx + 1] = std::string{"IZXY"}[idx];
    g.SwapXZ(idx);
    CHECK(Str(g) == expected);
    g.SwapXZ(idx);
    CHECK(Str(g) == "+IXZY");
}

TEST_CASE("PauliString XOR across block boundaries", "[pauli][xor]") {
    auto n{GENERATE(0, 1, 63, 64, 65, 129)};
    auto a{Pauli(std::string(n, 'X'))};
    const auto b{Pauli(std::string(n, 'Y'))};
    a.XorWith(b);
    CHECK(Str(a) == "+" + std::string(n, 'Z'));
    CHECK(Str(b) == "+" + std::string(n, 'Y'));
    a.XorWith(b);
    CHECK(Str(a) == "+" + std::string(n, 'X'));
    a.XorWith(a);
    CHECK(Str(a) == "+" + std::string(n, 'I'));
}

TEST_CASE("e: single-qubit table", "[pauli][phase]") {
    struct Case {
        const char* a;
        const char* b;
        std::int64_t e;
    };
    auto c{GENERATE(values<Case>({
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
    CHECK(Pauli(c.a).PhaseExponent(Pauli(c.b)) == c.e);
}

TEST_CASE("e: matches Pauli matrix multiplication", "[pauli][phase]") {
    const std::complex<double> imaginary_unit{0, 1};
    for (char a : std::string{"IXYZ"}) {
        for (char b : std::string{"IXYZ"}) {
            CAPTURE(a, b);
            PauliString ga{Pauli(std::string(1, a))};
            PauliString gb{Pauli(std::string(1, b))};
            PauliString prod{ga};
            prod.XorWith(gb);
            char p{Str(prod)[1]};

            std::int64_t e{((ga.PhaseExponent(gb) % 4) + 4) % 4};
            std::complex<double> factor{std::pow(imaginary_unit, e)};
            Mat lhs{Multiply(PauliMatrix(a), PauliMatrix(b))};
            Mat rhs{PauliMatrix(p)};
            for (std::size_t k{0}; k < 4; ++k) {
                CHECK(std::abs(lhs[k] - factor * rhs[k]) < 1e-12);
            }
        }
    }
}

TEST_CASE("e: is additive over tensor products", "[pauli][phase]") {
    auto pair{GENERATE(table<std::string, std::string>({
        {"XY", "ZZ"},
        {"YYYY", "XZXZ"},
        {"IXYZ", "ZYXI"},
        {"ZZZ", "XXX"},
        {"YIY", "IYI"},
        {std::string(65, 'X'), std::string(65, 'Z')},
        {std::string(129, 'Y'), std::string(129, 'Y')},
    }))};
    auto [a, b]{pair};
    std::int64_t sum{0};
    for (std::size_t i{0}; i < a.size(); ++i) {
        sum += Pauli(std::string(1, a[i])).PhaseExponent(Pauli(std::string(1, b[i])));
    }
    CAPTURE(a, b);
    CHECK(Pauli(a).PhaseExponent(Pauli(b)) == sum);
}

}  // namespace stabaal
