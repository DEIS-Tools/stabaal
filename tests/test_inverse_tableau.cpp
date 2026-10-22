#include <array>
#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators.hpp>
#include <cmath>
#include <complex>
#include <cstddef>
#include <cstdint>
#include <iostream>
#include <limits>
#include <random>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include "stabaal/detail/basis_table.hpp"
#include "stabaal/detail/inverse_tableau.hpp"
#include "stabaal/detail/types.hpp"

namespace stabaal::detail {

// Builds unsigned operators from a Pauli string, e.g. pauli("XIZY") = X⊗I⊗Z⊗Y.
// Character i is the Pauli acting on qubit i.
static PauliString pauli(const std::string& s) {
    PauliString g{s.size()};
    for (std::size_t i{0}; i < s.size(); ++i) {
        switch (s[i]) {
            case 'I':
                break;
            case 'X':
                g.setX(i);
                break;
            case 'Z':
                g.setZ(i);
                break;
            case 'Y':
                g.setX(i);
                g.setZ(i);
                break;
            default:
                throw std::invalid_argument{"bad pauli char"};
        }
    }
    return g;
}

// Renders a generator as a signed Pauli string, e.g. "-XIZY".
static std::string str(const PauliString& g, const std::size_t n, const bool sign = false) {
    std::string s{sign ? "-" : "+"};
    for (std::size_t i{0}; i < n; ++i) {
        const bool x{g.getX(i)};
        const bool z{g.getZ(i)};
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
static Mat pauliMatrix(const char p) {
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
static Mat multiply(const Mat& a, const Mat& b) {
    return {a[0] * b[0] + a[1] * b[2], a[0] * b[1] + a[1] * b[3], a[2] * b[0] + a[3] * b[2],
            a[2] * b[1] + a[3] * b[3]};
}

struct StateVector {
    std::vector<std::complex<double>> amplitudes;

    explicit StateVector(const QubitIdx n) : amplitudes(std::size_t{1} << n) { amplitudes[0] = 1; }

    void add() { amplitudes.resize(amplitudes.size() * 2); }

    void apply(const std::uint8_t gate, const QubitIdx q, const QubitIdx target) {
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

    void measureZ(const QubitIdx q, const bool result) {
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
                tableau.add();
                reference.add();
                ++n;
            }
            const auto gate{static_cast<std::uint8_t>(choices() % (n == 1 ? 4 : 5))};
            const QubitIdx q{static_cast<QubitIdx>(choices() % n)};
            const QubitIdx target{
                static_cast<QubitIdx>((q + 1 + choices() % (n == 1 ? 1 : n - 1)) % n)};
            CAPTURE(trial, step, gate, q, target);
            switch (gate) {
                case 0:
                    tableau.applyH(q);
                    break;
                case 1:
                    tableau.applyX(q);
                    break;
                case 2:
                    tableau.applyS(q);
                    break;
                case 3:
                    tableau.applyZ(q);
                    break;
                case 4:
                    tableau.applyCx(q, target);
                    break;
            }
            reference.apply(gate, q, target);
            if (step % 8 == 7) {
                const bool x_basis{(choices() & 1U) != 0};
                if (x_basis) reference.apply(0, q, 0);
                const bool result{x_basis ? tableau.measureX(q) : tableau.measureZ(q)};
                reference.measureZ(q, result);
                CHECK((x_basis ? tableau.measureX(q) : tableau.measureZ(q)) == result);
                if (x_basis) reference.apply(0, q, 0);
            }
        }
        for (QubitIdx q{0}; q < n; ++q) {
            reference.measureZ(q, tableau.measureZ(q));
        }
    }
}

TEST_CASE("New qubits start in zero", "[inverse_tableau][add]") {
    const QubitIdx n{GENERATE(QubitIdx{0}, QubitIdx{1}, QubitIdx{4})};
    InverseTableau t{n};
    for (QubitIdx q{0}; q < n; ++q) CHECK_FALSE(t.measureZ(q));
    t.add();
    CHECK_FALSE(t.measureZ(n));
    t.applyX(n);
    CHECK(t.measureZ(n));
}

TEST_CASE("X and Z measurements sample both outcomes and remain repeatable",
          "[inverse_tableau][measure]") {
    const bool x_basis{GENERATE(false, true)};
    std::set<bool> seen{};
    for (std::uint32_t trial{0}; trial < 64; ++trial) {
        InverseTableau t{1};
        if (!x_basis) t.applyH(0);
        const bool result{x_basis ? t.measureX(0) : t.measureZ(0)};
        seen.insert(result);
        CHECK((x_basis ? t.measureX(0) : t.measureZ(0)) == result);
    }
    CHECK(seen.size() == 2);
}

TEST_CASE("Bell state measurements remain correlated", "[inverse_tableau][entangle]") {
    const bool x_basis{GENERATE(false, true)};
    const bool reverse_order{GENERATE(false, true)};
    std::set<bool> seen{};
    for (std::uint32_t trial{0}; trial < 64; ++trial) {
        InverseTableau t{2};
        t.applyH(0);
        t.applyCx(0, 1);
        const QubitIdx first{reverse_order ? 1U : 0U};
        const QubitIdx second{1 - first};
        const bool result{x_basis ? t.measureX(first) : t.measureZ(first)};
        seen.insert(result);
        CHECK((x_basis ? t.measureX(second) : t.measureZ(second)) == result);
    }
    CHECK(seen.size() == 2);
}

TEST_CASE("Print", "[inverse_tableau][Print]") {
    SECTION("fresh 1-qubit tableau") {
        const InverseTableau t{1};
        CoutCapture cap{};
        t.print();
        CHECK(cap.ss.str() ==
              "\n    | x | z | r |"
              "\nX0: | 1 | 0 | 0 |"
              "\nZ0: | 0 | 1 | 0 |"
              "\n");
    }
    SECTION("fresh 2-qubit tableau") {
        const InverseTableau t{2};
        CoutCapture cap{};
        t.print();
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
        t.applyX(0);
        CoutCapture cap{};
        t.print();
        CHECK(cap.ss.str() ==
              "\n    | x | z | r |"
              "\nX0: | 1 | 0 | 0 |"
              "\nZ0: | 0 | 1 | 1 |"
              "\n");
    }
    SECTION("rows swapped after H") {
        InverseTableau t{1};
        t.applyH(0);
        CoutCapture cap{};
        t.print();
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
    t.applyX(q);
    t.applyH(q);
    CHECK(t.measureX(q));
    t.add();
    CHECK(t.measureX(q));
    CHECK_FALSE(t.measureZ(n));
}

TEST_CASE("BasisTable appends identity rows and positive signs", "[basis][add]") {
    BasisTable table{2};
    CHECK_FALSE(table.signs.any(2));
    table.gens[0].setX(0);
    table.signs.set(0, true);
    table.add(2);
    CHECK(table.signs.get(0));
    CHECK_FALSE(table.signs.get(1));
    CHECK_FALSE(table.signs.get(2));
    CHECK(str(table.gens[0], 3, table.signs.get(0)) == "-XII");
    CHECK(str(table.gens[1], 3, table.signs.get(1)) == "+III");
    CHECK(str(table.gens[2], 3, table.signs.get(2)) == "+III");
}

TEST_CASE("BasisTable growth reuses capacity and preserves rows", "[basis][growth]") {
    BasisTable table{0};
    CHECK_FALSE(table.gens);
    CHECK(BasisTable::allocatedBytes(0) == 0);
    for (std::size_t n{0}; n < 130; ++n) {
        const auto* previous{table.gens.data()};
        table.add(n);
        if (n != 0 && (n & (n - 1)) != 0) CHECK(table.gens.data() == previous);
        for (std::size_t row{0}; row <= n; ++row) {
            CHECK(table.signs.get(row) == (row < n));
            for (std::size_t col{0}; col <= n; ++col) {
                CHECK(table.gens[row].getX(col) == (row < n && row == col));
                CHECK_FALSE(table.gens[row].getZ(col));
            }
        }
        table.gens[n].setX(n);
        table.signs.set(n);
    }
}

TEST_CASE("InverseTableau copies remain independent after growth", "[inverse_tableau][copy]") {
    const auto n{GENERATE(QubitIdx{0}, QubitIdx{64}, QubitIdx{65})};
    InverseTableau original{n};
    original.add();
    original.applyX(n);
    InverseTableau copy{original};
    InverseTableau assigned{0};
    assigned = original;
    original.applyX(n);
    CHECK_FALSE(original.measureZ(n));
    CHECK(copy.measureZ(n));
    CHECK(assigned.measureZ(n));
    copy.add();
    assigned.add();
    CHECK_FALSE(copy.measureZ(n + 1));
    CHECK_FALSE(assigned.measureZ(n + 1));
    CHECK(copy.measureZ(n));
    CHECK(assigned.measureZ(n));
}

TEST_CASE("BasisTable copies preserve independent signs and rows", "[basis][copy]") {
    const auto n{GENERATE(QubitIdx{0}, QubitIdx{64}, QubitIdx{65})};
    BasisTable original{n};
    original.add(n);
    original.signs.set(n);
    original.gens[n].setX(n);

    BasisTable copy{original, n + 1};
    BasisTable assigned{0};
    assigned = BasisTable{original, n + 1};
    original.signs.set(n, false);
    original.gens[n].setX(n, false);
    CHECK(copy.signs.get(n));
    CHECK(copy.gens[n].getX(n));
    CHECK(assigned.signs.get(n));
    CHECK(assigned.gens[n].getX(n));

    BasisTable moved{std::move(copy)};
    CHECK(moved.signs.get(n));
    copy.add(0);
    CHECK_FALSE(copy.signs.get(0));
    assigned = std::move(moved);
    CHECK(assigned.gens[n].getX(n));
    moved.add(0);
    CHECK_FALSE(moved.signs.get(0));
}

TEST_CASE("PauliString copies and moves preserve storage ownership", "[pauli][copy]") {
    const auto n{GENERATE(QubitIdx{0}, QubitIdx{64}, QubitIdx{65}, QubitIdx{129})};
    auto original{pauli(std::string(n, 'Z'))};
    PauliString copy{original, n};
    PauliString assigned;
    assigned = PauliString{original, n};
    original.addI(n);
    CHECK(str(copy, n) == "+" + std::string(n, 'Z'));
    CHECK(str(assigned, n) == str(copy, n));

    PauliString moved{std::move(copy)};
    CHECK(str(moved, n) == str(assigned, n));
    copy.addI(0);
    CHECK(str(copy, 1) == "+I");
    assigned = std::move(moved);
    CHECK(str(assigned, n) == "+" + std::string(n, 'Z'));
    moved.addI(0);
    CHECK(str(moved, 1) == "+I");
    if (n != 0) {
        assigned.setX(n - 1);
        CHECK_FALSE(original.getX(n - 1));
    }
}

TEST_CASE("PauliString grows across word boundaries", "[pauli][growth]") {
    const auto n{GENERATE(0, 1, 31, 32, 33, 63, 64, 65, 127, 128, 129)};
    auto g{pauli(std::string(n, 'Z'))};
    CHECK_FALSE(g.anyX(n));
    g.addI(n);
    CHECK(str(g, n + 1) == "+" + std::string(n, 'Z') + "I");
    g.setX(n);
    CHECK(g.anyX(n + 1));
    g.swapXZ(n);
    CHECK_FALSE(g.anyX(n + 1));
    CHECK(g.getZ(n));
    g.setZ(n, false);
    CHECK(str(g, n + 1) == "+" + std::string(n, 'Z') + "I");
    g.setX(n);
    g.setX(n, false);
    CHECK_FALSE(g.anyX(n + 1));
}

TEST_CASE("PauliString swaps X/Z components", "[pauli][swap]") {
    const auto idx{GENERATE(0, 1, 2, 3)};
    auto g{pauli("IXZY")};
    auto expected{std::string{"+IXZY"}};
    expected[idx + 1] = std::string{"IZXY"}[idx];
    g.swapXZ(idx);
    CHECK(str(g, 4) == expected);
    g.swapXZ(idx);
    CHECK(str(g, 4) == "+IXZY");
}

TEST_CASE("PauliString multiplication across block boundaries", "[pauli][xor]") {
    const auto n{GENERATE(0, 1, 63, 64, 65, 129)};
    auto a{pauli(std::string(n, 'X'))};
    const auto b{pauli(std::string(n, 'Y'))};
    CHECK(a.xorWithPhase(b, n) == n);
    CHECK(str(a, n) == "+" + std::string(n, 'Z'));
    CHECK(str(b, n) == "+" + std::string(n, 'Y'));
    CHECK(a.xorWithPhase(b, n) == 3 * n);
    CHECK(str(a, n) == "+" + std::string(n, 'X'));
    CHECK(a.xorWithPhase(a, n) == 0);
    CHECK(str(a, n) == "+" + std::string(n, 'I'));
}

TEST_CASE("Fused XOR preserves phase, result, and source", "[pauli][phase][xor]") {
    constexpr std::array<std::array<std::int64_t, 4>, 4> k_phases{{
        {0, 0, 0, 0},
        {0, 0, 1, -1},
        {0, 3, 4, 1},
        {0, 1, 3, 0},
    }};
    constexpr std::array<const char*, 4> k_products{"IXYZ", "XIZY", "YZIX", "ZYXI"};
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
            expected[i] = k_products[a][b];
            exponent += k_phases[a][b];
            if (expected[i] == 'Y') self_exponent += 4;
        }
        auto a{pauli(lhs)};
        const auto b{pauli(rhs)};
        CHECK(a.xorWithPhase(b, n) == exponent);
        CHECK(str(a, n) == "+" + expected);
        CHECK(str(b, n) == "+" + rhs);

        CHECK(a.xorWithPhase(a, n) == self_exponent);
        CHECK(str(a, n) == "+" + std::string(n, 'I'));
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
    CHECK(pauli(c.a).xorWithPhase(pauli(c.b), 1) == c.e);
}

TEST_CASE("e: matches Pauli matrix multiplication", "[pauli][phase]") {
    const std::complex<double> imaginary_unit{0, 1};
    for (const char a : std::string{"IXYZ"}) {
        for (const char b : std::string{"IXYZ"}) {
            CAPTURE(a, b);
            const PauliString ga{pauli(std::string(1, a))};
            const PauliString gb{pauli(std::string(1, b))};
            PauliString prod{ga, 1};
            const auto exponent{prod.xorWithPhase(gb, 1)};
            const char p{str(prod, 1)[1]};

            const std::int64_t e{((exponent % 4) + 4) % 4};
            const std::complex<double> factor{std::pow(imaginary_unit, e)};
            const Mat lhs{multiply(pauliMatrix(a), pauliMatrix(b))};
            const Mat rhs{pauliMatrix(p)};
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
        sum += pauli(std::string(1, a[i])).xorWithPhase(pauli(std::string(1, b[i])), 1);
    }
    CAPTURE(a, b);
    CHECK(pauli(a).xorWithPhase(pauli(b), a.size()) == sum);
}

TEST_CASE("InverseTableau rejects invalid gate and measurement indices",
          "[inverse_tableau][bounds]") {
    const QubitIdx n{GENERATE(QubitIdx{0}, QubitIdx{1}, QubitIdx{65})};
    InverseTableau tableau{n};
    for (const auto q : {n, n + 1, std::numeric_limits<QubitIdx>::max()}) {
        CHECK_THROWS_AS(tableau.applyH(q), std::out_of_range);
        CHECK_THROWS_AS(tableau.applyX(q), std::out_of_range);
        CHECK_THROWS_AS(tableau.applyS(q), std::out_of_range);
        CHECK_THROWS_AS(tableau.applyZ(q), std::out_of_range);
        CHECK_THROWS_AS(tableau.applyCx(q, 0), std::out_of_range);
        CHECK_THROWS_AS(tableau.applyCx(0, q), std::out_of_range);
        CHECK_THROWS_AS(tableau.measureX(q), std::out_of_range);
        CHECK_THROWS_AS(tableau.measureZ(q), std::out_of_range);
    }
    for (QubitIdx q{0}; q < n; ++q) CHECK_FALSE(tableau.measureZ(q));
}

}  // namespace stabaal::detail
