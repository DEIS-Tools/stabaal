#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators.hpp>
#include <catch2/generators/catch_generators_range.hpp>
#include <cstddef>
#include <cstdint>
#include <iostream>
#include <random>
#include <sstream>
#include <stdexcept>
#include <string>

#include "gate_effect.hpp"
#include "inverse_tableau.hpp"
#include "types.hpp"

namespace stabaal {

// Renders a GenEffect as a signed Pauli string, e.g. "-XIZY".
template <std::size_t N>
static std::string str(const GenEffect<N> &g) {
    std::string s{g.r ? "-" : "+"};
    for (std::size_t i{0}; i < N; ++i) {
        s += g.x[i] ? (g.z[i] ? 'Y' : 'X') : (g.z[i] ? 'Z' : 'I');
    }
    return s;
}

// Every row and sign of the tableau, as printed.
static std::string rows(const InverseTableau &t) {
    std::ostringstream ss{};
    auto *const old{std::cout.rdbuf(ss.rdbuf())};
    t.print();
    std::cout.rdbuf(old);
    return ss.str();
}

// Puts a and b in the same pseudo-random state using the hand-written gates.
static void scramble(InverseTableau &a, InverseTableau &b, const QubitIdx n,
                     const std::uint32_t seed, const std::uint32_t steps = 20) {
    std::mt19937 rng{seed};
    for (std::uint32_t i{0}; i < steps; ++i) {
        const QubitIdx q{static_cast<QubitIdx>(rng() % n)};
        const QubitIdx t{static_cast<QubitIdx>((q + 1 + rng() % (n - 1)) % n)};
        switch (rng() % 5) {
            case 0:
                a.applyH(q);
                b.applyH(q);
                break;
            case 1:
                a.applyS(q);
                b.applyS(q);
                break;
            case 2:
                a.applyX(q);
                b.applyX(q);
                break;
            case 3:
                a.applyZ(q);
                b.applyZ(q);
                break;
            case 4:
                a.applyCx(q, t);
                b.applyCx(q, t);
                break;
        }
    }
}

TEST_CASE("FindGateEffect runs at compile time", "[gate_effect]") {
    constexpr auto k_h{findGateEffect<1>(kHMat)};
    STATIC_REQUIRE(k_h.computational[0].x[0]);
    STATIC_REQUIRE(!k_h.computational[0].z[0]);
    STATIC_REQUIRE(k_h.hadamard[0].z[0]);
    STATIC_REQUIRE(!k_h.hadamard[0].x[0]);
}

// computational[j] is the image of Z_j and hadamard[j] the image of X_j, both as U^dagger P U.
TEST_CASE("FindGateEffect single-qubit images", "[gate_effect]") {
    SECTION("I") {
        constexpr auto k_g{findGateEffect<1>(kIMat)};
        CHECK(str(k_g.computational[0]) == "+Z");
        CHECK(str(k_g.hadamard[0]) == "+X");
    }
    SECTION("H") {
        constexpr auto k_g{findGateEffect<1>(kHMat)};
        CHECK(str(k_g.computational[0]) == "+X");
        CHECK(str(k_g.hadamard[0]) == "+Z");
    }
    SECTION("S: X -> -Y") {
        constexpr auto k_g{findGateEffect<1>(kSMat)};
        CHECK(str(k_g.computational[0]) == "+Z");
        CHECK(str(k_g.hadamard[0]) == "-Y");
    }
    SECTION("S^dagger: X -> +Y") {
        constexpr auto k_g{findGateEffect<1>(kSMat.adjoint())};
        CHECK(str(k_g.computational[0]) == "+Z");
        CHECK(str(k_g.hadamard[0]) == "+Y");
    }
    SECTION("X") {
        constexpr auto k_g{findGateEffect<1>(kXMat)};
        CHECK(str(k_g.computational[0]) == "-Z");
        CHECK(str(k_g.hadamard[0]) == "+X");
    }
    SECTION("Z") {
        constexpr auto k_g{findGateEffect<1>(kZMat)};
        CHECK(str(k_g.computational[0]) == "+Z");
        CHECK(str(k_g.hadamard[0]) == "-X");
    }
    SECTION("Y") {
        constexpr auto k_g{findGateEffect<1>(kYMat)};
        CHECK(str(k_g.computational[0]) == "-Z");
        CHECK(str(k_g.hadamard[0]) == "-X");
    }
}

TEST_CASE("FindGateEffect CX images", "[gate_effect]") {
    constexpr auto k_g{findGateEffect<2>(kCxMat)};
    CHECK(str(k_g.computational[0]) == "+ZI");
    CHECK(str(k_g.computational[1]) == "+ZZ");
    CHECK(str(k_g.hadamard[0]) == "+XX");
    CHECK(str(k_g.hadamard[1]) == "+IX");
}

TEST_CASE("Generated methods on a fresh tableau", "[gate_effect][inverse_tableau]") {
    SECTION("H(0) gives |+>") {
        InverseTableau t{1};
        t.h(0);
        CHECK(t.measureX(0) == false);
    }
    SECTION("H(0) Z(0) gives |->") {
        InverseTableau t{1};
        t.h(0);
        t.z(0);
        CHECK(t.measureX(0) == true);
    }
    SECTION("exact rows after H(0) Cx(0, 1)") {
        InverseTableau t{2};
        t.h(0);
        t.cx(0, 1);
        CHECK(rows(t) ==
              "\n    | x x | z z | r |"
              "\nX0: | 0 1 | 1 0 | 0 |"
              "\nX1: | 0 1 | 0 0 | 0 |"
              "\nZ0: | 1 0 | 0 0 | 0 |"
              "\nZ1: | 1 0 | 0 1 | 0 |"
              "\n");
    }
    SECTION("out of range qubits are rejected") {
        InverseTableau t{2};
        CHECK_THROWS_AS(t.h(2), std::out_of_range);
        CHECK_THROWS_AS(t.cx(0, 2), std::out_of_range);
    }
}

TEST_CASE("Generated methods match the hand-written gates", "[gate_effect][inverse_tableau]") {
    const auto seed{GENERATE(range(0U, 32U))};
    CAPTURE(seed);
    InverseTableau a{4};
    InverseTableau b{4};
    scramble(a, b, 4, seed);
    REQUIRE(rows(a) == rows(b));

    SECTION("H") {
        for (QubitIdx q{0}; q < 4; ++q) {
            a.applyH(q);
            b.h(q);
        }
    }
    SECTION("S") {
        for (QubitIdx q{0}; q < 4; ++q) {
            a.applyS(q);
            b.s(q);
        }
    }
    SECTION("Z") {
        for (QubitIdx q{0}; q < 4; ++q) {
            a.applyZ(q);
            b.z(q);
        }
    }
    SECTION("Cx adjacent") {
        a.applyCx(0, 1);
        b.cx(0, 1);
    }
    SECTION("Cx reversed") {
        a.applyCx(3, 1);
        b.cx(3, 1);
    }
    SECTION("Cx non-adjacent") {
        a.applyCx(0, 3);
        b.cx(0, 3);
    }

    CHECK(rows(a) == rows(b));
}

TEST_CASE("apply<M> for gates without a generated method", "[gate_effect][inverse_tableau]") {
    const auto seed{GENERATE(range(0U, 32U))};
    CAPTURE(seed);
    InverseTableau a{3};
    InverseTableau b{3};
    scramble(a, b, 3, seed);

    SECTION("I leaves the tableau unchanged") { b.apply<kIMat>(1); }
    SECTION("X matches ApplyX") {
        a.applyX(2);
        b.apply<kXMat>(2);
    }
    SECTION("Y matches Z then X (equal up to global phase)") {
        a.applyZ(0);
        a.applyX(0);
        b.apply<kYMat>(0);
    }
    SECTION("S^dagger matches three S") {
        for (std::uint32_t i{0}; i < 3; ++i) a.applyS(1);
        b.apply<kSMat.adjoint()>(1);
    }

    CHECK(rows(a) == rows(b));
}

TEST_CASE("Generated methods match the hand-written gates on random circuits",
          "[gate_effect][inverse_tableau]") {
    const auto seed{GENERATE(range(0U, 64U))};
    CAPTURE(seed);
    std::mt19937 rng{seed};
    InverseTableau a{5};
    InverseTableau b{5};
    for (std::uint32_t i{0}; i < 40; ++i) {
        const QubitIdx q{static_cast<QubitIdx>(rng() % 5)};
        const QubitIdx t{static_cast<QubitIdx>((q + 1 + rng() % 4) % 5)};
        switch (rng() % 4) {
            case 0:
                a.applyH(q);
                b.h(q);
                break;
            case 1:
                a.applyS(q);
                b.s(q);
                break;
            case 2:
                a.applyZ(q);
                b.z(q);
                break;
            case 3:
                a.applyCx(q, t);
                b.cx(q, t);
                break;
        }
    }

    CHECK(rows(a) == rows(b));
}

}  // namespace stabaal
