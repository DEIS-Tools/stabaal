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

#include "stabaal/detail/inverse_tableau.hpp"
#include "stabaal/simulator.hpp"
#include "stabaal/types.hpp"

namespace stabaal::detail {

// Redirects std::cout for the lifetime of the object.
struct CoutCapture {
    std::ostringstream ss{};
    std::streambuf* old;
    CoutCapture() : old{std::cout.rdbuf(ss.rdbuf())} {}
    ~CoutCapture() { std::cout.rdbuf(old); }
};

struct StateVector {
    std::vector<std::complex<double>> amplitudes;

    explicit StateVector(const QubitSize n) : amplitudes(std::size_t{1} << n) { amplitudes[0] = 1; }

    void Add() { amplitudes.resize(amplitudes.size() * 2); }

    void Apply(const std::uint8_t gate, const QubitSize q, const QubitSize target) {
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
                case 4: {
                    const auto a{amplitudes[i]};
                    amplitudes[i] = std::complex<double>{0, -1} * amplitudes[j];
                    amplitudes[j] = std::complex<double>{0, 1} * a;
                    break;
                }
                case 5:
                    amplitudes[j] *= std::complex<double>{0, -1};
                    break;
                case 6:
                    break;
                case 7:
                    if ((j & (std::size_t{1} << target)) == 0) {
                        std::swap(amplitudes[j], amplitudes[j | (std::size_t{1} << target)]);
                    }
                    break;
                case 8:
                    if ((j & (std::size_t{1} << target)) != 0) amplitudes[j] = -amplitudes[j];
                    break;
                case 9:
                    if ((j & (std::size_t{1} << target)) == 0) {
                        const auto k{j | (std::size_t{1} << target)};
                        const auto a{amplitudes[j]};
                        amplitudes[j] = std::complex<double>{0, -1} * amplitudes[k];
                        amplitudes[k] = std::complex<double>{0, 1} * a;
                    }
                    break;
            }
        }
    }

    void MeasureZ(const QubitSize q, const bool result) {
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
    const QubitSize initial_qubits{
        GENERATE(QubitSize{1}, QubitSize{2}, QubitSize{3}, QubitSize{4})};
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
            const auto gate{static_cast<std::uint8_t>(choices() % (n == 1 ? 7 : 10))};
            const QubitSize q{static_cast<QubitSize>(choices() % n)};
            const QubitSize target{
                static_cast<QubitSize>((q + 1 + choices() % (n == 1 ? 1 : n - 1)) % n)};
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
                    tableau.ApplyY(q);
                    break;
                case 5:
                    tableau.ApplySdg(q);
                    break;
                case 6:
                    tableau.ApplyI(q);
                    break;
                case 7:
                    tableau.ApplyCx(q, target);
                    break;
                case 8:
                    tableau.ApplyCz(q, target);
                    break;
                case 9:
                    tableau.ApplyCy(q, target);
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
        for (QubitSize q{0}; q < n; ++q) {
            reference.MeasureZ(q, tableau.MeasureZ(q));
        }
    }
}

TEST_CASE("Y flips X and Z eigenvalues and preserves Y eigenvalues", "[simulator][Y]") {
    const auto n{GENERATE(QubitSize{1}, QubitSize{64}, QubitSize{65}, QubitSize{129})};
    const auto basis{GENERATE(0, 1, 2)};  // Z, X, Y
    const bool negative{GENERATE(false, true)};
    const auto q{n - 1};
    Simulator<> sim{n};
    if (negative) sim.X(q);
    if (basis != 0) sim.H(q);
    if (basis == 2) sim.S(q);
    sim.Y(q);
    if (basis == 2) {
        sim.S(q);
        sim.Z(q);
    }
    CHECK((basis == 0 ? sim.MeasureZ(q) : sim.MeasureX(q)) == (negative != (basis != 2)));
    if (q != 0) CHECK_FALSE(sim.MeasureZ(0));
}

TEST_CASE("Sdg maps X to negative Y, Y to X, and preserves Z", "[simulator][Sdg]") {
    const auto n{
        GENERATE(QubitSize{1}, QubitSize{63}, QubitSize{64}, QubitSize{65}, QubitSize{129})};
    const auto basis{GENERATE(0, 1, 2)};  // Z, X, Y
    const bool negative{GENERATE(false, true)};
    const auto q{n - 1};
    CAPTURE(n, basis, negative);
    Simulator<> sim{n};
    if (negative) sim.X(q);
    if (basis != 0) sim.H(q);
    if (basis == 2) sim.S(q);
    sim.Sdg(q);
    if (basis == 1) {
        sim.S(q);
        sim.Z(q);
    }
    CHECK((basis == 0 ? sim.MeasureZ(q) : sim.MeasureX(q)) == (negative != (basis == 1)));
    if (q != 0) CHECK_FALSE(sim.MeasureZ(0));
}

TEST_CASE("Sdg cancels an entangled phase across word boundaries", "[simulator][Sdg][entangle]") {
    const auto n{GENERATE(QubitSize{2}, QubitSize{64}, QubitSize{65}, QubitSize{129})};
    const bool reverse{GENERATE(false, true)};
    const bool negative{GENERATE(false, true)};
    const bool x_basis{GENERATE(false, true)};
    const QubitSize c{reverse ? n - 1 : 0};
    const QubitSize t{reverse ? 0 : n - 1};
    CAPTURE(n, c, t, negative, x_basis);
    Simulator<> sim{n};
    sim.H(c);
    sim.Cx(c, t);
    if (negative) sim.Z(c);
    sim.S(c);
    sim.Sdg(t);
    const auto first{x_basis ? sim.MeasureX(c) : sim.MeasureZ(c)};
    const auto second{x_basis ? sim.MeasureX(t) : sim.MeasureZ(t)};
    CHECK((first != second) == (x_basis && negative));
}

TEST_CASE("Sx preserves X, maps Y to Z, and maps Z to negative Y", "[simulator][Sx]") {
    const auto n{
        GENERATE(QubitSize{1}, QubitSize{63}, QubitSize{64}, QubitSize{65}, QubitSize{129})};
    const auto basis{GENERATE(0, 1, 2)};  // Z, X, Y
    const bool negative{GENERATE(false, true)};
    const auto q{n - 1};
    CAPTURE(n, basis, negative);
    Simulator<> sim{n};
    if (negative) sim.X(q);
    if (basis != 0) sim.H(q);
    if (basis == 2) sim.S(q);
    sim.Sx(q);
    if (basis == 0) sim.Sdg(q);
    CHECK((basis == 2 ? sim.MeasureZ(q) : sim.MeasureX(q)) == (negative != (basis == 0)));
    for (QubitSize spectator{0}; spectator < q; ++spectator) {
        CHECK_FALSE(sim.MeasureZ(spectator));
    }
}

TEST_CASE("Two Sx gates equal X and four equal identity", "[simulator][Sx]") {
    const auto n{GENERATE(QubitSize{1}, QubitSize{64}, QubitSize{65}, QubitSize{129})};
    const auto basis{GENERATE(0, 1, 2)};  // Z, X, Y
    const bool negative{GENERATE(false, true)};
    const auto repetitions{GENERATE(2, 4)};
    const auto q{n - 1};
    CAPTURE(n, basis, negative, repetitions);
    Simulator<> sim{n};
    if (negative) sim.X(q);
    if (basis != 0) sim.H(q);
    if (basis == 2) sim.S(q);
    for (int i{0}; i < repetitions; ++i) sim.Sx(q);
    if (basis == 2) sim.Sdg(q);
    CHECK((basis == 0 ? sim.MeasureZ(q) : sim.MeasureX(q)) ==
          (negative != (repetitions == 2 && basis != 1)));
}

TEST_CASE("Sx cancels H Sdg H on an entangled qubit across word boundaries",
          "[simulator][Sx][entangle]") {
    const auto n{GENERATE(QubitSize{2}, QubitSize{64}, QubitSize{65}, QubitSize{129})};
    const bool reverse{GENERATE(false, true)};
    const bool negative{GENERATE(false, true)};
    const bool flipped{GENERATE(false, true)};
    const QubitSize c{reverse ? n - 1 : 0};
    const QubitSize t{reverse ? 0 : n - 1};
    CAPTURE(n, c, t, negative, flipped);
    Simulator<> sim{n};
    if (negative) sim.X(c);
    if (flipped) sim.X(t);
    sim.H(c);
    sim.Cx(c, t);
    sim.Sx(t);
    sim.H(t);
    sim.Sdg(t);
    sim.H(t);
    sim.Cx(c, t);
    sim.H(c);
    CHECK(sim.MeasureZ(c) == negative);
    CHECK(sim.MeasureZ(t) == flipped);
}

TEST_CASE("Identity preserves an entangled state and its phase", "[simulator][I]") {
    const auto n{GENERATE(QubitSize{2}, QubitSize{65}, QubitSize{129})};
    const bool x_basis{GENERATE(false, true)};
    Simulator<> sim{n};
    sim.H(0);
    sim.Cx(0, n - 1);
    sim.Z(0);
    for (QubitSize q{0}; q < n; ++q) sim.I(q);
    const auto first{x_basis ? sim.MeasureX(0) : sim.MeasureZ(0)};
    const auto second{x_basis ? sim.MeasureX(n - 1) : sim.MeasureZ(n - 1)};
    CHECK((first != second) == x_basis);
}

TEST_CASE("CZ and CY apply the target Pauli only when the control is one", "[simulator][Cz][Cy]") {
    const auto n{GENERATE(QubitSize{2}, QubitSize{65}, QubitSize{129})};
    const bool reverse{GENERATE(false, true)};
    const bool cy{GENERATE(false, true)};
    const bool control{GENERATE(false, true)};
    const bool negative{GENERATE(false, true)};
    const auto basis{GENERATE(0, 1, 2)};  // Z, X, Y
    const QubitSize c{reverse ? n - 1 : 0};
    const QubitSize t{reverse ? 0 : n - 1};
    CAPTURE(n, cy, c, t, control, negative, basis);
    Simulator<> sim{n};
    if (control) sim.X(c);
    if (negative) sim.X(t);
    if (basis != 0) sim.H(t);
    if (basis == 2) sim.S(t);
    if (cy) sim.Cy(c, t);
    else sim.Cz(c, t);
    if (basis == 2) {
        sim.S(t);
        sim.Z(t);
    }
    const bool flips{cy ? basis != 2 : basis != 0};
    CHECK((basis == 0 ? sim.MeasureZ(t) : sim.MeasureX(t)) == (negative != (control && flips)));
    CHECK(sim.MeasureZ(c) == control);
}

TEST_CASE("CZ and CY preserve coherent control phases", "[simulator][Cz][Cy][entangle]") {
    const bool cy{GENERATE(false, true)};
    const bool x_basis{GENERATE(false, true)};
    const bool reverse{GENERATE(false, true)};
    const QubitSize c{reverse ? 1U : 0U};
    const QubitSize t{1 - c};
    CAPTURE(cy, x_basis, c, t);
    Simulator<> sim{2};
    sim.H(c);
    if (cy) {
        sim.Cy(c, t);
        sim.S(t);
        sim.Z(t);
    } else {
        sim.H(t);
        sim.Cz(c, t);
        sim.H(t);
    }
    const auto first{x_basis ? sim.MeasureX(c) : sim.MeasureZ(c)};
    CHECK((x_basis ? sim.MeasureX(t) : sim.MeasureZ(t)) == first);
}

TEST_CASE("New qubits start in zero", "[inverse_tableau][add]") {
    const QubitSize n{GENERATE(QubitSize{0}, QubitSize{1}, QubitSize{4})};
    InverseTableau t{n};
    for (QubitSize q{0}; q < n; ++q) CHECK_FALSE(t.MeasureZ(q));
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
        const QubitSize first{reverse_order ? 1U : 0U};
        const QubitSize second{1 - first};
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
    const QubitSize n{
        GENERATE(QubitSize{1}, QubitSize{63}, QubitSize{64}, QubitSize{65}, QubitSize{129})};
    const auto q{n - 1};
    InverseTableau t{n};
    t.ApplyX(q);
    t.ApplyH(q);
    CHECK(t.MeasureX(q));
    t.Add();
    CHECK(t.MeasureX(q));
    CHECK_FALSE(t.MeasureZ(n));
}

TEST_CASE("InverseTableau copies remain independent after growth", "[inverse_tableau][copy]") {
    const auto n{GENERATE(QubitSize{0}, QubitSize{64}, QubitSize{65})};
    InverseTableau original{n};
    original.Add();
    original.ApplyX(n);
    InverseTableau copy{original};
    InverseTableau assigned{0};
    assigned = original;
    original.ApplyX(n);
    CHECK_FALSE(original.MeasureZ(n));
    CHECK(copy.MeasureZ(n));
    CHECK(assigned.MeasureZ(n));
    copy.Add();
    assigned.Add();
    CHECK_FALSE(copy.MeasureZ(n + 1));
    CHECK_FALSE(assigned.MeasureZ(n + 1));
    CHECK(copy.MeasureZ(n));
    CHECK(assigned.MeasureZ(n));
}

TEST_CASE("Growing tableaus preserve row swaps and copies discard spare capacity",
          "[inverse_tableau][copy][growth]") {
    const auto initial{
        GENERATE(QubitSize{0}, QubitSize{32}, QubitSize{63}, QubitSize{64}, QubitSize{128})};
    InverseTableau original{initial};
    for (QubitSize q{0}; q < initial + 145; ++q) {
        if (q >= initial) original.Add();
        if (q % 3 == 0) original.ApplyX(q);
        if (q % 2 != 0) original.ApplyH(q);
        if (q < initial || q % 17 != 0) continue;

        InverseTableau copy{original};
        InverseTableau assigned{1};
        assigned = original;
        const InverseTableau exact{q + 1};
        CHECK(copy.MemoryUsage() == exact.MemoryUsage());
        CHECK(assigned.MemoryUsage() == exact.MemoryUsage());
        CHECK(copy.MemoryUsage() <= original.MemoryUsage());
        if (q % 2 != 0) original.ApplyZ(q);
        else original.ApplyX(q);
        const auto outcome{q % 2 != 0 ? original.MeasureX(q) : original.MeasureZ(q)};
        CHECK(outcome != (q % 3 == 0));
        for (QubitSize i{0}; i <= q; ++i) {
            CHECK((i % 2 != 0 ? copy.MeasureX(i) : copy.MeasureZ(i)) == (i % 3 == 0));
            CHECK((i % 2 != 0 ? assigned.MeasureX(i) : assigned.MeasureZ(i)) == (i % 3 == 0));
        }
        if (q % 2 != 0) original.ApplyZ(q);
        else original.ApplyX(q);
        copy.Add();
        CHECK_FALSE(copy.MeasureZ(q + 1));
    }
}

TEST_CASE("Entanglement survives repeated storage growth across word boundaries",
          "[inverse_tableau][measurement][growth]") {
    const auto n{GENERATE(QubitSize{65}, QubitSize{129}, QubitSize{257})};
    InverseTableau tableau{0};
    tableau.Add();
    tableau.ApplyH(0);
    for (QubitSize q{1}; q < n; ++q) {
        tableau.Add();
        tableau.ApplyCx(0, q);
    }
    auto z_copy{tableau};
    const auto outcome{z_copy.MeasureZ(0)};
    for (QubitSize q{1}; q < n; ++q) CHECK(z_copy.MeasureZ(q) == outcome);
    bool parity{false};
    for (QubitSize q{0}; q < n; ++q) parity = parity != tableau.MeasureX(q);
    CHECK_FALSE(parity);
}

TEST_CASE("Simulator rejects invalid gate and measurement indices", "[simulator][bounds]") {
    const QubitSize n{GENERATE(QubitSize{0}, QubitSize{1}, QubitSize{65})};
    Simulator<> tableau{n};
    for (const auto q : {n, n + 1, std::numeric_limits<QubitSize>::max()}) {
        CHECK_THROWS_AS(tableau.H(q), std::out_of_range);
        CHECK_THROWS_AS(tableau.X(q), std::out_of_range);
        CHECK_THROWS_AS(tableau.S(q), std::out_of_range);
        CHECK_THROWS_AS(tableau.Sdg(q), std::out_of_range);
        CHECK_THROWS_AS(tableau.Sx(q), std::out_of_range);
        CHECK_THROWS_AS(tableau.Z(q), std::out_of_range);
        CHECK_THROWS_AS(tableau.Y(q), std::out_of_range);
        CHECK_THROWS_AS(tableau.I(q), std::out_of_range);
        CHECK_THROWS_AS(tableau.Cx(q, 0), std::out_of_range);
        CHECK_THROWS_AS(tableau.Cx(0, q), std::out_of_range);
        CHECK_THROWS_AS(tableau.Cz(q, 0), std::out_of_range);
        CHECK_THROWS_AS(tableau.Cz(0, q), std::out_of_range);
        CHECK_THROWS_AS(tableau.Cy(q, 0), std::out_of_range);
        CHECK_THROWS_AS(tableau.Cy(0, q), std::out_of_range);
        CHECK_THROWS_AS(tableau.MeasureX(q), std::out_of_range);
        CHECK_THROWS_AS(tableau.MeasureZ(q), std::out_of_range);
        CHECK_THROWS_AS(tableau.IsDeterministicMeasureX(q), std::out_of_range);
        CHECK_THROWS_AS(tableau.IsDeterministicMeasureZ(q), std::out_of_range);
        CHECK_THROWS_AS(tableau.ForceMeasureX(q, false), std::out_of_range);
        CHECK_THROWS_AS(tableau.ForceMeasureZ(q, true), std::out_of_range);
    }
    for (QubitSize q{0}; q < n; ++q) CHECK_FALSE(tableau.MeasureZ(q));
}

TEST_CASE("Simulator validation follows growth and copying", "[simulator][bounds][copy]") {
    Simulator<> simulator{0};
    CHECK(simulator.NumQubits() == 0);
    CHECK_THROWS_AS(simulator.X(0), std::out_of_range);
    simulator.Add();
    CHECK(simulator.NumQubits() == 1);
    simulator.X(0);
    CHECK(simulator.MeasureZ(0));
    CHECK_THROWS_AS(simulator.X(1), std::out_of_range);
    auto copy{simulator};
    simulator.Add();
    CHECK_FALSE(simulator.MeasureZ(1));
    CHECK_THROWS_AS(copy.MeasureZ(1), std::out_of_range);
    copy = simulator;
    CHECK_FALSE(copy.MeasureZ(1));
}

TEST_CASE("Simulator rejects qubit count overflow before calling the backend",
          "[simulator][bounds]") {
    struct FullBackend : InverseTableau {
        explicit FullBackend(QubitSize /*num_qubits*/) : InverseTableau{0} {}
        [[nodiscard]] QubitSize NumQubits() const noexcept {
            return std::numeric_limits<QubitSize>::max();
        }
    };
    Simulator<FullBackend> simulator{0};
    CHECK_THROWS_AS(simulator.Add(), std::length_error);
}

TEST_CASE("GHZ measurement parity preserves phases across word boundaries",
          "[inverse_tableau][measurement]") {
    const QubitSize n{GENERATE(QubitSize{1}, QubitSize{2}, QubitSize{63}, QubitSize{64},
                               QubitSize{65}, QubitSize{129})};
    const unsigned phase{GENERATE(0U, 1U, 2U, 3U)};
    InverseTableau tableau{n};
    tableau.ApplyH(0);
    for (QubitSize q{1}; q < n; ++q) tableau.ApplyCx(0, q);
    for (unsigned i{0}; i < phase; ++i) tableau.ApplyS(0);
    if (phase % 2 != 0) {
        for (unsigned i{0}; i < 3; ++i) tableau.ApplyS(n - 1);
    }
    bool parity{false};
    for (QubitSize q{0}; q < n; ++q) {
        const bool result{tableau.MeasureX(q)};
        CHECK(tableau.MeasureX(q) == result);
        parity = parity != result;
    }
    CHECK(parity == (phase >= 2));
}

TEST_CASE("IsDeterministicMeasureZ returns true for 0 and 1 ket") {
    InverseTableau tableau{2};
    tableau.ApplyX(1);
    CHECK(tableau.IsDeterministicMeasureZ(0) == true);
    CHECK(tableau.IsDeterministicMeasureZ(1) == true);
}

TEST_CASE("IsDeterministicMeasureX returns true for + and - ket") {
    InverseTableau tableau{2};
    tableau.ApplyX(1);
    tableau.ApplyH(0);
    tableau.ApplyH(1);
    CHECK(tableau.IsDeterministicMeasureX(0) == true);
    CHECK(tableau.IsDeterministicMeasureX(1) == true);
}

TEST_CASE("IsDeterministicMeasureZ works when measuring bell pair") {
    InverseTableau tableau{2};
    tableau.ApplyH(0);
    tableau.ApplyCx(0, 1);
    CHECK_FALSE(tableau.IsDeterministicMeasureZ(0));
    CHECK_FALSE(tableau.IsDeterministicMeasureZ(1));
    tableau.MeasureZ(0);
    CHECK(tableau.IsDeterministicMeasureZ(1));
}

TEST_CASE("ForceMeasureZ works correctly") {
    for (std::size_t i{0}; i < 100; i++) {
        InverseTableau tableau{2};
        tableau.ApplyH(0);
        tableau.ApplyH(1);
        CHECK(tableau.ForceMeasureZ(0, 1) == 1);
        CHECK(tableau.ForceMeasureZ(1, 0) == 0);
    }
}

TEST_CASE("ForceMeasureX works correctly") {
    for (std::size_t i{0}; i < 100; i++) {
        InverseTableau tableau{2};
        CHECK(tableau.ForceMeasureX(0, 1) == 1);
        CHECK(tableau.ForceMeasureX(1, 0) == 0);
    }
}

}  // namespace stabaal::detail
