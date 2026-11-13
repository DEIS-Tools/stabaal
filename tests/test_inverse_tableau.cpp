#include <algorithm>
#include <array>
#include <bit>
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

#include "stabaal/backends/inverse_tableau.hpp"
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

    [[nodiscard]] double PauliExpectation(const std::size_t xs, const std::size_t zs) const {
        constexpr std::array<std::complex<double>, 4> phases{{1, {0, 1}, -1, {0, -1}}};
        const auto phase{phases[std::popcount(xs & zs) % 4]};
        std::complex<double> expectation{};
        for (std::size_t i{0}; i < amplitudes.size(); ++i) {
            const auto sign{std::popcount(i & zs) % 2 == 0 ? 1.0 : -1.0};
            expectation += std::conj(amplitudes[i ^ xs]) * phase * sign * amplitudes[i];
        }
        CHECK(std::abs(expectation.imag()) < 1e-9);
        return expectation.real();
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

TEST_CASE("Measurement RNG belongs to each tableau and is preserved by copies",
          "[inverse_tableau][random][copy]") {
    InverseTableau original{1};
    for (unsigned i{0}; i < 11; ++i) {
        original.ApplyH(0);
        original.MeasureZ(0);
    }

    auto copied{original};
    InverseTableau assigned{1};
    assigned = original;
    InverseTableau resized{0};
    resized = original;
    auto grown{original};
    const auto added{grown.Add()};
    grown.Discard(added);
    auto moved{std::move(copied)};
    InverseTableau unrelated{1};

    for (unsigned i{0}; i < 256; ++i) {
        original.ApplyH(0);
        const auto expected{original.MeasureZ(0)};
        for (auto* tableau : {&moved, &assigned, &resized, &grown}) {
            unrelated.ApplyH(0);
            unrelated.MeasureZ(0);
            tableau->ApplyH(0);
            CHECK(tableau->MeasureZ(0) == expected);
        }
    }
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
              "\nX0 [a]: | 1 | 0 | 0 |"
              "\nZ0 [a]: | 0 | 1 | 0 |"
              "\n");
    }
    SECTION("fresh 2-qubit tableau") {
        const InverseTableau t{2};
        CoutCapture cap{};
        t.Print();
        CHECK(cap.ss.str() ==
              "\n    | x x | z z | r |"
              "\nX0 [a]: | 1 0 | 0 0 | 0 |"
              "\nX1 [a]: | 0 1 | 0 0 | 0 |"
              "\nZ0 [a]: | 0 0 | 1 0 | 0 |"
              "\nZ1 [a]: | 0 0 | 0 1 | 0 |"
              "\n");
    }
    SECTION("sign column after X") {
        InverseTableau t{1};
        t.ApplyX(0);
        CoutCapture cap{};
        t.Print();
        CHECK(cap.ss.str() ==
              "\n    | x | z | r |"
              "\nX0 [a]: | 1 | 0 | 0 |"
              "\nZ0 [a]: | 0 | 1 | 1 |"
              "\n");
    }
    SECTION("rows swapped after H") {
        InverseTableau t{1};
        t.ApplyH(0);
        CoutCapture cap{};
        t.Print();
        CHECK(cap.ss.str() ==
              "\n    | x | z | r |"
              "\nX0 [a]: | 0 | 1 | 0 |"
              "\nZ0 [a]: | 1 | 0 | 0 |"
              "\n");
    }
    SECTION("discarded slots show their status") {
        InverseTableau t{4};
        t.ApplyH(0);
        t.ApplyCx(0, 1);
        t.Discard(1);
        t.Discard(3);
        CoutCapture cap{};
        t.Print();
        const auto output{cap.ss.str()};
        CHECK(output.find("X0 [a]:") != std::string::npos);
        CHECK(output.find("X1 [p]:") != std::string::npos);
        CHECK(output.find("X3 [†]:") != std::string::npos);
        CHECK(output.find("Z1 [p]:") != std::string::npos);
        CHECK(output.find("Z3 [†]:") != std::string::npos);
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
    CHECK(std::as_const(tableau).IsDeterministicMeasureZ(0) == true);
    CHECK(std::as_const(tableau).IsDeterministicMeasureZ(1) == true);
}

TEST_CASE("IsDeterministicMeasureX returns true for + and - ket") {
    InverseTableau tableau{2};
    tableau.ApplyX(1);
    tableau.ApplyH(0);
    tableau.ApplyH(1);
    CHECK(std::as_const(tableau).IsDeterministicMeasureX(0) == true);
    CHECK(std::as_const(tableau).IsDeterministicMeasureX(1) == true);
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

TEST_CASE("Single-qubit Pauli eigenstates are unentangled", "[simulator][entanglement]") {
    const auto n{GENERATE(QubitSize{1}, QubitSize{64}, QubitSize{65}, QubitSize{129})};
    const auto basis{GENERATE(0, 1, 2)};  // Z, X, Y
    const bool negative{GENERATE(false, true)};
    Simulator<> sim{n};
    const auto q{n - 1};
    if (n > 2) {
        sim.H(0);
        sim.Cx(0, 1);
    }
    if (negative) sim.X(q);
    if (basis != 0) sim.H(q);
    if (basis == 2) sim.S(q);
    if (basis == 2) {
        CHECK_FALSE(sim.IsDeterministicMeasureX(q));
        CHECK_FALSE(sim.IsDeterministicMeasureZ(q));
    }
    if (basis == 2) sim.Sdg(q);
    if (basis != 0) sim.H(q);
    CHECK(sim.IsDeterministicMeasureZ(q));
    CHECK(sim.MeasureZ(q) == negative);
}

TEST_CASE("Bell pairs remain intact across word boundaries", "[simulator][entanglement]") {
    const auto n{GENERATE(QubitSize{2}, QubitSize{64}, QubitSize{65}, QubitSize{129})};
    const auto rotation{GENERATE(0, 1, 2)};
    const auto target{n - 1};
    Simulator<> sim{n};
    sim.H(0);
    sim.Cx(0, target);
    if (rotation == 1) {
        sim.H(0);
        sim.H(target);
    } else if (rotation == 2) {
        sim.S(0);
        sim.Sdg(target);
    }
    const auto memory{sim.MemoryUsage()};
    for (unsigned repeat{0}; repeat < 3; ++repeat) {
        CHECK_FALSE(sim.IsDeterministicMeasureZ(0));
        CHECK_FALSE(sim.IsDeterministicMeasureZ(target));
    }
    CHECK(sim.MemoryUsage() == memory);
    if (rotation == 1) {
        sim.H(0);
        sim.H(target);
    } else if (rotation == 2) {
        sim.Sdg(0);
        sim.S(target);
    }
    sim.Cx(0, target);
    sim.H(0);
    CHECK(sim.IsDeterministicMeasureZ(0));
    CHECK(sim.IsDeterministicMeasureZ(target));
    CHECK_FALSE(sim.MeasureZ(0));
    CHECK_FALSE(sim.MeasureZ(target));
}

TEST_CASE("GHZ measurements remain correlated after collapse", "[simulator][entanglement]") {
    Simulator<> sim{4};
    sim.H(0);
    sim.Cx(0, 1);
    sim.Cx(0, 2);
    const auto outcome{sim.MeasureZ(0)};
    CHECK(sim.MeasureZ(1) == outcome);
    CHECK(sim.MeasureZ(2) == outcome);
}

TEST_CASE("Bell pairs using added qubits survive growth and copying", "[simulator][entanglement]") {
    const auto initial{GENERATE(QubitSize{0}, QubitSize{63}, QubitSize{64}, QubitSize{128})};
    Simulator<> sim{initial};
    const auto first{sim.Add()};
    const auto second{sim.Add()};
    sim.H(first);
    sim.Cx(first, second);
    for (unsigned i{0}; i < 150; ++i) sim.Add();
    auto copy{sim};
    copy.Cx(first, second);
    copy.H(first);
    CHECK_FALSE(copy.MeasureZ(first));
    CHECK_FALSE(copy.MeasureZ(second));
    CHECK_FALSE(sim.IsDeterministicMeasureZ(first));
    CHECK_FALSE(sim.IsDeterministicMeasureZ(second));
}

TEST_CASE("Discard resets every pure product eigenstate for storage reuse",
          "[simulator][discard]") {
    const auto n{
        GENERATE(QubitSize{1}, QubitSize{3}, QubitSize{64}, QubitSize{65}, QubitSize{129})};
    const auto basis{GENERATE(0, 1, 2)};  // Z, X, Y
    const bool negative{GENERATE(false, true)};
    const auto q{n - 1};
    Simulator<> sim{n};
    if (n > 2) {
        sim.H(0);
        sim.Cx(0, 1);
    }
    if (negative) sim.X(q);
    if (basis != 0) sim.H(q);
    if (basis == 2) sim.S(q);
    sim.Discard(q);
    CHECK(sim.NumQubits() == n - 1);
    auto copy{sim};
    const auto memory{sim.MemoryUsage()};
    const auto added{sim.Add()};
    CHECK(added != q);
    CHECK(sim.NumQubits() == n);
    if (n > 1) CHECK(sim.MemoryUsage() == memory + std::size_t{n} * sizeof(std::uint32_t));
    CHECK(sim.IsDeterministicMeasureZ(added));
    CHECK_FALSE(sim.MeasureZ(added));
    CHECK(copy.Add() == added);
    CHECK_FALSE(copy.MeasureZ(added));
    if (n > 2) {
        sim.Cx(0, 1);
        sim.H(0);
        CHECK(sim.IsDeterministicMeasureZ(0));
        CHECK(sim.IsDeterministicMeasureZ(1));
        CHECK_FALSE(sim.MeasureZ(0));
        CHECK_FALSE(sim.MeasureZ(1));
    }
}

TEST_CASE("Discarding a Bell partner produces an exact maximally mixed qubit",
          "[simulator][discard][mixed]") {
    const auto n{GENERATE(QubitSize{2}, QubitSize{64}, QubitSize{65}, QubitSize{129})};
    Simulator<> sim{n};
    sim.H(0);
    sim.Cx(0, n - 1);
    sim.Discard(n - 1);
    CHECK(sim.NumQubits() == n - 1);
    CHECK_FALSE(sim.IsDeterministicMeasureX(0));
    CHECK_FALSE(sim.IsDeterministicMeasureZ(0));
    auto copy{sim};
    const auto added{sim.Add()};
    CHECK(added != n - 1);
    CHECK(sim.IsDeterministicMeasureZ(added));
    CHECK_FALSE(sim.MeasureZ(added));
    sim.H(0);
    sim.S(0);
    CHECK_FALSE(sim.IsDeterministicMeasureX(0));
    CHECK_FALSE(sim.IsDeterministicMeasureZ(0));
    sim.Cx(0, added);
    // The reduced state is a classical mixture of |00> and |11>.
    const auto outcome{sim.MeasureZ(0)};
    CHECK(sim.IsDeterministicMeasureZ(added));
    CHECK(sim.MeasureZ(added) == outcome);
    CHECK_FALSE(copy.IsDeterministicMeasureZ(0));
    copy.Discard(0);  // Locally mixed but not entangled with another active qubit.
    const auto replacement{copy.Add()};
    CHECK(copy.IsDeterministicMeasureZ(replacement));
    CHECK_FALSE(copy.MeasureZ(replacement));
}

TEST_CASE("A partial trace distinguishes classical GHZ correlations from entanglement",
          "[simulator][discard][mixed]") {
    Simulator<> sim{4};
    sim.H(0);
    sim.Cx(0, 1);
    sim.Cx(0, 3);
    sim.Discard(3);
    for (const QubitId q : {0, 1}) {
        CHECK_FALSE(sim.IsDeterministicMeasureX(q));
        CHECK_FALSE(sim.IsDeterministicMeasureZ(q));
    }
    auto zz{sim};
    zz.Cx(0, 1);
    CHECK(zz.IsDeterministicMeasureZ(1));
    CHECK_FALSE(zz.MeasureZ(1));
    auto xx{sim};
    xx.H(0);
    xx.H(1);
    xx.Cx(0, 1);
    CHECK_FALSE(xx.IsDeterministicMeasureZ(1));
    // A remaining product qubit can still be reused without changing the mixed subsystem.
    sim.Discard(2);
    const auto memory{sim.MemoryUsage()};
    const auto added{sim.Add()};
    CHECK(sim.MemoryUsage() <= memory + 4 * sizeof(std::uint32_t));
    CHECK_FALSE(sim.MeasureZ(added));
    CHECK_FALSE(sim.IsDeterministicMeasureZ(1));
}

TEST_CASE("Partial trace preserves independent Bell entanglement through growth and copies",
          "[simulator][discard][mixed][copy]") {
    Simulator<> sim{4};
    sim.H(0);
    sim.Cx(0, 1);
    sim.H(2);
    sim.Cx(2, 3);
    sim.Discard(3);
    std::vector<QubitId> added;
    for (unsigned i{0}; i < 70; ++i) added.push_back(sim.Add());
    CHECK_FALSE(sim.IsDeterministicMeasureZ(2));
    auto copy{sim};
    sim.Discard(0);
    CHECK_FALSE(sim.IsDeterministicMeasureZ(1));
    copy.Cx(0, 1);
    copy.H(0);
    CHECK_FALSE(copy.MeasureZ(0));
    CHECK_FALSE(copy.MeasureZ(1));
    for (const auto id : added) CHECK_FALSE(copy.MeasureZ(id));

    InverseTableau assigned{2};
    assigned.ApplyH(0);
    assigned.ApplyCx(0, 1);
    assigned.Discard(1);
    auto other{assigned};
    other = assigned;
    auto moved{std::move(other)};
    CHECK_FALSE(moved.IsDeterministicMeasureZ(0));
    moved.Discard(0);
    CHECK(moved.NumQubits() == 0);
    auto empty_copy{moved};
    const auto fresh{empty_copy.Add()};
    CHECK_FALSE(empty_copy.MeasureZ(fresh));
    CHECK(empty_copy.IsDeterministicMeasureZ(fresh));
}

TEST_CASE("Mixed entanglement reduction handles hidden slots across multiple words",
          "[simulator][discard][mixed]") {
    const auto n{GENERATE(QubitSize{65}, QubitSize{129}, QubitSize{257})};
    const auto half{n / 2};
    Simulator<> sim{n};
    sim.H(0);
    sim.Cx(0, 1);
    for (QubitSize q{2}; q < half; ++q) {
        sim.H(q);
        sim.Cx(q, half + q);
    }
    for (QubitSize q{half}; q < n; ++q) sim.Discard(q);
    CHECK(sim.NumQubits() == half);
    CHECK_FALSE(sim.IsDeterministicMeasureX(half - 1));
    CHECK_FALSE(sim.IsDeterministicMeasureZ(half - 1));
    auto copy{sim};
    copy.Cx(0, 1);
    copy.H(0);
    CHECK_FALSE(copy.MeasureZ(0));
    CHECK_FALSE(copy.MeasureZ(1));
    CHECK_FALSE(copy.IsDeterministicMeasureZ(half - 1));
    sim.Discard(0);
    const auto fresh{sim.Add()};
    CHECK(sim.IsDeterministicMeasureZ(fresh));
    CHECK_FALSE(sim.MeasureZ(fresh));
    CHECK_FALSE(sim.IsDeterministicMeasureZ(1));
}

struct ActiveReferenceQubit {
    QubitId id;
    QubitSize reference_index;
};

// Compare every Pauli expectation, which completely specifies the active reduced density matrix.
// The reference never resets or reuses discarded slots; tracing them out means ignoring them.
void CheckReducedState(const InverseTableau& tableau, const StateVector& reference,
                       const std::vector<ActiveReferenceQubit>& active) {
    REQUIRE(tableau.NumQubits() == active.size());
    const auto pauli_count{std::size_t{1} << (2 * active.size())};
    for (std::size_t pauli{1}; pauli < pauli_count; ++pauli) {
        auto measured{tableau};
        std::size_t xs{0};
        std::size_t zs{0};
        std::size_t target{0};
        for (std::size_t q{0}; q < active.size(); ++q) {
            const auto type{(pauli >> (2 * q)) & 3U};  // I, X, Y, Z
            const auto bit{std::size_t{1} << active[q].reference_index};
            if (type == 1 || type == 2) xs |= bit;
            if (type == 2 || type == 3) zs |= bit;
            if (type == 0) continue;
            target = q;
            if (type == 2) measured.ApplySdg(active[q].id);
            if (type != 3) measured.ApplyH(active[q].id);
        }
        for (std::size_t q{0}; q < active.size(); ++q) {
            if (q != target && ((pauli >> (2 * q)) & 3U) != 0) {
                measured.ApplyCx(active[q].id, active[target].id);
            }
        }
        const auto expectation{reference.PauliExpectation(xs, zs)};
        const bool deterministic{std::abs(expectation) > 1e-9};
        CAPTURE(pauli, expectation);
        REQUIRE(measured.IsDeterministicMeasureZ(active[target].id) == deterministic);
        if (deterministic) {
            CHECK(std::abs(std::abs(expectation) - 1) < 1e-9);
            CHECK(measured.MeasureZ(active[target].id) == (expectation < 0));
        }
    }
}

TEST_CASE("Exact mixed states agree with state-vector partial traces during id reuse",
          "[inverse_tableau][discard][mixed][reference]") {
    std::mt19937 choices{71421};
    constexpr std::array<InverseTableau::UnaryGateFunction, 7> unary{
        &InverseTableau::ApplyH, &InverseTableau::ApplyX, &InverseTableau::ApplyS,
        &InverseTableau::ApplyZ, &InverseTableau::ApplyY, &InverseTableau::ApplySdg,
        &InverseTableau::ApplyI};
    constexpr std::array<InverseTableau::BinaryGateFunction, 3> binary{
        &InverseTableau::ApplyCx, &InverseTableau::ApplyCz, &InverseTableau::ApplyCy};
    for (unsigned trial{0}; trial < 24; ++trial) {
        InverseTableau tableau{3};
        StateVector reference{3};
        QubitSize reference_size{3};
        std::vector<ActiveReferenceQubit> active{{0, 0}, {1, 1}, {2, 2}};
        tableau.ApplyH(0);
        tableau.ApplyCx(0, 1);
        reference.Apply(0, 0, 0);
        reference.Apply(7, 0, 1);
        for (unsigned step{0}; step < 40; ++step) {
            CAPTURE(trial, step);
            if (active.empty() && reference_size == 7) break;
            const auto action{choices() % 6};
            if (active.empty() || (action == 0 && active.size() < 3 && reference_size < 7)) {
                active.push_back({tableau.Add(), reference_size++});
                reference.Add();
            } else {
                const auto q{choices() % active.size()};
                const auto [id, ref]{active[q]};
                if (action == 1) {
                    tableau.Discard(id);
                    active.erase(active.begin() + static_cast<std::ptrdiff_t>(q));
                    CHECK_THROWS_AS(tableau.ApplyH(id), std::out_of_range);
                } else if (action == 2) {
                    const bool x_basis{(choices() & 1U) != 0};
                    if (x_basis) reference.Apply(0, ref, 0);
                    reference.MeasureZ(ref, x_basis ? tableau.MeasureX(id) : tableau.MeasureZ(id));
                    if (x_basis) reference.Apply(0, ref, 0);
                } else if (action == 3 && active.size() > 1) {
                    const auto t{(q + 1 + choices() % (active.size() - 1)) % active.size()};
                    const auto gate{choices() % binary.size()};
                    (tableau.*binary[gate])(id, active[t].id);
                    reference.Apply(static_cast<std::uint8_t>(7 + gate), ref,
                                    active[t].reference_index);
                } else {
                    const auto gate{choices() % unary.size()};
                    (tableau.*unary[gate])(id);
                    reference.Apply(static_cast<std::uint8_t>(gate), ref, 0);
                }
            }
            CheckReducedState(tableau, reference, active);
        }
    }
}

TEST_CASE("Geometric contraction preserves dense inverse frames and stable IDs",
          "[inverse_tableau][discard][contract][reference]") {
    std::mt19937 choices{92183};
    for (unsigned trial{0}; trial < 32; ++trial) {
        InverseTableau tableau{6};
        StateVector reference{6};
        std::vector<ActiveReferenceQubit> active{{0, 0}, {1, 1}, {2, 2}, {3, 3}, {4, 4}, {5, 5}};
        for (unsigned step{0}; step < 50; ++step) {
            const auto q{static_cast<QubitSize>(choices() % 6)};
            const auto t{static_cast<QubitSize>((q + 1 + choices() % 5) % 6)};
            const auto gate{choices() % 3};
            if (gate == 0) {
                tableau.ApplyH(q);
                reference.Apply(0, q, 0);
            } else if (gate == 1) {
                tableau.ApplyS(q);
                reference.Apply(2, q, 0);
            } else {
                tableau.ApplyCx(q, t);
                reference.Apply(7, q, t);
            }
        }
        for (const QubitId id : {1, 3, 4}) {
            const auto outcome{tableau.MeasureZ(id)};
            reference.MeasureZ(static_cast<QubitSize>(id), outcome);
            tableau.Discard(id);
            std::erase_if(active, [id](const auto& q) { return q.id == id; });
            CHECK_THROWS_AS(tableau.ApplyH(id), std::out_of_range);
        }
        CAPTURE(trial);
        CheckReducedState(tableau, reference, active);
        const auto id{tableau.Add()};
        reference.Add();
        active.push_back({id, 6});
        tableau.ApplyH(id);
        reference.Apply(0, 6, 0);
        tableau.ApplyCx(id, active.front().id);
        reference.Apply(7, 6, active.front().reference_index);
        CheckReducedState(tableau, reference, active);
        auto copy{tableau};
        CheckReducedState(copy, reference, active);
    }
}

TEST_CASE("Hidden purification agrees with exact partial traces",
          "[inverse_tableau][discard][mixed][contract][reference]") {
    std::mt19937 choices{28314};
    for (unsigned trial{0}; trial < 24; ++trial) {
        InverseTableau tableau{7};
        StateVector reference{7};
        std::vector<ActiveReferenceQubit> active{{0, 0}, {1, 1}, {2, 2}, {3, 3},
                                                 {4, 4}, {5, 5}, {6, 6}};
        // A connected graph state guarantees every initial slot is locally mixed.
        for (QubitId q{0}; q < 7; ++q) {
            tableau.ApplyH(q);
            reference.Apply(0, q, 0);
        }
        for (QubitId q{0}; q < 6; ++q) {
            tableau.ApplyCz(q, q + 1);
            reference.Apply(8, q, q + 1);
        }
        for (unsigned step{0}; step < 40; ++step) {
            const auto q{static_cast<QubitSize>(choices() % 7)};
            const auto t{static_cast<QubitSize>((q + 1 + choices() % 6) % 7)};
            if (choices() % 2 == 0) {
                tableau.ApplyS(q);
                reference.Apply(2, q, 0);
            } else {
                tableau.ApplyCx(q, t);
                reference.Apply(7, q, t);
            }
        }
        for (const QubitId id : {0, 2, 3, 4, 6}) {
            tableau.Discard(id);
            std::erase_if(active, [id](const auto& q) { return q.id == id; });
            CHECK_THROWS_AS(tableau.ApplyH(id), std::out_of_range);
        }
        CAPTURE(trial);
        CheckReducedState(tableau, reference, active);
        // Verify the reduced state remains exact as active gates change it.
        for (unsigned step{0}; step < 10; ++step) {
            tableau.ApplyH(active[0].id);
            reference.Apply(0, active[0].reference_index, 0);
            tableau.ApplyCx(active[0].id, active[1].id);
            reference.Apply(7, active[0].reference_index, active[1].reference_index);
            CheckReducedState(tableau, reference, active);
        }
        const auto outcome{tableau.MeasureX(active[0].id)};
        reference.Apply(0, active[0].reference_index, 0);
        reference.MeasureZ(active[0].reference_index, outcome);
        reference.Apply(0, active[0].reference_index, 0);
        CheckReducedState(tableau, reference, active);
    }
}

TEST_CASE("Repeated entangled discard preserves exact mixed states",
          "[simulator][discard][mixed][contract][memory]") {
    Simulator<> sim{2};
    QubitId first{0};
    QubitId second{1};
    std::size_t largest{sim.MemoryUsage()};
    sim.H(first);
    sim.Cx(first, second);
    sim.Discard(second);
    second = sim.Add();
    for (unsigned step{0}; step < 100; ++step) {
        // Transfer maximal local mixedness to the new qubit, without measuring either partner.
        sim.H(second);
        sim.Cx(second, first);
        const auto old{first};
        sim.Discard(first);
        first = second;
        second = sim.Add();
        CHECK_FALSE(sim.IsDeterministicMeasureX(first));
        CHECK_FALSE(sim.IsDeterministicMeasureZ(first));
        CHECK_THROWS_AS(sim.H(old), std::out_of_range);
        CHECK(sim.IsDeterministicMeasureZ(second));
        largest = std::max(largest, sim.MemoryUsage());
    }
    CHECK(largest > 4096);
}

TEST_CASE("Contraction triggers geometrically and preserves high-valued live IDs",
          "[simulator][discard][contract][memory]") {
    Simulator<> sim{128};
    sim.H(126);
    sim.Cx(126, 127);
    const auto memory{sim.MemoryUsage()};
    for (QubitId q{0}; q < 63; ++q) {
        sim.Discard(q);
        CHECK(sim.MemoryUsage() == memory);
    }
    sim.Discard(63);
    CHECK(sim.MemoryUsage() < memory);
    const auto contracted{sim.MemoryUsage()};
    for (QubitId q{64}; q < 95; ++q) {
        sim.Discard(q);
        CHECK(sim.MemoryUsage() == contracted);
    }
    sim.Discard(95);
    CHECK(sim.MemoryUsage() < contracted);
    sim.Cx(126, 127);
    sim.H(126);
    CHECK_FALSE(sim.MeasureZ(126));
    CHECK_FALSE(sim.MeasureZ(127));
    const auto added{sim.Add()};
    CHECK_FALSE(sim.MeasureZ(added));
    CHECK_THROWS_AS(sim.H(63), std::out_of_range);
}

TEST_CASE("Hidden purification preserves independent mixed qubits across word boundaries",
          "[simulator][discard][mixed][contract]") {
    Simulator<> sim{129};
    for (QubitId q{0}; q < 64; ++q) {
        sim.H(q);
        sim.Cx(q, q + 64);
    }
    const auto initial{sim.MemoryUsage()};
    for (QubitId q{64}; q < 129; ++q) sim.Discard(q);
    for (QubitId q{0}; q < 48; ++q) sim.Discard(q);
    CHECK(sim.NumQubits() == 16);
    CHECK(sim.MemoryUsage() == initial);
    for (QubitId q{48}; q < 64; ++q) {
        sim.H(q);
        sim.S(q);
        CHECK_FALSE(sim.IsDeterministicMeasureX(q));
        CHECK_FALSE(sim.IsDeterministicMeasureZ(q));
    }
    for (QubitId q{48}; q < 63; ++q) {
        auto parity{sim};
        parity.Cx(q, q + 1);
        CHECK_FALSE(parity.IsDeterministicMeasureZ(q + 1));
    }
    const auto added{sim.Add()};
    CHECK(sim.IsDeterministicMeasureZ(added));
    sim.Cx(63, added);
    sim.Cx(63, added);
    CHECK(sim.IsDeterministicMeasureZ(added));
    CHECK_FALSE(sim.MeasureZ(added));
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

TEST_CASE("Test gate information") {
    InverseTableau tab{1};
    auto unary_gates = tab.SupportedUnaryGates();
    auto binary_gates = tab.SupportedBinaryGates();
    for (InverseTableau::UnaryGateFunction gate : unary_gates) {
        CHECK(tab.GetGateInformation(gate).arity == 1);
    }
    for (InverseTableau::BinaryGateFunction gate : binary_gates) {
        CHECK(tab.GetGateInformation(gate).arity == 2);
    }
}

}  // namespace stabaal::detail
