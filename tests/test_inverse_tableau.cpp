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

TEST_CASE("InverseTableau copies remain independent after growth", "[inverse_tableau][copy]") {
    const auto n{GENERATE(QubitIdx{0}, QubitIdx{64}, QubitIdx{65})};
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
        GENERATE(QubitIdx{0}, QubitIdx{32}, QubitIdx{63}, QubitIdx{64}, QubitIdx{128})};
    InverseTableau original{initial};
    for (QubitIdx q{0}; q < initial + 145; ++q) {
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
        if (q % 2 != 0)
            original.ApplyZ(q);
        else
            original.ApplyX(q);
        const auto outcome{q % 2 != 0 ? original.MeasureX(q) : original.MeasureZ(q)};
        CHECK(outcome != (q % 3 == 0));
        for (QubitIdx i{0}; i <= q; ++i) {
            CHECK((i % 2 != 0 ? copy.MeasureX(i) : copy.MeasureZ(i)) == (i % 3 == 0));
            CHECK((i % 2 != 0 ? assigned.MeasureX(i) : assigned.MeasureZ(i)) == (i % 3 == 0));
        }
        if (q % 2 != 0)
            original.ApplyZ(q);
        else
            original.ApplyX(q);
        copy.Add();
        CHECK_FALSE(copy.MeasureZ(q + 1));
    }
}

TEST_CASE("Entanglement survives repeated storage growth across word boundaries",
          "[inverse_tableau][measurement][growth]") {
    const auto n{GENERATE(QubitIdx{65}, QubitIdx{129}, QubitIdx{257})};
    InverseTableau tableau{0};
    tableau.Add();
    tableau.ApplyH(0);
    for (QubitIdx q{1}; q < n; ++q) {
        tableau.Add();
        tableau.ApplyCx(0, q);
    }
    auto z_copy{tableau};
    const auto outcome{z_copy.MeasureZ(0)};
    for (QubitIdx q{1}; q < n; ++q) CHECK(z_copy.MeasureZ(q) == outcome);
    bool parity{false};
    for (QubitIdx q{0}; q < n; ++q) parity = parity != tableau.MeasureX(q);
    CHECK_FALSE(parity);
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

TEST_CASE("GHZ measurement parity preserves phases across word boundaries",
          "[inverse_tableau][measurement]") {
    const QubitIdx n{GENERATE(QubitIdx{1}, QubitIdx{2}, QubitIdx{63}, QubitIdx{64}, QubitIdx{65},
                              QubitIdx{129})};
    const unsigned phase{GENERATE(0U, 1U, 2U, 3U)};
    InverseTableau tableau{n};
    tableau.ApplyH(0);
    for (QubitIdx q{1}; q < n; ++q) tableau.ApplyCx(0, q);
    for (unsigned i{0}; i < phase; ++i) tableau.ApplyS(0);
    if (phase % 2 != 0) {
        for (unsigned i{0}; i < 3; ++i) tableau.ApplyS(n - 1);
    }
    bool parity{false};
    for (QubitIdx q{0}; q < n; ++q) {
        const bool result{tableau.MeasureX(q)};
        CHECK(tableau.MeasureX(q) == result);
        parity = parity != result;
    }
    CHECK(parity == (phase >= 2));
}

}  // namespace stabaal::detail
