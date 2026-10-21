#include <array>
#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators.hpp>
#include <cstddef>
#include <cstdint>
#include <random>
#include <vector>

#include "stabaal/detail/bit_word.hpp"
#include "stabaal/detail/pauli_string.hpp"

using stabaal::detail::kBitsPerWord;
using stabaal::detail::PauliStringView;
using stabaal::detail::WordCount;

namespace {

struct PauliBuffer {
    explicit PauliBuffer(const std::size_t num_qubits)
        : xs(WordCount(num_qubits)), zs(WordCount(num_qubits)) {}

    [[nodiscard]] PauliStringView View() noexcept { return {xs, zs}; }

    std::vector<std::uint64_t> xs;
    std::vector<std::uint64_t> zs;
};

}  // namespace

TEST_CASE("PauliStringView products agree with single-qubit Pauli multiplication", "[pauli]") {
    const auto n{GENERATE(0U, 1U, 63U, 64U, 65U, 127U, 128U, 129U, 10000U)};
    constexpr std::array<std::array<int, 4>, 4> phases{
        {{0, 0, 0, 0}, {0, 0, 3, 1}, {0, 1, 0, 3}, {0, 3, 1, 0}}};
    std::mt19937 rng{42};
    for (std::uint32_t trial{0}; trial < 4; ++trial) {
        PauliBuffer a_buffer{n};
        auto a{a_buffer.View()};
        PauliBuffer b_buffer{n};
        auto b{b_buffer.View()};
        std::vector<unsigned> expected(n);
        int phase{0};
        bool any_x{false};
        for (std::uint32_t q{0}; q < n; ++q) {
            const auto left{static_cast<unsigned>(rng() % 4)};
            const auto right{static_cast<unsigned>(rng() % 4)};
            a.SetX(q, (left & 1U) != 0);
            a.SetZ(q, (left & 2U) != 0);
            b.SetX(q, (right & 1U) != 0);
            b.SetZ(q, (right & 2U) != 0);
            expected[q] = left ^ right;
            phase = (phase + phases[left][right]) % 4;
            any_x = any_x || (expected[q] & 1U) != 0;
        }
        const auto exponent{a.XorWithPhase(b, n)};
        CHECK((exponent % 4 + 4) % 4 == phase);
        CHECK(a.AnyX(n) == any_x);
        for (std::uint32_t q{0}; q < n; ++q) {
            CHECK(a.GetX(q) == ((expected[q] & 1U) != 0));
            CHECK(a.GetZ(q) == ((expected[q] & 2U) != 0));
        }
        CHECK(a.XorWithPhase(a, n) % 4 == 0);
        CHECK_FALSE(a.AnyX(n));
        for (std::uint32_t q{0}; q < n; ++q) CHECK_FALSE(a.GetZ(q));
    }
}

TEST_CASE("Word cleanup matches ordered Cx and S dagger updates", "[pauli][measurement]") {
    const auto n{GENERATE(1U, 2U, 63U, 64U, 65U, 127U, 128U, 129U, 256U, 1024U, 10000U)};
    std::mt19937 rng{1234};
    for (unsigned trial{0}; trial < 16; ++trial) {
        const auto pivot{static_cast<std::uint32_t>(rng() % n)};
        PauliBuffer measured_buffer{n};
        auto measured{measured_buffer.View()};
        PauliBuffer row_buffer{n};
        auto row{row_buffer.View()};
        std::vector<bool> expected_x(n);
        std::vector<bool> expected_z(n);
        std::vector<std::uint64_t> targets(WordCount(n));
        for (std::uint32_t q{0}; q < n; ++q) {
            measured.SetX(q, q == pivot || (q > pivot && (rng() & 1U) != 0));
            measured.SetZ(q, (rng() & 1U) != 0);
            expected_x[q] = (rng() & 1U) != 0;
            expected_z[q] = (rng() & 1U) != 0;
            row.SetX(q, expected_x[q]);
            row.SetZ(q, expected_z[q]);
        }
        bool apply_s_dagger{measured.GetZ(pivot)};
        bool flip{false};
        const bool x{expected_x[pivot]};
        bool z{expected_z[pivot]};
        for (std::uint32_t q{pivot + 1}; q < n; ++q) {
            if (!measured.GetX(q)) continue;
            apply_s_dagger = apply_s_dagger != measured.GetZ(q);
            flip = flip != (x && expected_z[q] && (expected_x[q] == z));
            expected_x[q] = expected_x[q] != x;
            z = z != expected_z[q];
        }
        if (apply_s_dagger) {
            flip = flip != (x && !z);
            z = z != x;
        }
        expected_z[pivot] = z;
        CHECK(measured.MeasurementTargets(pivot, targets) == apply_s_dagger);
        for (std::uint32_t q{0}; q < n; ++q) {
            CHECK(((targets[q / kBitsPerWord] >> (q % kBitsPerWord)) & 1U) ==
                  static_cast<unsigned>(q != pivot && measured.GetX(q)));
        }
        CHECK(row.CleanUp(pivot, targets, apply_s_dagger) == flip);
        for (std::uint32_t q{0}; q < n; ++q) {
            CHECK(row.GetX(q) == expected_x[q]);
            CHECK(row.GetZ(q) == expected_z[q]);
        }
    }
}
