#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators.hpp>
#include <cstddef>
#include <cstdint>
#include <random>
#include <utility>
#include <vector>

#include "dynamic_bitset.hpp"

using stabaal::DynamicBitset;

TEST_CASE("DynamicBitset starts zero with packed storage", "[bitset]") {
    const auto n{GENERATE(0, 1, 63, 64, 65, 127, 128, 129, 1024)};
    const DynamicBitset bits{static_cast<std::size_t>(n)};
    CHECK(bits.Words().size() == (n + 63) / 64);
    CHECK_FALSE(bits.Any());
    for (std::size_t i{0}; i < static_cast<std::size_t>(n); ++i) CHECK_FALSE(bits.Get(i));
}

TEST_CASE("DynamicBitset updates match independent bit values", "[bitset]") {
    const auto n{GENERATE(0, 1, 63, 64, 65, 129)};
    DynamicBitset bits{static_cast<std::size_t>(n)};
    std::vector<bool> expected(static_cast<std::size_t>(n), false);
    std::mt19937 rng{42};
    for (std::uint32_t step{0}; step < 300; ++step) {
        const bool value{(rng() & 1U) != 0};
        if (expected.empty() || step % 3 == 0) {
            bits.Resize(expected.size() + 1);
            bits.Set(expected.size(), value);
            expected.push_back(value);
        } else {
            const auto idx{rng() % expected.size()};
            bits.Set(idx, value);
            expected[idx] = value;
        }
        REQUIRE(bits.Words().size() == (expected.size() + 63) / 64);
        bool any{false};
        for (std::size_t i{0}; i < expected.size(); ++i) {
            REQUIRE(bits.Get(i) == expected[i]);
            any = any || expected[i];
        }
        CHECK(bits.Any() == any);
        const auto used{expected.size() % 64};
        if (used != 0) CHECK((bits.Words().back() >> used) == 0);
    }
}

TEST_CASE("DynamicBitset copies own their storage and moved-from objects can grow", "[bitset]") {
    DynamicBitset original{65};
    original.Set(64);
    DynamicBitset copy{original};
    copy.Set(64, false);
    CHECK(original.Get(64));
    CHECK_FALSE(copy.Any());

    DynamicBitset moved{std::move(original)};
    REQUIRE(moved.Words().size() == 2);
    CHECK(moved.Get(64));
    original.Resize(0);
    original.Resize(1);
    original.Set(0);
    REQUIRE(original.Words().size() == 1);
    CHECK(original.Get(0));

    copy = std::move(moved);
    CHECK(copy.Get(64));
    moved.Resize(0);
    moved.Resize(1);
    CHECK_FALSE(moved.Any());
    moved = copy;
    moved.Set(64, false);
    CHECK(copy.Get(64));
}

TEST_CASE("DynamicBitset growth preserves bits and zeroes new words", "[bitset]") {
    DynamicBitset bits{1};
    bits.Set(0);
    bits.Resize(129);
    CHECK(bits.Get(0));
    CHECK_FALSE(bits.Get(64));
    CHECK_FALSE(bits.Get(128));
    bits.Set(128);
    bits.Resize(64);
    bits.Resize(129);
    CHECK(bits.Get(0));
    CHECK_FALSE(bits.Get(128));
}
