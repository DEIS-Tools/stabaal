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
    REQUIRE(bits.size() == n);
    CHECK(bits.words().size() == (n + 63) / 64);
    CHECK_FALSE(bits.any());
    for (std::size_t i{0}; i < bits.size(); ++i) CHECK_FALSE(bits.get(i));
}

TEST_CASE("DynamicBitset updates match independent bit values", "[bitset]") {
    const auto n{GENERATE(0, 1, 63, 64, 65, 129)};
    DynamicBitset bits{static_cast<std::size_t>(n)};
    std::vector<bool> expected(static_cast<std::size_t>(n), false);
    std::mt19937 rng{42};
    for (std::uint32_t step{0}; step < 300; ++step) {
        const bool value{(rng() & 1U) != 0};
        if (expected.empty() || step % 3 == 0) {
            bits.pushBack(value);
            expected.push_back(value);
        } else {
            const auto idx{rng() % expected.size()};
            bits.set(idx, value);
            expected[idx] = value;
        }
        REQUIRE(bits.size() == expected.size());
        bool any{false};
        for (std::size_t i{0}; i < expected.size(); ++i) {
            REQUIRE(bits.get(i) == expected[i]);
            any = any || expected[i];
        }
        CHECK(bits.any() == any);
        const auto used{bits.size() % 64};
        if (used != 0) CHECK((bits.words().back() >> used) == 0);
    }
}

TEST_CASE("DynamicBitset copies own their storage and moved-from objects can grow", "[bitset]") {
    DynamicBitset original{65};
    original.set(64);
    DynamicBitset copy{original};
    copy.set(64, false);
    CHECK(original.get(64));
    CHECK_FALSE(copy.any());

    DynamicBitset moved{std::move(original)};
    REQUIRE(moved.size() == 65);
    CHECK(moved.get(64));
    CHECK(original.size() == 0);
    original.pushBack(true);
    REQUIRE(original.size() == 1);
    CHECK(original.get(0));

    copy = std::move(moved);
    CHECK(copy.get(64));
    CHECK(moved.size() == 0);
    moved.pushBack();
    CHECK_FALSE(moved.any());
    moved = copy;
    moved.set(64, false);
    CHECK(copy.get(64));
}
