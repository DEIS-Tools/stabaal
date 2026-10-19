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
    REQUIRE(bits.Size() == n);
    CHECK(bits.Words().size() == (n + 63) / 64);
    CHECK_FALSE(bits.Any());
    for (std::size_t i{0}; i < bits.Size(); ++i) CHECK_FALSE(bits.Get(i));
}

TEST_CASE("DynamicBitset updates match independent bit values", "[bitset]") {
    const auto n{GENERATE(0, 1, 63, 64, 65, 129)};
    DynamicBitset bits{static_cast<std::size_t>(n)};
    std::vector<bool> expected(static_cast<std::size_t>(n), false);
    std::mt19937 rng{42};
    for (std::uint32_t step{0}; step < 300; ++step) {
        const bool value{(rng() & 1U) != 0};
        if (expected.empty() || step % 3 == 0) {
            bits.PushBack(value);
            expected.push_back(value);
        } else {
            const auto idx{rng() % expected.size()};
            bits.Set(idx, value);
            expected[idx] = value;
        }
        REQUIRE(bits.Size() == expected.size());
        bool any{false};
        for (std::size_t i{0}; i < expected.size(); ++i) {
            REQUIRE(bits.Get(i) == expected[i]);
            any = any || expected[i];
        }
        CHECK(bits.Any() == any);
        const auto used{bits.Size() % 64};
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
    REQUIRE(moved.Size() == 65);
    CHECK(moved.Get(64));
    CHECK(original.Size() == 0);
    original.PushBack(true);
    REQUIRE(original.Size() == 1);
    CHECK(original.Get(0));

    copy = std::move(moved);
    CHECK(copy.Get(64));
    CHECK(moved.Size() == 0);
    moved.PushBack();
    CHECK_FALSE(moved.Any());
    moved = copy;
    moved.Set(64, false);
    CHECK(copy.Get(64));
}
