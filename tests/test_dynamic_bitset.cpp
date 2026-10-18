#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators.hpp>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <random>
#include <stdexcept>
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

TEST_CASE("DynamicBitset XOR handles partial words and self aliasing", "[bitset]") {
    const auto n{GENERATE(0, 1, 63, 64, 65, 129)};
    DynamicBitset a{static_cast<std::size_t>(n)};
    DynamicBitset b{static_cast<std::size_t>(n)};
    for (std::size_t i{0}; i < a.Size(); ++i) {
        a.Set(i, i % 2 == 0);
        b.Set(i, i % 3 == 0);
    }
    a ^= b;
    for (std::size_t i{0}; i < a.Size(); ++i) {
        CHECK(a.Get(i) == ((i % 2 == 0) != (i % 3 == 0)));
        CHECK(b.Get(i) == (i % 3 == 0));
    }
    a ^= a;
    CHECK_FALSE(a.Any());
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

TEST_CASE("DynamicBitset rejects indices outside its logical size", "[bitset][bounds]") {
    const auto n{GENERATE(std::size_t{0}, std::size_t{1}, std::size_t{63}, std::size_t{64},
                          std::size_t{65})};
    DynamicBitset bits{n};
    for (const auto idx : {n, n + 1, std::numeric_limits<std::size_t>::max()}) {
        CHECK_THROWS_AS(bits.Get(idx), std::out_of_range);
        CHECK_THROWS_AS(bits.Set(idx), std::out_of_range);
        CHECK_FALSE(bits.Any());
        CHECK(bits.Size() == n);
    }
}
