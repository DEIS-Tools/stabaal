#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators.hpp>
#include <cstddef>
#include <cstdint>
#include <random>
#include <utility>
#include <vector>

#include "stabaal/detail/dynamic_bitset.hpp"

using stabaal::detail::DynamicBitset;

TEST_CASE("DynamicBitset starts zero with packed storage", "[bitset]") {
    const auto n{GENERATE(0, 1, 63, 64, 65, 127, 128, 129, 193, 257, 513, 1024, 1025)};
    const DynamicBitset bits{static_cast<std::size_t>(n)};
    CHECK(bits.words(n).size() == (n + 63) / 64);
    const auto used_bytes{bits.words(n).size() * sizeof(std::uint64_t)};
    CHECK(DynamicBitset::allocatedBytes(n) >= used_bytes);
    if (n == 0) {
        CHECK(DynamicBitset::allocatedBytes(n) == 0);
    } else {
        CHECK(DynamicBitset::allocatedBytes(n) < 2 * used_bytes);
    }
    CHECK_FALSE(bits.any(n));
    for (std::size_t i{0}; i < static_cast<std::size_t>(n); ++i) CHECK_FALSE(bits.get(i));
}

TEST_CASE("DynamicBitset updates match independent bit values", "[bitset]") {
    const auto n{GENERATE(0, 1, 63, 64, 65, 129)};
    DynamicBitset bits{static_cast<std::size_t>(n)};
    std::vector<bool> expected(static_cast<std::size_t>(n), false);
    std::mt19937 rng{42};
    for (std::uint32_t step{0}; step < 300; ++step) {
        const bool value{(rng() & 1U) != 0};
        if (expected.empty() || step % 3 == 0) {
            bits.resize(expected.size(), expected.size() + 1);
            bits.set(expected.size(), value);
            expected.push_back(value);
        } else {
            const auto idx{rng() % expected.size()};
            bits.set(idx, value);
            expected[idx] = value;
        }
        REQUIRE(bits.words(expected.size()).size() == (expected.size() + 63) / 64);
        bool any{false};
        for (std::size_t i{0}; i < expected.size(); ++i) {
            REQUIRE(bits.get(i) == expected[i]);
            any = any || expected[i];
        }
        CHECK(bits.any(expected.size()) == any);
        const auto used{expected.size() % 64};
        if (used != 0) CHECK((bits.words(expected.size()).back() >> used) == 0);
    }
}

TEST_CASE("DynamicBitset copies own their storage and moved-from objects can grow", "[bitset]") {
    DynamicBitset original{65};
    original.set(64);
    DynamicBitset copy{original, 65};
    copy.set(64, false);
    CHECK(original.get(64));
    CHECK_FALSE(copy.any(65));

    DynamicBitset moved{std::move(original)};
    REQUIRE(moved.words(65).size() == 2);
    CHECK(moved.get(64));
    CHECK(original.words(0).data() == nullptr);
    CHECK(DynamicBitset::allocatedBytes(0) == 0);
    original.resize(0, 1);
    original.set(0);
    REQUIRE(original.words(1).size() == 1);
    CHECK(original.get(0));

    copy = std::move(moved);
    CHECK(copy.get(64));
    CHECK(moved.words(0).data() == nullptr);
    moved.resize(0, 1);
    CHECK_FALSE(moved.any(1));
    moved = DynamicBitset{copy, 65};
    moved.set(64, false);
    CHECK(copy.get(64));
}

TEST_CASE("DynamicBitset growth preserves bits and zeroes new words", "[bitset]") {
    DynamicBitset bits{1};
    bits.set(0);
    const auto* storage{bits.words(1).data()};
    bits.resize(1, 64);
    CHECK(bits.words(64).data() == storage);
    bits.resize(64, 129);
    CHECK(bits.get(0));
    CHECK_FALSE(bits.get(64));
    CHECK_FALSE(bits.get(128));
    bits.set(128);
    bits.resize(129, 64);
    CHECK(DynamicBitset::allocatedBytes(64) == sizeof(std::uint64_t));
    bits.resize(64, 129);
    CHECK(bits.get(0));
    CHECK_FALSE(bits.get(128));
    bits.resize(129, 0);
    CHECK(DynamicBitset::allocatedBytes(0) == 0);
    CHECK_FALSE(bits.any(0));
    const DynamicBitset empty_copy{bits, 0};
    CHECK(empty_copy.words(0).empty());
}

TEST_CASE("DynamicBitset reuses geometric capacity and clears truncated words", "[bitset]") {
    DynamicBitset bits{129};
    bits.set(128);
    const auto* storage{bits.words(129).data()};
    CHECK(DynamicBitset::allocatedBytes(129) == 32);
    REQUIRE(bits.words(129).size() == 3);

    bits.resize(129, 193);
    CHECK(bits.words(193).data() == storage);
    CHECK(DynamicBitset::allocatedBytes(193) == 32);
    CHECK(bits.get(128));
    CHECK_FALSE(bits.get(192));
    bits.set(192);

    bits.resize(193, 129);
    CHECK(bits.words(129).data() == storage);
    bits.resize(129, 193);
    CHECK_FALSE(bits.get(192));

    DynamicBitset copy{bits, 193};
    copy.resize(193, 257);
    CHECK(DynamicBitset::allocatedBytes(257) == 64);
    CHECK(copy.get(128));
    CHECK_FALSE(copy.get(192));
    CHECK_FALSE(copy.get(256));
    CHECK(bits.get(128));
}
