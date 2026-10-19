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
    const auto n{GENERATE(0, 1, 63, 64, 65, 127, 128, 129, 193, 257, 513, 1024, 1025)};
    const DynamicBitset bits{static_cast<std::size_t>(n)};
    CHECK(bits.Words(n).size() == (n + 63) / 64);
    const auto used_bytes{bits.Words(n).size() * sizeof(std::uint64_t)};
    CHECK(DynamicBitset::AllocatedBytes(n) >= used_bytes);
    if (n == 0) {
        CHECK(DynamicBitset::AllocatedBytes(n) == 0);
    } else {
        CHECK(DynamicBitset::AllocatedBytes(n) < 2 * used_bytes);
    }
    CHECK_FALSE(bits.Any(n));
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
            bits.Resize(expected.size(), expected.size() + 1);
            bits.Set(expected.size(), value);
            expected.push_back(value);
        } else {
            const auto idx{rng() % expected.size()};
            bits.Set(idx, value);
            expected[idx] = value;
        }
        REQUIRE(bits.Words(expected.size()).size() == (expected.size() + 63) / 64);
        bool any{false};
        for (std::size_t i{0}; i < expected.size(); ++i) {
            REQUIRE(bits.Get(i) == expected[i]);
            any = any || expected[i];
        }
        CHECK(bits.Any(expected.size()) == any);
        const auto used{expected.size() % 64};
        if (used != 0) CHECK((bits.Words(expected.size()).back() >> used) == 0);
    }
}

TEST_CASE("DynamicBitset copies own their storage and moved-from objects can grow", "[bitset]") {
    DynamicBitset original{65};
    original.Set(64);
    DynamicBitset copy{original, 65};
    copy.Set(64, false);
    CHECK(original.Get(64));
    CHECK_FALSE(copy.Any(65));

    DynamicBitset moved{std::move(original)};
    REQUIRE(moved.Words(65).size() == 2);
    CHECK(moved.Get(64));
    CHECK(original.Words(0).data() == nullptr);
    CHECK(DynamicBitset::AllocatedBytes(0) == 0);
    original.Resize(0, 1);
    original.Set(0);
    REQUIRE(original.Words(1).size() == 1);
    CHECK(original.Get(0));

    copy = std::move(moved);
    CHECK(copy.Get(64));
    CHECK(moved.Words(0).data() == nullptr);
    moved.Resize(0, 1);
    CHECK_FALSE(moved.Any(1));
    moved = DynamicBitset{copy, 65};
    moved.Set(64, false);
    CHECK(copy.Get(64));
}

TEST_CASE("DynamicBitset growth preserves bits and zeroes new words", "[bitset]") {
    DynamicBitset bits{1};
    bits.Set(0);
    const auto* storage{bits.Words(1).data()};
    bits.Resize(1, 64);
    CHECK(bits.Words(64).data() == storage);
    bits.Resize(64, 129);
    CHECK(bits.Get(0));
    CHECK_FALSE(bits.Get(64));
    CHECK_FALSE(bits.Get(128));
    bits.Set(128);
    bits.Resize(129, 64);
    CHECK(DynamicBitset::AllocatedBytes(64) == sizeof(std::uint64_t));
    bits.Resize(64, 129);
    CHECK(bits.Get(0));
    CHECK_FALSE(bits.Get(128));
    bits.Resize(129, 0);
    CHECK(DynamicBitset::AllocatedBytes(0) == 0);
    CHECK_FALSE(bits.Any(0));
    const DynamicBitset empty_copy{bits, 0};
    CHECK(empty_copy.Words(0).empty());
}

TEST_CASE("DynamicBitset reuses geometric capacity and clears truncated words", "[bitset]") {
    DynamicBitset bits{129};
    bits.Set(128);
    const auto* storage{bits.Words(129).data()};
    CHECK(DynamicBitset::AllocatedBytes(129) == 32);
    REQUIRE(bits.Words(129).size() == 3);

    bits.Resize(129, 193);
    CHECK(bits.Words(193).data() == storage);
    CHECK(DynamicBitset::AllocatedBytes(193) == 32);
    CHECK(bits.Get(128));
    CHECK_FALSE(bits.Get(192));
    bits.Set(192);

    bits.Resize(193, 129);
    CHECK(bits.Words(129).data() == storage);
    bits.Resize(129, 193);
    CHECK_FALSE(bits.Get(192));

    DynamicBitset copy{bits, 193};
    copy.Resize(193, 257);
    CHECK(DynamicBitset::AllocatedBytes(257) == 64);
    CHECK(copy.Get(128));
    CHECK_FALSE(copy.Get(192));
    CHECK_FALSE(copy.Get(256));
    CHECK(bits.Get(128));
}
