#include "stabaal/detail/dynamic_bitset.hpp"

#include <algorithm>
#include <bit>
#include <cassert>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <span>
#include <utility>

#include "stabaal/detail/dynamic_array.hpp"

namespace stabaal::detail {

namespace {

constexpr std::size_t kBitsPerWord{std::numeric_limits<std::uint64_t>::digits};

constexpr std::size_t wordCount(const std::size_t num_bits) noexcept {
    return num_bits / kBitsPerWord + (num_bits % kBitsPerWord != 0 ? 1 : 0);
}

constexpr std::size_t capacityWords(const std::size_t num_bits) noexcept {
    return num_bits == 0 ? 0 : std::bit_ceil(wordCount(num_bits));
}

}  // namespace

DynamicBitset::DynamicBitset() noexcept = default;

DynamicBitset::DynamicBitset(const std::size_t num_bits) : words_{capacityWords(num_bits)} {}

DynamicBitset::DynamicBitset(const DynamicBitset& other, const std::size_t num_bits)
    : DynamicBitset{num_bits} {
    if (num_bits != 0) std::copy_n(other.words_.data(), wordCount(num_bits), words_.data());
}

DynamicBitset::~DynamicBitset() noexcept = default;
DynamicBitset::DynamicBitset(DynamicBitset&& other) noexcept = default;
DynamicBitset& DynamicBitset::operator=(DynamicBitset&& other) noexcept = default;

std::size_t DynamicBitset::allocatedBytes(const std::size_t num_bits) noexcept {
    return capacityWords(num_bits) * sizeof(std::uint64_t);
}

bool DynamicBitset::any(const std::size_t num_bits) const noexcept {
    const auto used{words(num_bits)};
    return std::any_of(used.begin(), used.end(),
                       [](const std::uint64_t word) noexcept { return word != 0; });
}

bool DynamicBitset::get(const std::size_t idx) const noexcept {
    assert(words_);
    const auto word_idx = idx / kBitsPerWord;
    const auto bit_offset = idx % kBitsPerWord;
    const auto mask = uint64_t{1} << bit_offset;
    return (words_[word_idx] & mask) != 0;
}

void DynamicBitset::set(const std::size_t idx, const bool value) noexcept {
    assert(words_);
    const auto word_idx = idx / kBitsPerWord;
    const auto bit_offset = idx % kBitsPerWord;
    const auto mask{std::uint64_t{1} << bit_offset};
    auto& word{words_[word_idx]};
    if (value) {
        word |= mask;
        return;
    }

    word &= ~mask;
}

void DynamicBitset::resize(const std::size_t old_num_bits, const std::size_t new_num_bits) {
    const auto old_num_words{wordCount(old_num_bits)};
    const auto new_num_words{wordCount(new_num_bits)};
    if (old_num_words == new_num_words) return;
    if (capacityWords(old_num_bits) == capacityWords(new_num_bits)) {
        if (new_num_words < old_num_words) {
            const auto truncated{words(old_num_bits).subspan(new_num_words)};
            std::fill(truncated.begin(), truncated.end(), 0);
        }

        return;
    }

    DynamicBitset resized{new_num_bits};
    const auto retained{std::min(old_num_words, new_num_words)};
    if (retained != 0) std::copy_n(words_.data(), retained, resized.words_.data());
    *this = std::move(resized);
}

void DynamicBitset::copyFrom(const DynamicBitset& other, const std::size_t num_bits) noexcept {
    if (num_bits != 0) std::copy_n(other.words_.data(), wordCount(num_bits), words_.data());
}

std::span<std::uint64_t> DynamicBitset::words(const std::size_t num_bits) noexcept {
    return {words_.data(), wordCount(num_bits)};
}

std::span<const std::uint64_t> DynamicBitset::words(const std::size_t num_bits) const noexcept {
    return {words_.data(), wordCount(num_bits)};
}

}  // namespace stabaal::detail
