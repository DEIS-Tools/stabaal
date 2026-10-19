#include "dynamic_bitset.hpp"

#include <algorithm>
#include <bit>
#include <cassert>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <span>
#include <utility>

using namespace stabaal;

namespace {

constexpr std::size_t kBitsPerWord{std::numeric_limits<std::uint64_t>::digits};

constexpr std::size_t WordCount(const std::size_t num_bits) noexcept {
    return num_bits / kBitsPerWord + (num_bits % kBitsPerWord != 0 ? 1 : 0);
}

constexpr std::size_t CapacityWords(const std::size_t num_bits) noexcept {
    return num_bits == 0 ? 0 : std::bit_ceil(WordCount(num_bits));
}

}  // namespace

DynamicBitset::DynamicBitset() noexcept = default;

DynamicBitset::DynamicBitset(const std::size_t num_bits) {
    const auto num_words{CapacityWords(num_bits)};
    if (num_words != 0) words_ = std::make_unique<std::uint64_t[]>(num_words);
}

DynamicBitset::DynamicBitset(const DynamicBitset& other, const std::size_t num_bits)
    : DynamicBitset{num_bits} {
    if (num_bits != 0) std::copy_n(other.words_.get(), WordCount(num_bits), words_.get());
}

DynamicBitset::~DynamicBitset() noexcept = default;
DynamicBitset::DynamicBitset(DynamicBitset&& other) noexcept = default;
DynamicBitset& DynamicBitset::operator=(DynamicBitset&& other) noexcept = default;

std::size_t DynamicBitset::AllocatedBytes(const std::size_t num_bits) const noexcept {
    return CapacityWords(num_bits) * sizeof(std::uint64_t);
}

bool DynamicBitset::Any(const std::size_t num_bits) const noexcept {
    const auto words{Words(num_bits)};
    return std::any_of(words.begin(), words.end(),
                       [](const std::uint64_t word) noexcept { return word != 0; });
}

bool DynamicBitset::Get(const std::size_t idx) const noexcept {
    assert(words_);
    const auto word_idx = idx / kBitsPerWord;
    const auto bit_offset = idx % kBitsPerWord;
    const auto mask = uint64_t{1} << bit_offset;
    return (words_[word_idx] & mask) != 0;
}

void DynamicBitset::Set(const std::size_t idx, const bool value) noexcept {
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

void DynamicBitset::Resize(const std::size_t old_num_bits, const std::size_t new_num_bits) {
    const auto old_num_words{WordCount(old_num_bits)};
    const auto new_num_words{WordCount(new_num_bits)};
    if (old_num_words == new_num_words) return;
    if (CapacityWords(old_num_bits) == CapacityWords(new_num_bits)) {
        if (new_num_words < old_num_words) {
            std::fill_n(words_.get() + new_num_words, old_num_words - new_num_words, 0);
        }
        
        return;
    }

    DynamicBitset resized{new_num_bits};
    const auto retained{std::min(old_num_words, new_num_words)};
    if (retained != 0) std::copy_n(words_.get(), retained, resized.words_.get());
    *this = std::move(resized);
}

std::span<std::uint64_t> DynamicBitset::Words(const std::size_t num_bits) noexcept {
    return {words_.get(), WordCount(num_bits)};
}

std::span<const std::uint64_t> DynamicBitset::Words(const std::size_t num_bits) const noexcept {
    return {words_.get(), WordCount(num_bits)};
}
