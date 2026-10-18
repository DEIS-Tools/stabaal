#include "dynamic_bitset.hpp"

#include <algorithm>
#include <cassert>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <span>
#include <stdexcept>
#include <utility>

using namespace stabaal;

namespace {

constexpr std::size_t kBitsPerWord{std::numeric_limits<std::uint64_t>::digits};

}  // namespace

DynamicBitset::DynamicBitset() noexcept : num_bits_{0} {}

DynamicBitset::DynamicBitset(const std::size_t num_bits)
    : words_(num_bits / kBitsPerWord + (num_bits % kBitsPerWord != 0 ? 1 : 0)),
      num_bits_{num_bits} {}

DynamicBitset::DynamicBitset(const DynamicBitset& other) = default;
DynamicBitset& DynamicBitset::operator=(const DynamicBitset& other) = default;
DynamicBitset::~DynamicBitset() noexcept = default;

DynamicBitset::DynamicBitset(DynamicBitset&& other) noexcept
    : words_{std::move(other.words_)}, num_bits_{std::exchange(other.num_bits_, 0)} {}

DynamicBitset& DynamicBitset::operator=(DynamicBitset&& other) noexcept {
    if (this != &other) {
        words_ = std::move(other.words_);
        num_bits_ = std::exchange(other.num_bits_, 0);
    }

    return *this;
}

std::size_t DynamicBitset::Size() const noexcept { return num_bits_; }

bool DynamicBitset::Any() const noexcept {
    return std::any_of(words_.begin(), words_.end(),
                       [](const std::uint64_t word) noexcept { return word != 0; });
}

bool DynamicBitset::Get(const std::size_t idx) const {
    if (idx >= num_bits_) [[unlikely]] {
        throw std::out_of_range{"Bit index out of range"};
    }

    return (words_[idx / kBitsPerWord] & (std::uint64_t{1} << (idx % kBitsPerWord))) != 0;
}

void DynamicBitset::Set(const std::size_t idx, const bool value) {
    if (idx >= num_bits_) [[unlikely]] {
        throw std::out_of_range{"Bit index out of range"};
    }

    const auto mask{std::uint64_t{1} << (idx % kBitsPerWord)};
    auto& word{words_[idx / kBitsPerWord]};
    word = (word & ~mask) | (value ? mask : 0);
}

void DynamicBitset::PushBack(const bool value) {
    if (num_bits_ % kBitsPerWord == 0) words_.push_back(0);
    ++num_bits_;
    if (value) Set(num_bits_ - 1);
}

DynamicBitset& DynamicBitset::operator^=(const DynamicBitset& other) noexcept {
    assert(num_bits_ == other.num_bits_);
    for (std::size_t i{0}; i < words_.size(); ++i) {
        words_[i] ^= other.words_[i];
    }

    return *this;
}

std::span<const std::uint64_t> DynamicBitset::Words() const noexcept { return words_; }
