#include "dynamic_bitset.hpp"

#include <algorithm>
#include <cassert>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <span>

using namespace stabaal;

namespace {

constexpr std::size_t kBitsPerWord{std::numeric_limits<std::uint64_t>::digits};

}  // namespace

DynamicBitset::DynamicBitset() noexcept = default;

DynamicBitset::DynamicBitset(const std::size_t num_bits)
    : words_(num_bits / kBitsPerWord + (num_bits % kBitsPerWord != 0 ? 1 : 0)) {}

DynamicBitset::DynamicBitset(const DynamicBitset& other) = default;
DynamicBitset& DynamicBitset::operator=(const DynamicBitset& other) = default;
DynamicBitset::~DynamicBitset() noexcept = default;

DynamicBitset::DynamicBitset(DynamicBitset&& other) noexcept = default;
DynamicBitset& DynamicBitset::operator=(DynamicBitset&& other) noexcept = default;

std::size_t DynamicBitset::AllocatedBytes() const noexcept {
    return words_.capacity() * sizeof(std::uint64_t);
}

bool DynamicBitset::Any() const noexcept {
    return std::any_of(words_.begin(), words_.end(),
                       [](const std::uint64_t word) noexcept { return word != 0; });
}

bool DynamicBitset::Get(const std::size_t idx) const noexcept {
    assert(idx / kBitsPerWord < words_.size());
    return (words_[idx / kBitsPerWord] & (std::uint64_t{1} << (idx % kBitsPerWord))) != 0;
}

void DynamicBitset::Set(const std::size_t idx, const bool value) noexcept {
    assert(idx / kBitsPerWord < words_.size());
    const auto mask{std::uint64_t{1} << (idx % kBitsPerWord)};
    auto& word{words_[idx / kBitsPerWord]};
    word = (word & ~mask) | (value ? mask : 0);
}

void DynamicBitset::Resize(const std::size_t num_bits) {
    words_.resize(num_bits / kBitsPerWord + (num_bits % kBitsPerWord != 0 ? 1 : 0));
}

std::span<std::uint64_t> DynamicBitset::Words() noexcept { return words_; }

std::span<const std::uint64_t> DynamicBitset::Words() const noexcept { return words_; }
