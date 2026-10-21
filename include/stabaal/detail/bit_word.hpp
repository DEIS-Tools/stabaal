#ifndef STABAAL_DETAIL_BIT_WORD_HPP
#define STABAAL_DETAIL_BIT_WORD_HPP

#include <cstddef>
#include <cstdint>
#include <limits>
#include <span>

namespace stabaal::detail {

inline constexpr std::size_t kBitsPerWord{std::numeric_limits<std::uint64_t>::digits};

[[nodiscard]] constexpr std::size_t WordCount(const std::size_t num_bits) noexcept {
    return num_bits / kBitsPerWord + (num_bits % kBitsPerWord != 0 ? 1 : 0);
}

[[nodiscard]] constexpr bool GetBit(const std::span<const std::uint64_t> words,
                                    const std::size_t idx) noexcept {
    return ((words[idx / kBitsPerWord] >> (idx % kBitsPerWord)) & 1U) != 0;
}

constexpr void SetBit(const std::span<std::uint64_t> words, const std::size_t idx,
                      const bool value) noexcept {
    const auto mask{std::uint64_t{1} << (idx % kBitsPerWord)};
    auto& word{words[idx / kBitsPerWord]};
    word = (word & ~mask) | (value ? mask : 0);
}

constexpr void FlipBit(const std::span<std::uint64_t> words, const std::size_t idx,
                       const bool flip = true) noexcept {
    words[idx / kBitsPerWord] ^= static_cast<std::uint64_t>(flip) << (idx % kBitsPerWord);
}

}  // namespace stabaal::detail

#endif
