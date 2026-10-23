#include "stabaal/detail/bit_ops.hpp"

#include <cstddef>

#include "stabaal/types.hpp"

namespace stabaal::detail {

namespace {

/**
 * Offset from beginning of a word
 * @param index Bit index.
 * @return Bit position within its word.
 */
[[nodiscard]] std::size_t GetWordOffset(const std::size_t index) noexcept {
    return index % kBitsPerWord;
}

/**
 * Mask with a 1 at position index
 * @param index Bit index.
 * @return A word with only the selected bit set.
 */
[[nodiscard]] BitWord GetMask(const std::size_t index) noexcept {
    return BitWord{1} << GetWordOffset(index);
}

/**
 * Find the word containing the selected bit.
 * @param words Words containing the bit.
 * @param index Bit index across all words, not a word index.
 * @return A read-only reference to the word containing the bit.
 */
[[nodiscard]] const BitWord& GetWord(const ConstBitWordSpan words,
                                     const std::size_t index) noexcept {
    return words[index / kBitsPerWord];
}

/**
 * Find the word containing the selected bit.
 * @param words Words containing the bit.
 * @param index Bit index across all words, not a word index.
 * @return A writable reference to the word containing the bit.
 */
[[nodiscard]] BitWord& GetWord(const BitWordSpan words, const std::size_t index) noexcept {
    return words[index / kBitsPerWord];
}

}  // namespace

std::size_t WordCount(const std::size_t num_bits) noexcept {
    return num_bits / kBitsPerWord + (num_bits % kBitsPerWord != 0 ? 1 : 0);
}

bool GetBit(const ConstBitWordSpan words, const std::size_t index) noexcept {
    return (GetWord(words, index) & GetMask(index)) != 0;
}

void SetBit(const BitWordSpan words, const std::size_t index, const bool value) noexcept {
    const auto mask{GetMask(index)};
    auto& word{GetWord(words, index)};
    word = (word & ~mask) | (BitWord{value} << GetWordOffset(index));
}

void FlipBit(const BitWordSpan words, const std::size_t index, const bool flip) noexcept {
    GetWord(words, index) ^= BitWord{flip} << GetWordOffset(index);
}

}  // namespace stabaal::detail
