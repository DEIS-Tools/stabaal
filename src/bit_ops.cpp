#include "stabaal/detail/bit_ops.hpp"

#include <cstddef>

#include "stabaal/types.hpp"

namespace stabaal::detail {

namespace {

[[nodiscard]] BitWord GetMask(const std::size_t idx) noexcept {
    return BitWord{1} << (idx % kBitsPerWord);
}

[[nodiscard]] const BitWord& GetWord(const ConstBitWordSpan words, const std::size_t idx) noexcept {
    return words[idx / kBitsPerWord];
}

[[nodiscard]] BitWord& GetWord(const BitWordSpan words, const std::size_t idx) noexcept {
    return words[idx / kBitsPerWord];
}

}  // namespace

std::size_t WordCount(const std::size_t num_bits) noexcept {
    return num_bits / kBitsPerWord + (num_bits % kBitsPerWord != 0 ? 1 : 0);
}

bool GetBit(const ConstBitWordSpan words, const std::size_t idx) noexcept {
    return (GetWord(words, idx) & GetMask(idx)) != 0;
}

void SetBit(const BitWordSpan words, const std::size_t idx, const bool value) noexcept {
    const auto mask{GetMask(idx)};
    auto& word{GetWord(words, idx)};
    if (value) {
        word |= mask;
        return;
    }
    
    word &= ~mask;
}

void FlipBit(const BitWordSpan words, const std::size_t idx, const bool flip) noexcept {
    GetWord(words, idx) ^= flip ? GetMask(idx) : BitWord{0};
}

}  // namespace stabaal::detail
