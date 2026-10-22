#ifndef STABAAL_DETAIL_BIT_WORD_HPP
#define STABAAL_DETAIL_BIT_WORD_HPP

#include <cstddef>

#include "stabaal/types.hpp"

namespace stabaal::detail {

[[nodiscard]] std::size_t WordCount(std::size_t num_bits) noexcept;
[[nodiscard]] bool GetBit(ConstBitWordSpan words, std::size_t idx) noexcept;
void SetBit(BitWordSpan words, std::size_t idx, bool value) noexcept;
void FlipBit(BitWordSpan words, std::size_t idx, bool flip = true) noexcept;

}  // namespace stabaal::detail

#endif
