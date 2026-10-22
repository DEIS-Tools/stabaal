#ifndef STABAAL_DETAIL_BIT_OPS_HPP
#define STABAAL_DETAIL_BIT_OPS_HPP

#include <cstddef>

#include "stabaal/types.hpp"

namespace stabaal::detail {

/**
 * Number of words required to store num_bits.
 * @param num_bits Number of bits to store.
 * @return Number of words needed to represent num_bits.
 */
[[nodiscard]] std::size_t WordCount(std::size_t num_bits) noexcept;

/**
 * Get the bit at index.
 * @param words Words containing the bit.
 * @param index Bit index.
 * @return The bit value at index.
 */
[[nodiscard]] bool GetBit(ConstBitWordSpan words, std::size_t index) noexcept;

/**
 * Set the bit at index to value.
 * @param words Words containing the bit.
 * @param index Bit index.
 * @param value Value to write.
 */
void SetBit(BitWordSpan words, std::size_t index, bool value) noexcept;

/**
 * Flip the bit at index if flip is true.
 * @param words Words containing the bit.
 * @param index Bit index.
 * @param flip Whether to flip the bit or leave it unchanged.
 */
void FlipBit(BitWordSpan words, std::size_t index, bool flip = true) noexcept;

}  // namespace stabaal::detail

#endif
