#ifndef STABAAL_TYPES_HPP
#define STABAAL_TYPES_HPP

#include <cstddef>
#include <cstdint>
#include <limits>
#include <span>

namespace stabaal {

using QubitIdx = std::uint32_t;
using BitWord = std::uint64_t;
using BitWordSpan = std::span<BitWord>;
using ConstBitWordSpan = std::span<const BitWord>;

constexpr std::size_t kBitsPerWord{std::numeric_limits<BitWord>::digits};

}  // namespace stabaal

#endif
