#ifndef STABAAL_TYPES_HPP
#define STABAAL_TYPES_HPP

#include <cstddef>
#include <cstdint>
#include <limits>
#include <span>

namespace stabaal {

using QubitIdx = std::uint32_t;
using BitWord = std::uint64_t;
using PhaseExponent = std::int64_t;
using BitWordSpan = std::span<BitWord>;
using ConstBitWordSpan = std::span<const BitWord>;

// Keep phase arithmetic wide enough for the full QubitIdx range
static_assert((std::numeric_limits<PhaseExponent>::max() - 3) / 4 >=
              std::numeric_limits<QubitIdx>::max());

constexpr std::size_t kBitsPerWord{std::numeric_limits<BitWord>::digits};

}  // namespace stabaal

#endif
