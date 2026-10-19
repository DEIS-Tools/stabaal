#ifndef STABAAL_TYPES_HPP
#define STABAAL_TYPES_HPP

#include <cstdint>
#include <memory>

namespace stabaal {

using QubitIdx = std::uint32_t;

template <typename T>
using DynamicArray = std::unique_ptr<T[]>;  // NOLINT(cppcoreguidelines-avoid-c-arrays)

}  // namespace stabaal

#endif
