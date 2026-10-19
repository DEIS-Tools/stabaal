#ifndef STABAAL_TYPES_HPP
#define STABAAL_TYPES_HPP

#include <cstddef>
#include <cstdint>
#include <memory>

namespace stabaal {

using QubitIdx = std::uint32_t;

template <typename T>
using DynamicArray = std::unique_ptr<T[]>;  // NOLINT(cppcoreguidelines-avoid-c-arrays)

template <typename T>
DynamicArray<T> MakeDynamicArray(const std::size_t num_elements) {
    if (num_elements == 0) return {};
    return std::make_unique<T[]>(num_elements);  // NOLINT(cppcoreguidelines-avoid-c-arrays)
}

}  // namespace stabaal

#endif
