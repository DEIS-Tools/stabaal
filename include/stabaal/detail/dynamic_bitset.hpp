#ifndef STABAAL_DETAIL_DYNAMIC_BITSET_HPP
#define STABAAL_DETAIL_DYNAMIC_BITSET_HPP

#include <cstddef>
#include <cstdint>
#include <span>

#include "dynamic_array.hpp"

namespace stabaal::detail {

class DynamicBitset {
   public:
    DynamicBitset() noexcept;
    explicit DynamicBitset(std::size_t num_bits);
    DynamicBitset(const DynamicBitset& other, std::size_t num_bits);
    DynamicBitset(const DynamicBitset& other) = delete;
    DynamicBitset& operator=(const DynamicBitset& other) = delete;
    DynamicBitset(DynamicBitset&& other) noexcept;
    DynamicBitset& operator=(DynamicBitset&& other) noexcept;
    ~DynamicBitset() noexcept;

    static std::size_t allocatedBytes(std::size_t num_bits) noexcept;
    bool any(std::size_t num_bits) const noexcept;
    bool get(std::size_t idx) const noexcept;
    void set(std::size_t idx, bool value = true) noexcept;
    void resize(std::size_t old_num_bits, std::size_t new_num_bits);
    void copyFrom(const DynamicBitset& other, std::size_t num_bits) noexcept;
    std::span<std::uint64_t> words(std::size_t num_bits) noexcept;
    std::span<const std::uint64_t> words(std::size_t num_bits) const noexcept;

   private:
    DynamicArray<std::uint64_t> words_;
};

}  // namespace stabaal::detail

#endif
