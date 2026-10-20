#ifndef STABAAL_DYNAMIC_BITSET_HPP
#define STABAAL_DYNAMIC_BITSET_HPP

#include <cstddef>
#include <cstdint>
#include <span>

#include "dynamic_array.hpp"

namespace stabaal {

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

    static std::size_t AllocatedBytes(std::size_t num_bits) noexcept;
    bool Any(std::size_t num_bits) const noexcept;
    bool Get(std::size_t idx) const noexcept;
    void Set(std::size_t idx, bool value = true) noexcept;
    void Resize(std::size_t old_num_bits, std::size_t new_num_bits);
    std::span<std::uint64_t> Words(std::size_t num_bits) noexcept;
    std::span<const std::uint64_t> Words(std::size_t num_bits) const noexcept;

   private:
    DynamicArray<std::uint64_t> words_;
};

}  // namespace stabaal

#endif
