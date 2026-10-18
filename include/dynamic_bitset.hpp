#ifndef STABAAL_DYNAMIC_BITSET_HPP
#define STABAAL_DYNAMIC_BITSET_HPP

#include <cstddef>
#include <cstdint>
#include <span>
#include <vector>

namespace stabaal {

class DynamicBitset {
   public:
    DynamicBitset() noexcept;
    explicit DynamicBitset(std::size_t num_bits);
    DynamicBitset(const DynamicBitset& other);
    DynamicBitset& operator=(const DynamicBitset& other);
    DynamicBitset(DynamicBitset&& other) noexcept;
    DynamicBitset& operator=(DynamicBitset&& other) noexcept;
    ~DynamicBitset() noexcept;

    std::size_t Size() const noexcept;
    std::size_t AllocatedBytes() const noexcept;
    bool Any() const noexcept;
    bool Get(std::size_t idx) const noexcept;
    void Set(std::size_t idx, bool value = true) noexcept;
    void PushBack(bool value = false);
    std::span<std::uint64_t> Words() noexcept;
    std::span<const std::uint64_t> Words() const noexcept;

   private:
    std::vector<std::uint64_t> words_;
    std::size_t num_bits_;
};

}  // namespace stabaal

#endif
