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
    DynamicBitset(const DynamicBitset &other);
    DynamicBitset &operator=(const DynamicBitset &other);
    DynamicBitset(DynamicBitset &&other) noexcept;
    DynamicBitset &operator=(DynamicBitset &&other) noexcept;
    ~DynamicBitset() noexcept;

    std::size_t size() const noexcept;
    std::size_t allocatedBytes() const noexcept;
    bool any() const noexcept;
    bool get(std::size_t idx) const noexcept;
    void set(std::size_t idx, bool value = true) noexcept;
    void pushBack(bool value = false);
    std::span<std::uint64_t> words() noexcept;
    std::span<const std::uint64_t> words() const noexcept;

   private:
    std::vector<std::uint64_t> words_;
    std::size_t num_bits_;
};

}  // namespace stabaal

#endif
