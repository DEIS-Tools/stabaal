#ifndef STABAAL_DYNAMIC_ARRAY_HPP
#define STABAAL_DYNAMIC_ARRAY_HPP

#include <concepts>
#include <cstddef>
#include <memory>

namespace stabaal {

// Stores only the array data to avoid overhead, and size must be tracked by callers.
template <std::default_initializable T>
class DynamicArray {
   public:
    DynamicArray() noexcept = default;

    explicit DynamicArray(const std::size_t num_elements)
        // NOLINTNEXTLINE(cppcoreguidelines-avoid-c-arrays)
        : data_{num_elements == 0 ? nullptr : std::make_unique<T[]>(num_elements)} {}

    DynamicArray(const DynamicArray&) = delete;
    DynamicArray& operator=(const DynamicArray&) = delete;
    DynamicArray(DynamicArray&&) noexcept = default;
    DynamicArray& operator=(DynamicArray&&) noexcept = default;
    ~DynamicArray() noexcept = default;

    T& operator[](const std::size_t index) noexcept { return data_[index]; }
    const T& operator[](const std::size_t index) const noexcept { return data_[index]; }

    T* Data() noexcept { return data_.get(); }
    const T* Data() const noexcept { return data_.get(); }

    explicit operator bool() const noexcept { return static_cast<bool>(data_); }

   private:
    std::unique_ptr<T[]> data_; // NOLINT(cppcoreguidelines-avoid-c-arrays)
};

}  // namespace stabaal

#endif
