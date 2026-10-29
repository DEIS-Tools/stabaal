#ifndef STABAAL_DETAIL_RANDOM_HPP
#define STABAAL_DETAIL_RANDOM_HPP

#include <bit>
#include <concepts>
#include <cstdint>
#include <limits>
#include <random>

namespace stabaal::detail {
// From https://martin.ankerl.com/2018/12/08/fast-random-bool/
class Sfc64 {
   public:
    // NOLINTNEXTLINE(readability-identifier-naming): standard engine API.
    using result_type = std::uint64_t;

    // NOLINTNEXTLINE(readability-identifier-naming): required by uniform_random_bit_generator.
    [[nodiscard]] static constexpr result_type min() noexcept { return 0; }
    // NOLINTNEXTLINE(readability-identifier-naming): required by uniform_random_bit_generator.
    [[nodiscard]] static constexpr result_type max() noexcept {
        return std::numeric_limits<result_type>::max();
    }

    Sfc64() : Sfc64(std::random_device{}()) {}
    explicit Sfc64(const std::uint64_t seed) noexcept
        : m_a_{seed}, m_b_{seed}, m_c_{seed}, m_counter_{seed} {
        for (auto i = 0; i < 12; ++i) {
            static_cast<void>(operator()());
        }
    }

    [[nodiscard]] uint64_t operator()() noexcept {
        auto const tmp = m_a_ + m_b_ + m_counter_++;
        m_a_ = m_b_ ^ (m_b_ >> right_shift_);
        m_b_ = m_c_ + (m_c_ << left_shift_);
        m_c_ = std::rotl(m_c_, rotation_) + tmp;
        return tmp;
    }

   private:
    static constexpr int rotation_{24};
    static constexpr unsigned right_shift_{11};
    static constexpr unsigned left_shift_{3};

    uint64_t m_a_;
    uint64_t m_b_;
    uint64_t m_c_;
    uint64_t m_counter_;
};

template <std::unsigned_integral U>
class Random {
   public:
    template <std::uniform_random_bit_generator Rng>
    [[nodiscard]] bool operator()(Rng& rng) {
        if (rand_vals_ == 1) [[unlikely]] {
            rand_vals_ = std::uniform_int_distribution<U>{}(rng) | sentinel_mask;
        }

        const bool bit = rand_vals_ & U{1};
        rand_vals_ >>= 1U;
        return bit;
    }

   private:
    static constexpr U sentinel_mask = U{1} << (sizeof(U) * 8 - 1);
    U rand_vals_{1};
};
}  // namespace stabaal::detail

#endif
