#ifndef STABAAL_DETAIL_PAULI_NOISE_HPP
#define STABAAL_DETAIL_PAULI_NOISE_HPP

#include <cstdint>

#include "stabaal/detail/random.hpp"

namespace stabaal::detail {
enum class PauliError : std::uint8_t { kI, kX, kY, kZ };

class PauliNoise {
   public:
    /**
     * PauliNoise channel
     * @param wi Weight of no error occuring.
     * @param wx Weight of a bit flip occuring.
     * @param wy Weight of bit flip and phase flip occuring.
     * @param wz Weight of phase flip occuring.
     */
    explicit PauliNoise(double wi, double wx, double wy, double wz);

    /**
     * Sample a PauliError
     * @return PauliError sampled
     */
    [[nodiscard]] PauliError operator()(Sfc64& rng) const noexcept;

   private:
    double i_end_;
    double x_end_;
    double y_end_;
};
}  // namespace stabaal::detail

#endif