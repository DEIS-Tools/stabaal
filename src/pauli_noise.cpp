#include "stabaal/detail/pauli_noise.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <random>
#include <stdexcept>

#include "stabaal/detail/random.hpp"

namespace stabaal::detail {
PauliNoise::PauliNoise(double wi, double wx, double wy, double wz) {
    for (const double weight : {wi, wx, wy, wz}) {
        if (!std::isfinite(weight) || weight < 0.0) {
            throw std::invalid_argument("Weights must be finite and nonnegative");
        }
    }

    const double largest = std::max({wi, wx, wy, wz});
    if (largest == 0.0) {
        throw std::invalid_argument("At least one weight must be positive");
    }

    wi /= largest;
    wx /= largest;
    wy /= largest;
    wz /= largest;

    const double total = wi + wx + wy + wz;
    i_end_ = wi / total;
    x_end_ = (wi + wx) / total;
    y_end_ = (wi + wx + wy) / total;
}

PauliError PauliNoise::operator()(Sfc64& rng) const noexcept {
    const double val = std::generate_canonical<double, std::numeric_limits<double>::digits>(rng);
    if (val < i_end_) {
        return PauliError::kI;
    }

    if (val < x_end_) {
        return PauliError::kX;
    }

    if (val < y_end_) {
        return PauliError::kY;
    }

    return PauliError::kZ;
}
}  // namespace stabaal::detail
