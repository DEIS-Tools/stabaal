#ifndef STABAAL_NOISE_NOISE_CHANNEL_HPP
#define STABAAL_NOISE_NOISE_CHANNEL_HPP

#include <concepts>

#include "stabaal/backends/simulation_backend.hpp"
#include "stabaal/random.hpp"
#include "stabaal/types.hpp"

namespace stabaal {
/** A channel that samples and applies noise through a backend. */
template <typename T, typename Backend>
concept NoiseChannel = SimulationBackend<Backend> &&
                       requires(const T& noise, Backend& backend, QubitId q, RandomEngine& rng) {
                           { noise.Apply(backend, q, rng) } -> std::same_as<void>;
                       };

}  // namespace stabaal

#endif
