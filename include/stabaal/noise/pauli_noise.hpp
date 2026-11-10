#ifndef STABAAL_NOISE_PAULI_NOISE_HPP
#define STABAAL_NOISE_PAULI_NOISE_HPP

#include <cstdint>

#include "stabaal/backends/simulation_backend.hpp"
#include "stabaal/random.hpp"

namespace stabaal {
enum class PauliError : std::uint8_t { kI, kX, kY, kZ };

class PauliNoise {
   public:
    /**
     * Create a Pauli noise channel from relative weights, normalized to probabilities.
     * @param wi Weight of no error (I).
     * @param wx Weight of a bit flip (X).
     * @param wy Weight of a combined bit and phase flip (Y).
     * @param wz Weight of a phase flip (Z).
     * @throws std::invalid_argument If a weight is nonfinite or negative, or all weights are zero.
     */
    explicit PauliNoise(double wi, double wx, double wy, double wz);

    /**
     * Sample one Pauli error
     * @param rng Random engine to advance when sampling.
     * @return The sampled error: I, X, Y, or Z.
     */
    [[nodiscard]] PauliError operator()(RandomEngine& rng) const noexcept;

    /**
     * Sample one Pauli error and apply it to qubit q through backend.
     * An identity sample leaves the backend unchanged.
     * @tparam Backend Simulator backend used to apply the sampled error.
     * @param backend Backend whose state receives the sampled error.
     * @param q Active qubit ID in backend.
     * @param rng Random engine to advance when sampling.
     */
    template <SimulationBackend Backend>
    void Apply(Backend& backend, QubitId q, RandomEngine& rng) const {
        switch ((*this)(rng)) {
            case PauliError::kI:
                break;
            case PauliError::kX:
                backend.ApplyX(q);
                break;
            case PauliError::kY:
                backend.ApplyY(q);
                break;
            case PauliError::kZ:
                backend.ApplyZ(q);
                break;
        }
    }

   private:
    double i_end_;
    double x_end_;
    double y_end_;
};
}  // namespace stabaal

#endif
