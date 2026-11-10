#ifndef STABAAL_DETAIL_PAULI_NOISE_BACKEND_HPP
#define STABAAL_DETAIL_PAULI_NOISE_BACKEND_HPP

#include "simulation_backend.hpp"

namespace stabaal::detail {
template <SimulationBackend Backend>
class PauliNoiseBackend {
   private:
    Backend backend_;
};
}  // namespace stabaal::detail

#endif