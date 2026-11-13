#ifndef STABAAL_BACKENDS_PAULI_NOISE_BACKEND_HPP
#define STABAAL_BACKENDS_PAULI_NOISE_BACKEND_HPP

#include <cstddef>
#include <cstdint>
#include <optional>

#include "stabaal/backends/simulation_backend.hpp"
#include "stabaal/noise/pauli_noise.hpp"

namespace stabaal {

/**
 * Forward ideal operations and optionally apply a Pauli channel after each gate.
 * Noise is sampled independently for each qubit. Copies preserve the noise RNG state.
 * @tparam Backend Backend used to store and update the quantum state.
 */
template <SimulationBackend Backend>
class PauliNoiseBackend {
   public:
    /**
     * Create a backend with num_qubits qubits initialized to |0>.
     * @param num_qubits Number of qubits.
     */
    explicit PauliNoiseBackend(const QubitSize num_qubits) : backend_(num_qubits) {}

    /**
     * Create a backend with num_qubits qubits initialized to |0>.
     * Measurement and noise sampling use independent random engines.
     * @param num_qubits Number of qubits.
     * @param seed Measurement random seed, or std::nullopt for automatic seeding.
     * @param noise_seed Noise random seed, or std::nullopt for automatic seeding.
     */
    explicit PauliNoiseBackend(const QubitSize num_qubits, const std::optional<std::uint64_t> seed,
                               const std::optional<std::uint64_t> noise_seed = std::nullopt)
        requires std::constructible_from<Backend, QubitSize, std::optional<std::uint64_t>>
        : backend_(num_qubits, seed),
          rng_(noise_seed ? RandomEngine{*noise_seed} : RandomEngine{}) {}

    /**
     * Return the number of active qubits.
     * @return Number of active qubits.
     */
    [[nodiscard]] QubitSize NumQubits() const noexcept { return backend_.NumQubits(); }

    /**
     * Allocate a qubit initialized to |0>.
     * @return The stable ID of the allocated qubit.
     */
    QubitId Add() { return backend_.Add(); }

    /**
     * Trace out qubit q, preserving the remaining qubits' reduced state.
     * @param q Active qubit ID.
     */
    void Discard(const QubitId q) { backend_.Discard(q); }

    /**
     * Measure qubit q in the X basis and update the state to match the outcome.
     * @param q Active qubit ID.
     * @return The outcome: false for eigenvalue +1, true for eigenvalue -1.
     */
    bool MeasureX(const QubitId q) { return backend_.MeasureX(q); }

    /**
     * Measure qubit q in the Z basis and update the state to match the outcome.
     * @param q Active qubit ID.
     * @return The outcome: false for eigenvalue +1, true for eigenvalue -1.
     */
    bool MeasureZ(const QubitId q) { return backend_.MeasureZ(q); }

    /**
     * Measure qubit q in the X basis, then discard it.
     * @param q Active qubit ID.
     * @return The outcome: false for eigenvalue +1, true for eigenvalue -1.
     */
    bool MeasureAndDiscardX(const QubitId q) { return backend_.MeasureAndDiscardX(q); }

    /**
     * Measure qubit q in the Z basis, then discard it.
     * @param q Active qubit ID.
     * @return The outcome: false for eigenvalue +1, true for eigenvalue -1.
     */
    bool MeasureAndDiscardZ(const QubitId q) { return backend_.MeasureAndDiscardZ(q); }

    /**
     * Check whether measuring qubit q in the X basis is deterministic.
     * @param q Active qubit ID.
     * @return True if the outcome is certain, false otherwise.
     */
    [[nodiscard]] bool IsDeterministicMeasureX(const QubitId q) const
        noexcept(noexcept(backend_.IsDeterministicMeasureX(q))) {
        return backend_.IsDeterministicMeasureX(q);
    }

    /**
     * Force a nondeterministic X measurement of qubit q to the given outcome.
     * @param q Active qubit ID.
     * @param measurement Desired outcome: false for eigenvalue +1, true for eigenvalue -1.
     * @return The forced measurement outcome.
     * @pre The X measurement of q is nondeterministic.
     */
    bool ForceMeasureX(const QubitId q, const bool measurement) {
        return backend_.ForceMeasureX(q, measurement);
    }

    /**
     * Check whether measuring qubit q in the Z basis is deterministic.
     * @param q Active qubit ID.
     * @return True if the outcome is certain, false otherwise.
     */
    [[nodiscard]] bool IsDeterministicMeasureZ(const QubitId q) const
        noexcept(noexcept(backend_.IsDeterministicMeasureZ(q))) {
        return backend_.IsDeterministicMeasureZ(q);
    }

    /**
     * Force a nondeterministic Z measurement of qubit q to the given outcome.
     * @param q Active qubit ID.
     * @param measurement Desired outcome: false for eigenvalue +1, true for eigenvalue -1.
     * @return The forced measurement outcome.
     * @pre The Z measurement of q is nondeterministic.
     */
    bool ForceMeasureZ(const QubitId q, const bool measurement) {
        return backend_.ForceMeasureZ(q, measurement);
    }

    /**
     * Apply a Hadamard gate to qubit q.
     * @param q Active qubit ID.
     */
    void ApplyH(const QubitId q) noexcept(noexcept(backend_.ApplyH(q))) { backend_.ApplyH(q); }

    /**
     * Apply a Hadamard gate to qubit q.
     * Sample one Pauli error and apply it to q after the gate.
     * @param q Active qubit ID.
     * @param noise Pauli noise channel applied to qubit q.
     */
    void ApplyH(const QubitId q,
                const PauliNoise& noise) noexcept(noexcept(backend_.ApplyH(q)) &&
                                                  noexcept(noise.Apply(backend_, q, rng_))) {
        backend_.ApplyH(q);
        noise.Apply(backend_, q, rng_);
    }

    /**
     * Apply an S phase gate to qubit q.
     * @param q Active qubit ID.
     */
    void ApplyS(const QubitId q) noexcept(noexcept(backend_.ApplyS(q))) { backend_.ApplyS(q); }

    /**
     * Apply an S phase gate to qubit q.
     * Sample one Pauli error and apply it to q after the gate.
     * @param q Active qubit ID.
     * @param noise Pauli noise channel applied to qubit q.
     */
    void ApplyS(const QubitId q,
                const PauliNoise& noise) noexcept(noexcept(backend_.ApplyS(q)) &&
                                                  noexcept(noise.Apply(backend_, q, rng_))) {
        backend_.ApplyS(q);
        noise.Apply(backend_, q, rng_);
    }

    /**
     * Apply an S dagger gate to qubit q.
     * @param q Active qubit ID.
     */
    void ApplySdg(const QubitId q) noexcept(noexcept(backend_.ApplySdg(q))) {
        backend_.ApplySdg(q);
    }

    /**
     * Apply an S dagger gate to qubit q.
     * Sample one Pauli error and apply it to q after the gate.
     * @param q Active qubit ID.
     * @param noise Pauli noise channel applied to qubit q.
     */
    void ApplySdg(const QubitId q,
                  const PauliNoise& noise) noexcept(noexcept(backend_.ApplySdg(q)) &&
                                                    noexcept(noise.Apply(backend_, q, rng_))) {
        backend_.ApplySdg(q);
        noise.Apply(backend_, q, rng_);
    }

    /**
     * Apply a square-root X gate to qubit q.
     * @param q Active qubit ID.
     */
    void ApplySx(const QubitId q) noexcept(noexcept(backend_.ApplySx(q))) { backend_.ApplySx(q); }

    /**
     * Apply a square-root X gate to qubit q.
     * Sample one Pauli error and apply it to q after the gate.
     * @param q Active qubit ID.
     * @param noise Pauli noise channel applied to qubit q.
     */
    void ApplySx(const QubitId q,
                 const PauliNoise& noise) noexcept(noexcept(backend_.ApplySx(q)) &&
                                                   noexcept(noise.Apply(backend_, q, rng_))) {
        backend_.ApplySx(q);
        noise.Apply(backend_, q, rng_);
    }

    /**
     * Apply the identity gate to qubit q.
     * @param q Active qubit ID.
     */
    void ApplyI(const QubitId q) noexcept(noexcept(backend_.ApplyI(q))) { backend_.ApplyI(q); }

    /**
     * Apply the identity gate to qubit q.
     * Sample one Pauli error and apply it to q after the gate.
     * @param q Active qubit ID.
     * @param noise Pauli noise channel applied to qubit q.
     */
    void ApplyI(const QubitId q,
                const PauliNoise& noise) noexcept(noexcept(backend_.ApplyI(q)) &&
                                                  noexcept(noise.Apply(backend_, q, rng_))) {
        backend_.ApplyI(q);
        noise.Apply(backend_, q, rng_);
    }

    /**
     * Apply a Pauli X gate to qubit q.
     * @param q Active qubit ID.
     */
    void ApplyX(const QubitId q) noexcept(noexcept(backend_.ApplyX(q))) { backend_.ApplyX(q); }

    /**
     * Apply a Pauli X gate to qubit q.
     * Sample one Pauli error and apply it to q after the gate.
     * @param q Active qubit ID.
     * @param noise Pauli noise channel applied to qubit q.
     */
    void ApplyX(const QubitId q,
                const PauliNoise& noise) noexcept(noexcept(backend_.ApplyX(q)) &&
                                                  noexcept(noise.Apply(backend_, q, rng_))) {
        backend_.ApplyX(q);
        noise.Apply(backend_, q, rng_);
    }

    /**
     * Apply a Pauli Z gate to qubit q.
     * @param q Active qubit ID.
     */
    void ApplyZ(const QubitId q) noexcept(noexcept(backend_.ApplyZ(q))) { backend_.ApplyZ(q); }

    /**
     * Apply a Pauli Z gate to qubit q.
     * Sample one Pauli error and apply it to q after the gate.
     * @param q Active qubit ID.
     * @param noise Pauli noise channel applied to qubit q.
     */
    void ApplyZ(const QubitId q,
                const PauliNoise& noise) noexcept(noexcept(backend_.ApplyZ(q)) &&
                                                  noexcept(noise.Apply(backend_, q, rng_))) {
        backend_.ApplyZ(q);
        noise.Apply(backend_, q, rng_);
    }

    /**
     * Apply a Pauli Y gate to qubit q.
     * @param q Active qubit ID.
     */
    void ApplyY(const QubitId q) noexcept(noexcept(backend_.ApplyY(q))) { backend_.ApplyY(q); }

    /**
     * Apply a Pauli Y gate to qubit q.
     * Sample one Pauli error and apply it to q after the gate.
     * @param q Active qubit ID.
     * @param noise Pauli noise channel applied to qubit q.
     */
    void ApplyY(const QubitId q,
                const PauliNoise& noise) noexcept(noexcept(backend_.ApplyY(q)) &&
                                                  noexcept(noise.Apply(backend_, q, rng_))) {
        backend_.ApplyY(q);
        noise.Apply(backend_, q, rng_);
    }

    /**
     * Apply a controlled X gate.
     * @param control Active control qubit ID.
     * @param target Active target qubit ID.
     */
    void ApplyCx(const QubitId control,
                 const QubitId target) noexcept(noexcept(backend_.ApplyCx(control, target))) {
        backend_.ApplyCx(control, target);
    }

    /**
     * Apply a controlled X gate.
     * Sample and apply independent Pauli errors to control and target after the gate,
     * in that order.
     * @param control Active control qubit ID.
     * @param target Active target qubit ID.
     * @param control_noise Pauli noise channel applied to the control qubit.
     * @param target_noise Pauli noise channel applied to the target qubit.
     */
    void ApplyCx(const QubitId control, const QubitId target, const PauliNoise& control_noise,
                 const PauliNoise&
                     target_noise) noexcept(noexcept(backend_.ApplyCx(control, target)) &&
                                            noexcept(control_noise.Apply(backend_, control,
                                                                         rng_)) &&
                                            noexcept(target_noise.Apply(backend_, target, rng_))) {
        backend_.ApplyCx(control, target);
        control_noise.Apply(backend_, control, rng_);
        target_noise.Apply(backend_, target, rng_);
    }

    /**
     * Apply a controlled Z gate.
     * @param control Active control qubit ID.
     * @param target Active target qubit ID.
     */
    void ApplyCz(const QubitId control,
                 const QubitId target) noexcept(noexcept(backend_.ApplyCz(control, target))) {
        backend_.ApplyCz(control, target);
    }

    /**
     * Apply a controlled Z gate.
     * Sample and apply independent Pauli errors to control and target after the gate,
     * in that order.
     * @param control Active control qubit ID.
     * @param target Active target qubit ID.
     * @param control_noise Pauli noise channel applied to the control qubit.
     * @param target_noise Pauli noise channel applied to the target qubit.
     */
    void ApplyCz(const QubitId control, const QubitId target, const PauliNoise& control_noise,
                 const PauliNoise&
                     target_noise) noexcept(noexcept(backend_.ApplyCz(control, target)) &&
                                            noexcept(control_noise.Apply(backend_, control,
                                                                         rng_)) &&
                                            noexcept(target_noise.Apply(backend_, target, rng_))) {
        backend_.ApplyCz(control, target);
        control_noise.Apply(backend_, control, rng_);
        target_noise.Apply(backend_, target, rng_);
    }

    /**
     * Apply a controlled Y gate.
     * @param control Active control qubit ID.
     * @param target Active target qubit ID.
     */
    void ApplyCy(const QubitId control,
                 const QubitId target) noexcept(noexcept(backend_.ApplyCy(control, target))) {
        backend_.ApplyCy(control, target);
    }

    /**
     * Apply a controlled Y gate.
     * Sample and apply independent Pauli errors to control and target after the gate,
     * in that order.
     * @param control Active control qubit ID.
     * @param target Active target qubit ID.
     * @param control_noise Pauli noise channel applied to the control qubit.
     * @param target_noise Pauli noise channel applied to the target qubit.
     */
    void ApplyCy(const QubitId control, const QubitId target, const PauliNoise& control_noise,
                 const PauliNoise&
                     target_noise) noexcept(noexcept(backend_.ApplyCy(control, target)) &&
                                            noexcept(control_noise.Apply(backend_, control,
                                                                         rng_)) &&
                                            noexcept(target_noise.Apply(backend_, target, rng_))) {
        backend_.ApplyCy(control, target);
        control_noise.Apply(backend_, control, rng_);
        target_noise.Apply(backend_, target, rng_);
    }

    /**
     * Print the underlying backend's state representation to standard output.
     */
    void Print() const noexcept(noexcept(backend_.Print())) { backend_.Print(); }

    /**
     * Return persistent memory usage, including the wrapper and backend allocations.
     * @return Memory usage in bytes.
     */
    [[nodiscard]] std::size_t MemoryUsage() const noexcept {
        return sizeof(*this) - sizeof(Backend) + backend_.MemoryUsage();
    }

   private:
    Backend backend_;
    RandomEngine rng_;
};

}  // namespace stabaal

#endif
