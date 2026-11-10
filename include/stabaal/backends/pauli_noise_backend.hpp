#ifndef STABAAL_BACKENDS_PAULI_NOISE_BACKEND_HPP
#define STABAAL_BACKENDS_PAULI_NOISE_BACKEND_HPP

#include <cstddef>

#include "stabaal/backends/simulation_backend.hpp"

namespace stabaal {

/**
 * Forward normal simulator operations to the underlying backend.
 * @tparam Backend Backend used to store and update the quantum state.
 */
template <SimulationBackend Backend>
class PauliNoiseBackend {
   public:
    /**
     * Create a backend with num_qubits qubits initialized to |0>.
     * @param num_qubits Number of qubits.
     */
    explicit PauliNoiseBackend(QubitSize num_qubits) : backend_(num_qubits) {}

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
    void Discard(QubitId q) { backend_.Discard(q); }

    /**
     * Measure qubit q in the X basis and update the state to match the outcome.
     * @param q Active qubit ID.
     * @return The outcome: false for eigenvalue +1, true for eigenvalue -1.
     */
    bool MeasureX(QubitId q) { return backend_.MeasureX(q); }

    /**
     * Measure qubit q in the Z basis and update the state to match the outcome.
     * @param q Active qubit ID.
     * @return The outcome: false for eigenvalue +1, true for eigenvalue -1.
     */
    bool MeasureZ(QubitId q) { return backend_.MeasureZ(q); }

    /**
     * Measure qubit q in the X basis, then discard it.
     * @param q Active qubit ID.
     * @return The outcome: false for eigenvalue +1, true for eigenvalue -1.
     */
    bool MeasureAndDiscardX(QubitId q) { return backend_.MeasureAndDiscardX(q); }

    /**
     * Measure qubit q in the Z basis, then discard it.
     * @param q Active qubit ID.
     * @return The outcome: false for eigenvalue +1, true for eigenvalue -1.
     */
    bool MeasureAndDiscardZ(QubitId q) { return backend_.MeasureAndDiscardZ(q); }

    /**
     * Check whether measuring qubit q in the X basis is deterministic.
     * @param q Active qubit ID.
     * @return True if the outcome is certain, false otherwise.
     */
    [[nodiscard]] bool IsDeterministicMeasureX(QubitId q) const {
        return backend_.IsDeterministicMeasureX(q);
    }

    /**
     * Force a nondeterministic X measurement of qubit q to the given outcome.
     * @param q Active qubit ID.
     * @param measurement Desired outcome: false for eigenvalue +1, true for eigenvalue -1.
     * @return The forced measurement outcome.
     * @pre The X measurement of q is nondeterministic.
     */
    bool ForceMeasureX(QubitId q, bool measurement) {
        return backend_.ForceMeasureX(q, measurement);
    }

    /**
     * Check whether measuring qubit q in the Z basis is deterministic.
     * @param q Active qubit ID.
     * @return True if the outcome is certain, false otherwise.
     */
    [[nodiscard]] bool IsDeterministicMeasureZ(QubitId q) const {
        return backend_.IsDeterministicMeasureZ(q);
    }

    /**
     * Force a nondeterministic Z measurement of qubit q to the given outcome.
     * @param q Active qubit ID.
     * @param measurement Desired outcome: false for eigenvalue +1, true for eigenvalue -1.
     * @return The forced measurement outcome.
     * @pre The Z measurement of q is nondeterministic.
     */
    bool ForceMeasureZ(QubitId q, bool measurement) {
        return backend_.ForceMeasureZ(q, measurement);
    }

    /**
     * Apply a Hadamard gate to qubit q.
     * @param q Active qubit ID.
     */
    void ApplyH(QubitId q) { backend_.ApplyH(q); }

    /**
     * Apply an S phase gate to qubit q.
     * @param q Active qubit ID.
     */
    void ApplyS(QubitId q) { backend_.ApplyS(q); }

    /**
     * Apply an S dagger gate to qubit q.
     * @param q Active qubit ID.
     */
    void ApplySdg(QubitId q) { backend_.ApplySdg(q); }

    /**
     * Apply a square-root X gate to qubit q.
     * @param q Active qubit ID.
     */
    void ApplySx(QubitId q) { backend_.ApplySx(q); }

    /**
     * Apply the identity gate to qubit q.
     * @param q Active qubit ID.
     */
    void ApplyI(QubitId q) { backend_.ApplyI(q); }

    /**
     * Apply a Pauli X gate to qubit q.
     * @param q Active qubit ID.
     */
    void ApplyX(QubitId q) { backend_.ApplyX(q); }

    /**
     * Apply a Pauli Z gate to qubit q.
     * @param q Active qubit ID.
     */
    void ApplyZ(QubitId q) { backend_.ApplyZ(q); }

    /**
     * Apply a Pauli Y gate to qubit q.
     * @param q Active qubit ID.
     */
    void ApplyY(QubitId q) { backend_.ApplyY(q); }

    /**
     * Apply a controlled X gate.
     * @param control Active control qubit ID.
     * @param target Active target qubit ID.
     */
    void ApplyCx(QubitId control, QubitId target) { backend_.ApplyCx(control, target); }

    /**
     * Apply a controlled Z gate.
     * @param control Active control qubit ID.
     * @param target Active target qubit ID.
     */
    void ApplyCz(QubitId control, QubitId target) { backend_.ApplyCz(control, target); }

    /**
     * Apply a controlled Y gate.
     * @param control Active control qubit ID.
     * @param target Active target qubit ID.
     */
    void ApplyCy(QubitId control, QubitId target) { backend_.ApplyCy(control, target); }

    /**
     * Print the underlying backend's state representation to standard output.
     */
    void Print() const { backend_.Print(); }

    /**
     * Return persistent memory usage, including the wrapper and backend allocations.
     * @return Memory usage in bytes.
     */
    [[nodiscard]] std::size_t MemoryUsage() const noexcept {
        return sizeof(*this) - sizeof(Backend) + backend_.MemoryUsage();
    }

   private:
    Backend backend_;
};

}  // namespace stabaal

#endif
