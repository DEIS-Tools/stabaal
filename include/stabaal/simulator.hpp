#ifndef STABAAL_SIMULATOR_HPP
#define STABAAL_SIMULATOR_HPP

#include <cstddef>
#include <limits>
#include <stdexcept>

#include "stabaal/backends/inverse_tableau.hpp"
#include "stabaal/backends/simulation_backend.hpp"
#include "types.hpp"

namespace stabaal {

template <SimulationBackend Backend = InverseTableau>
class Simulator {
   public:
    /**
     * Create a simulator with num_qubits qubits initialized to |0>^num_qubits.
     * Initial IDs are 0 through num_qubits - 1. IDs belong to this simulator;
     * copies preserve existing IDs and then evolve independently.
     * @param num_qubits Number of qubits.
     */
    explicit Simulator(QubitSize num_qubits) : backend_(num_qubits) {}

    /**
     * Return the number of active qubits.
     * @return Number of active qubits.
     */
    [[nodiscard]] QubitSize NumQubits() const noexcept { return backend_.NumQubits(); }

    /**
     * Allocate a qubit initialized to |0>, growing storage geometrically if needed.
     * @return A stable id, distinct from any previously discarded id even when storage is reused.
     * @throws std::length_error If the qubit count would exceed the maximum QubitSize.
     */
    QubitId Add() {
        if (backend_.NumQubits() == std::numeric_limits<QubitSize>::max()) [[unlikely]] {
            throw std::length_error{"Too many qubits"};
        }

        return backend_.Add();
    }

    /**
     * Trace out qubit q exactly, preserving the remaining qubits' reduced state.
     * @param q Qubit id.
     */
    void Discard(QubitId q) { backend_.Discard(q); }

    /**
     * Measure qubit q in the X basis and update the state to match the outcome.
     * @param q Qubit id.
     * @return The outcome: false for eigenvalue +1, true for eigenvalue -1.
     */
    bool MeasureX(QubitId q) { return backend_.MeasureX(q); }

    /**
     * Measure qubit q in the Z basis and update the state to match the outcome.
     * @param q Qubit id.
     * @return The outcome: false for eigenvalue +1, true for eigenvalue -1.
     */
    bool MeasureZ(QubitId q) { return backend_.MeasureZ(q); }

    /**
     * Measure and discard qubit q in the X basis.
     * @param id Qubit id.
     * @return The outcome: false for eigenvalue +1, true for eigenvalue -1.
     */
    bool MeasureAndDiscardX(QubitId id) { return backend_.MeasureAndDiscardX(id); }

    /**
     * Measure and discard qubit q in the Z basis.
     * @param id Qubit id.
     * @return The outcome: false for eigenvalue +1, true for eigenvalue -1.
     */
    bool MeasureAndDiscardZ(QubitId id) { return backend_.MeasureAndDiscardZ(id); }

    /**
     * Check whether measuring qubit q in the X basis is deterministic.
     * @param q Qubit id.
     * @return True if the measurement is deterministic, false otherwise.
     */
    [[nodiscard]] bool IsDeterministicMeasureX(QubitId q) const {
        return backend_.IsDeterministicMeasureX(q);
    }

    /**
     * Check whether measuring qubit q in the Z basis is deterministic.
     * @param q Qubit id.
     * @return True if the measurement is deterministic, false otherwise.
     */
    [[nodiscard]] bool IsDeterministicMeasureZ(QubitId q) const {
        return backend_.IsDeterministicMeasureZ(q);
    }

    /**
     * Force a non-deterministic X measurement of qubit q to the given outcome.
     * @param q Qubit id.
     * @param measurement Desired outcome: false for eigenvalue +1, true for eigenvalue -1.
     * @return The forced measurement outcome.
     */
    bool ForceMeasureX(QubitId q, bool measurement) {
        return backend_.ForceMeasureX(q, measurement);
    }

    /**
     * Force a non-deterministic Z measurement of qubit q to the given outcome.
     * @param q Qubit id.
     * @param measurement Desired outcome: false for eigenvalue +1, true for eigenvalue -1.
     * @return The forced measurement outcome.
     */
    bool ForceMeasureZ(QubitId q, bool measurement) {
        return backend_.ForceMeasureZ(q, measurement);
    }

    /**
     * Apply a Hadamard gate to qubit q.
     * @param q Stable qubit id.
     */
    void H(QubitId q) { backend_.ApplyH(q); }

    /**
     * Apply an S phase gate to qubit q.
     * @param q Stable qubit id.
     */
    void S(QubitId q) { backend_.ApplyS(q); }

    /**
     * Apply an S† gate to qubit q.
     * @param q Qubit id.
     */
    void Sdg(QubitId q) { backend_.ApplySdg(q); }

    /**
     * Apply an √X gate to qubit q.
     * @param q Qubit id.
     */
    void Sx(QubitId q) { backend_.ApplySx(q); }

    /**
     * Apply the identity gate to qubit q.
     * @param q Qubit id.
     */
    void I(QubitId q) { backend_.ApplyI(q); }

    /**
     * Apply a Pauli X gate to qubit q.
     * @param q Qubit id.
     */
    void X(QubitId q) { backend_.ApplyX(q); }

    /**
     * Apply a Pauli Z gate to qubit q.
     * @param q Qubit id.
     */
    void Z(QubitId q) { backend_.ApplyZ(q); }

    /**
     * Apply a Pauli Y gate to qubit q.
     * @param q Qubit id.
     */
    void Y(QubitId q) { backend_.ApplyY(q); }

    /**
     * Apply a controlled X gate from control to target.
     * @param control Control qubit id.
     * @param target Target qubit id.
     */
    void Cx(QubitId control, QubitId target) { backend_.ApplyCx(control, target); }

    /**
     * Apply a controlled Z gate from control to target.
     * @param control Control qubit id.
     * @param target Target qubit id.
     */
    void Cz(QubitId control, QubitId target) { backend_.ApplyCz(control, target); }

    /**
     * Apply a controlled Y gate from control to target.
     * @param control Control qubit id.
     * @param target Target qubit id.
     */
    void Cy(QubitId control, QubitId target) { backend_.ApplyCy(control, target); }

    /**
     * Print the backend's state representation to standard output.
     */
    void Print() const { backend_.Print(); }

    /**
     * Return backend memory usage in bytes.
     * @return Memory usage in bytes.
     */
    [[nodiscard]] std::size_t MemoryUsage() const noexcept { return backend_.MemoryUsage(); }

   private:
    Backend backend_;
};
}  // namespace stabaal

#endif
