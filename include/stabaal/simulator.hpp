#ifndef STABAAL_SIMULATOR_HPP
#define STABAAL_SIMULATOR_HPP

#include <concepts>
#include <cstddef>
#include <limits>
#include <stdexcept>

#include "detail/inverse_tableau.hpp"
#include "types.hpp"

namespace stabaal {

/**
 * Operations required of a Simulator backend.
 */
template <typename T>
concept SimulationBackend =
    std::constructible_from<T, QubitSize> && requires(T& backend, const T& const_backend, QubitId q,
                                                      QubitId c, QubitId t, bool measurement) {
        { const_backend.NumQubits() } noexcept -> std::same_as<QubitSize>;
        { backend.Add() } -> std::same_as<QubitId>;
        { backend.Discard(q) } -> std::same_as<void>;
        { backend.MeasureX(q) } -> std::same_as<bool>;
        { backend.MeasureZ(q) } -> std::same_as<bool>;
        { const_backend.IsDeterministicMeasureX(q) } -> std::same_as<bool>;
        { const_backend.IsDeterministicMeasureZ(q) } -> std::same_as<bool>;
        { backend.ForceMeasureX(q, measurement) } -> std::same_as<bool>;
        { backend.ForceMeasureZ(q, measurement) } -> std::same_as<bool>;
        { backend.ApplyH(q) } -> std::same_as<void>;
        { backend.ApplyS(q) } -> std::same_as<void>;
        { backend.ApplySdg(q) } -> std::same_as<void>;
        { backend.ApplySx(q) } -> std::same_as<void>;
        { backend.ApplyI(q) } -> std::same_as<void>;
        { backend.ApplyX(q) } -> std::same_as<void>;
        { backend.ApplyZ(q) } -> std::same_as<void>;
        { backend.ApplyY(q) } -> std::same_as<void>;
        { backend.ApplyCx(c, t) } -> std::same_as<void>;
        { backend.ApplyCz(c, t) } -> std::same_as<void>;
        { backend.ApplyCy(c, t) } -> std::same_as<void>;
        { const_backend.Print() } -> std::same_as<void>;
        { const_backend.MemoryUsage() } noexcept -> std::same_as<std::size_t>;
    };

template <SimulationBackend Backend = detail::InverseTableau>
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
     * Apply a controlled X gate with control c and target t.
     * @param c Control qubit id.
     * @param t Target qubit id.
     */
    void Cx(QubitId c, QubitId t) { backend_.ApplyCx(c, t); }

    /**
     * Apply a controlled Z gate with control c and target t.
     * @param c Control qubit id.
     * @param t Target qubit id.
     */
    void Cz(QubitId c, QubitId t) { backend_.ApplyCz(c, t); }

    /**
     * Apply a controlled Y gate with control c and target t.
     * @param c Control qubit id.
     * @param t Target qubit id.
     */
    void Cy(QubitId c, QubitId t) { backend_.ApplyCy(c, t); }

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
