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
    std::constructible_from<T, QubitSize> &&
    requires(T& backend, const T& const_backend, QubitSize q, QubitSize c, QubitSize t) {
        { const_backend.NumQubits() } noexcept -> std::same_as<QubitSize>;
        { backend.Add() } -> std::same_as<void>;
        { backend.MeasureX(q) } -> std::same_as<bool>;
        { backend.MeasureZ(q) } -> std::same_as<bool>;
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
     * @param num_qubits Number of qubits.
     */
    explicit Simulator(QubitSize num_qubits) : backend_(num_qubits) {}

    /**
     * Append a qubit initialized to |0>, growing storage if needed.
     * @throws std::length_error If the qubit count would exceed the maximum QubitSize.
     */
    void Add() {
        if (backend_.NumQubits() == std::numeric_limits<QubitSize>::max()) [[unlikely]] {
            throw std::length_error{"Too many qubits"};
        }

        backend_.Add();
    }

    /**
     * Measure qubit q in the X basis and update the state to match the outcome.
     * @param q Qubit index.
     * @return The outcome: false for eigenvalue +1, true for eigenvalue -1.
     */
    bool MeasureX(QubitSize q) {
        CheckQubit(q);
        return backend_.MeasureX(q);
    }

    /**
     * Measure qubit q in the Z basis and update the state to match the outcome.
     * @param q Qubit index.
     * @return The outcome: false for eigenvalue +1, true for eigenvalue -1.
     */
    bool MeasureZ(QubitSize q) {
        CheckQubit(q);
        return backend_.MeasureZ(q);
    }

    /**
     * Apply a Hadamard gate to qubit q.
     * @param q Zero-based qubit index.
     */
    void H(QubitSize q) {
        CheckQubit(q);
        backend_.ApplyH(q);
    }

    /**
     * Apply an S phase gate to qubit q.
     * @param q Qubit index.
     */
    void S(QubitSize q) {
        CheckQubit(q);
        backend_.ApplyS(q);
    }

    /**
     * Apply an S† gate to qubit q.
     * @param q Qubit index.
     */
    void Sdg(QubitSize q) {
        CheckQubit(q);
        backend_.ApplySdg(q);
    }

    /**
     * Apply an SX gate to qubit q.
     * @param q Qubit index.
     */
    void Sx(QubitSize q) {
        CheckQubit(q);
        backend_.ApplySx(q);
    }

    /**
     * Apply the identity gate to qubit q.
     * @param q Qubit index.
     */
    void I(QubitSize q) {
        CheckQubit(q);
        backend_.ApplyI(q);
    }

    /**
     * Apply a Pauli X gate to qubit q.
     * @param q Qubit index.
     */
    void X(QubitSize q) {
        CheckQubit(q);
        backend_.ApplyX(q);
    }

    /**
     * Apply a Pauli Z gate to qubit q.
     * @param q Qubit index.
     */
    void Z(QubitSize q) {
        CheckQubit(q);
        backend_.ApplyZ(q);
    }

    /**
     * Apply a Pauli Y gate to qubit q.
     * @param q Qubit index.
     */
    void Y(QubitSize q) {
        CheckQubit(q);
        backend_.ApplyY(q);
    }

    /**
     * Apply a controlled X gate with control c and target t.
     * @param c Control qubit index.
     * @param t Target qubit index.
     */
    void Cx(QubitSize c, QubitSize t) {
        CheckQubit(c);
        CheckQubit(t);
        backend_.ApplyCx(c, t);
    }

    /**
     * Apply a controlled Z gate with control c and target t.
     * @param c Control qubit index.
     * @param t Target qubit index.
     */
    void Cz(QubitSize c, QubitSize t) {
        CheckQubit(c);
        CheckQubit(t);
        backend_.ApplyCz(c, t);
    }

    /**
     * Apply a controlled Y gate with control c and target t.
     * @param c Control qubit index.
     * @param t Target qubit index.
     */
    void Cy(QubitSize c, QubitSize t) {
        CheckQubit(c);
        CheckQubit(t);
        backend_.ApplyCy(c, t);
    }

    /**
     * Print the backend's state representation to standard output.
     */
    void Print() const { backend_.Print(); }

    /**
     * Return the memory usage in bytes reported by the backend.
     * @return Memory usage in bytes.
     */
    [[nodiscard]] std::size_t MemoryUsage() const noexcept { return backend_.MemoryUsage(); }

   private:
    /**
     * Check that qubit q is active.
     * @param q Qubit index.
     * @throws std::out_of_range If q is outside the active qubit range.
     */
    void CheckQubit(QubitSize q) const {
        if (q >= backend_.NumQubits()) [[unlikely]] {
            throw std::out_of_range{"Qubit index out of range"};
        }
    }

    Backend backend_;
};
}  // namespace stabaal

#endif
