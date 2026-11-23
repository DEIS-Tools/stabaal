#ifndef STABAAL_SIMULATOR_HPP
#define STABAAL_SIMULATOR_HPP

#include <cstddef>
#include <cstdint>
#include <limits>
#include <optional>
#include <stdexcept>

#include "stabaal/backends/inverse_tableau.hpp"
#include "stabaal/backends/simulation_backend.hpp"
#include "stabaal/noise/pauli_noise.hpp"
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
    explicit Simulator(const QubitSize num_qubits) : backend_(num_qubits) {}

    /**
     * Create a simulator with num_qubits qubits initialized to |0>^num_qubits.
     * Initial IDs are 0 through num_qubits - 1. IDs belong to this simulator;
     * copies preserve existing IDs and then evolve independently.
     * @param num_qubits Number of qubits.
     * @param seed Measurement random seed, or std::nullopt for automatic seeding.
     */
    explicit Simulator(const QubitSize num_qubits, const std::optional<std::uint64_t> seed)
        requires std::constructible_from<Backend, QubitSize, std::optional<std::uint64_t>>
        : backend_(num_qubits, seed) {}

    /**
     * Create a simulator with num_qubits qubits initialized to |0>^num_qubits.
     * Initial IDs are 0 through num_qubits - 1. IDs belong to this simulator;
     * copies preserve existing IDs and then evolve independently.
     * Measurement and noise sampling use independent random engines.
     * @param num_qubits Number of qubits.
     * @param seed Measurement random seed, or std::nullopt for automatic seeding.
     * @param noise_seed Noise random seed, or std::nullopt for automatic seeding.
     */
    explicit Simulator(const QubitSize num_qubits, const std::optional<std::uint64_t> seed,
                       const std::optional<std::uint64_t> noise_seed)
        requires std::constructible_from<Backend, QubitSize, std::optional<std::uint64_t>,
                                         std::optional<std::uint64_t>>
        : backend_(num_qubits, seed, noise_seed) {}

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
    void Discard(const QubitId q) { backend_.Discard(q); }

    /**
     * Measure qubit q in the X basis and update the state to match the outcome.
     * @param q Qubit id.
     * @return The outcome: false for eigenvalue +1, true for eigenvalue -1.
     */
    bool MeasureX(const QubitId q) { return backend_.MeasureX(q); }

    /**
     * Measure qubit q in the Z basis and update the state to match the outcome.
     * @param q Qubit id.
     * @return The outcome: false for eigenvalue +1, true for eigenvalue -1.
     */
    bool MeasureZ(const QubitId q) { return backend_.MeasureZ(q); }

    /**
     * Measure and discard qubit q in the X basis.
     * @param id Qubit id.
     * @return The outcome: false for eigenvalue +1, true for eigenvalue -1.
     */
    bool MeasureAndDiscardX(const QubitId id) { return backend_.MeasureAndDiscardX(id); }

    /**
     * Measure and discard qubit q in the Z basis.
     * @param id Qubit id.
     * @return The outcome: false for eigenvalue +1, true for eigenvalue -1.
     */
    bool MeasureAndDiscardZ(const QubitId id) { return backend_.MeasureAndDiscardZ(id); }

    /**
     * Check whether measuring qubit q in the X basis is deterministic.
     * @param q Qubit id.
     * @return True if the measurement is deterministic, false otherwise.
     */
    [[nodiscard]] bool IsDeterministicMeasureX(const QubitId q) const
        noexcept(noexcept(backend_.IsDeterministicMeasureX(q))) {
        return backend_.IsDeterministicMeasureX(q);
    }

    /**
     * Check whether measuring qubit q in the Z basis is deterministic.
     * @param q Qubit id.
     * @return True if the measurement is deterministic, false otherwise.
     */
    [[nodiscard]] bool IsDeterministicMeasureZ(const QubitId q) const
        noexcept(noexcept(backend_.IsDeterministicMeasureZ(q))) {
        return backend_.IsDeterministicMeasureZ(q);
    }

    /**
     * Force a non-deterministic X measurement of qubit q to the given outcome.
     * @param q Qubit id.
     * @param measurement Desired outcome: false for eigenvalue +1, true for eigenvalue -1.
     * @return The forced measurement outcome.
     */
    bool ForceMeasureX(const QubitId q, const bool measurement) {
        return backend_.ForceMeasureX(q, measurement);
    }

    /**
     * Force a non-deterministic Z measurement of qubit q to the given outcome.
     * @param q Qubit id.
     * @param measurement Desired outcome: false for eigenvalue +1, true for eigenvalue -1.
     * @return The forced measurement outcome.
     */
    bool ForceMeasureZ(const QubitId q, const bool measurement) {
        return backend_.ForceMeasureZ(q, measurement);
    }

    /**
     * Apply a Hadamard gate to qubit q.
     * @param q Stable qubit id.
     */
    void H(const QubitId q) noexcept(noexcept(backend_.ApplyH(q))) { backend_.ApplyH(q); }

    /**
     * Apply a Hadamard gate to qubit q.
     * Sample one Pauli error and apply it to q after the gate.
     * @param q Stable qubit id.
     * @param noise Pauli noise channel applied to qubit q.
     */
    void H(const QubitId q, const PauliNoise& noise) noexcept(requires {
        { backend_.ApplyH(q, noise) } noexcept;
    })
        requires requires(Backend& backend) { backend.ApplyH(q, noise); }
    {
        backend_.ApplyH(q, noise);
    }

    /**
     * Apply an S phase gate to qubit q.
     * @param q Stable qubit id.
     */
    void S(const QubitId q) noexcept(noexcept(backend_.ApplyS(q))) { backend_.ApplyS(q); }

    /**
     * Apply an S phase gate to qubit q.
     * Sample one Pauli error and apply it to q after the gate.
     * @param q Stable qubit id.
     * @param noise Pauli noise channel applied to qubit q.
     */
    void S(const QubitId q, const PauliNoise& noise) noexcept(requires {
        { backend_.ApplyS(q, noise) } noexcept;
    })
        requires requires(Backend& backend) { backend.ApplyS(q, noise); }
    {
        backend_.ApplyS(q, noise);
    }

    /**
     * Apply an S† gate to qubit q.
     * @param q Qubit id.
     */
    void Sdg(const QubitId q) noexcept(noexcept(backend_.ApplySdg(q))) { backend_.ApplySdg(q); }

    /**
     * Apply an S† gate to qubit q.
     * Sample one Pauli error and apply it to q after the gate.
     * @param q Qubit id.
     * @param noise Pauli noise channel applied to qubit q.
     */
    void Sdg(const QubitId q, const PauliNoise& noise) noexcept(requires {
        { backend_.ApplySdg(q, noise) } noexcept;
    })
        requires requires(Backend& backend) { backend.ApplySdg(q, noise); }
    {
        backend_.ApplySdg(q, noise);
    }

    /**
     * Apply an √X gate to qubit q.
     * @param q Qubit id.
     */
    void Sx(const QubitId q) noexcept(noexcept(backend_.ApplySx(q))) { backend_.ApplySx(q); }

    /**
     * Apply an √X gate to qubit q.
     * Sample one Pauli error and apply it to q after the gate.
     * @param q Qubit id.
     * @param noise Pauli noise channel applied to qubit q.
     */
    void Sx(const QubitId q, const PauliNoise& noise) noexcept(requires {
        { backend_.ApplySx(q, noise) } noexcept;
    })
        requires requires(Backend& backend) { backend.ApplySx(q, noise); }
    {
        backend_.ApplySx(q, noise);
    }

    /**
     * Apply the identity gate to qubit q.
     * @param q Qubit id.
     */
    void I(const QubitId q) noexcept(noexcept(backend_.ApplyI(q))) { backend_.ApplyI(q); }

    /**
     * Apply the identity gate to qubit q.
     * Sample one Pauli error and apply it to q after the gate.
     * @param q Qubit id.
     * @param noise Pauli noise channel applied to qubit q.
     */
    void I(const QubitId q, const PauliNoise& noise) noexcept(requires {
        { backend_.ApplyI(q, noise) } noexcept;
    })
        requires requires(Backend& backend) { backend.ApplyI(q, noise); }
    {
        backend_.ApplyI(q, noise);
    }

    /**
     * Apply a Pauli X gate to qubit q.
     * @param q Qubit id.
     */
    void X(const QubitId q) noexcept(noexcept(backend_.ApplyX(q))) { backend_.ApplyX(q); }

    /**
     * Apply a Pauli X gate to qubit q.
     * Sample one Pauli error and apply it to q after the gate.
     * @param q Qubit id.
     * @param noise Pauli noise channel applied to qubit q.
     */
    void X(const QubitId q, const PauliNoise& noise) noexcept(requires {
        { backend_.ApplyX(q, noise) } noexcept;
    })
        requires requires(Backend& backend) { backend.ApplyX(q, noise); }
    {
        backend_.ApplyX(q, noise);
    }

    /**
     * Apply a Pauli Z gate to qubit q.
     * @param q Qubit id.
     */
    void Z(const QubitId q) noexcept(noexcept(backend_.ApplyZ(q))) { backend_.ApplyZ(q); }

    /**
     * Apply a Pauli Z gate to qubit q.
     * Sample one Pauli error and apply it to q after the gate.
     * @param q Qubit id.
     * @param noise Pauli noise channel applied to qubit q.
     */
    void Z(const QubitId q, const PauliNoise& noise) noexcept(requires {
        { backend_.ApplyZ(q, noise) } noexcept;
    })
        requires requires(Backend& backend) { backend.ApplyZ(q, noise); }
    {
        backend_.ApplyZ(q, noise);
    }

    /**
     * Apply a Pauli Y gate to qubit q.
     * @param q Qubit id.
     */
    void Y(const QubitId q) noexcept(noexcept(backend_.ApplyY(q))) { backend_.ApplyY(q); }

    /**
     * Apply a Pauli Y gate to qubit q.
     * Sample one Pauli error and apply it to q after the gate.
     * @param q Qubit id.
     * @param noise Pauli noise channel applied to qubit q.
     */
    void Y(const QubitId q, const PauliNoise& noise) noexcept(requires {
        { backend_.ApplyY(q, noise) } noexcept;
    })
        requires requires(Backend& backend) { backend.ApplyY(q, noise); }
    {
        backend_.ApplyY(q, noise);
    }

    /**
     * Apply a controlled X gate from control to target.
     * @param control Control qubit id.
     * @param target Target qubit id.
     */
    void Cx(const QubitId control,
            const QubitId target) noexcept(noexcept(backend_.ApplyCx(control, target))) {
        backend_.ApplyCx(control, target);
    }

    /**
     * Apply a controlled X gate from control to target.
     * Sample and apply independent Pauli errors to control and target after the gate,
     * in that order.
     * @param control Control qubit id.
     * @param target Target qubit id.
     * @param control_noise Pauli noise channel applied to the control qubit.
     * @param target_noise Pauli noise channel applied to the target qubit.
     */
    void Cx(const QubitId control, const QubitId target, const PauliNoise& control_noise,
            const PauliNoise& target_noise) noexcept(requires {
        { backend_.ApplyCx(control, target, control_noise, target_noise) } noexcept;
    })
        requires requires(Backend& backend) {
            backend.ApplyCx(control, target, control_noise, target_noise);
        }
    {
        backend_.ApplyCx(control, target, control_noise, target_noise);
    }

    /**
     * Apply a controlled Z gate from control to target.
     * @param control Control qubit id.
     * @param target Target qubit id.
     */
    void Cz(const QubitId control,
            const QubitId target) noexcept(noexcept(backend_.ApplyCz(control, target))) {
        backend_.ApplyCz(control, target);
    }

    /**
     * Apply a controlled Z gate from control to target.
     * Sample and apply independent Pauli errors to control and target after the gate,
     * in that order.
     * @param control Control qubit id.
     * @param target Target qubit id.
     * @param control_noise Pauli noise channel applied to the control qubit.
     * @param target_noise Pauli noise channel applied to the target qubit.
     */
    void Cz(const QubitId control, const QubitId target, const PauliNoise& control_noise,
            const PauliNoise& target_noise) noexcept(requires {
        { backend_.ApplyCz(control, target, control_noise, target_noise) } noexcept;
    })
        requires requires(Backend& backend) {
            backend.ApplyCz(control, target, control_noise, target_noise);
        }
    {
        backend_.ApplyCz(control, target, control_noise, target_noise);
    }

    /**
     * Apply a controlled Y gate from control to target.
     * @param control Control qubit id.
     * @param target Target qubit id.
     */
    void Cy(const QubitId control,
            const QubitId target) noexcept(noexcept(backend_.ApplyCy(control, target))) {
        backend_.ApplyCy(control, target);
    }

    /**
     * Apply a controlled Y gate from control to target.
     * Sample and apply independent Pauli errors to control and target after the gate,
     * in that order.
     * @param control Control qubit id.
     * @param target Target qubit id.
     * @param control_noise Pauli noise channel applied to the control qubit.
     * @param target_noise Pauli noise channel applied to the target qubit.
     */
    void Cy(const QubitId control, const QubitId target, const PauliNoise& control_noise,
            const PauliNoise& target_noise) noexcept(requires {
        { backend_.ApplyCy(control, target, control_noise, target_noise) } noexcept;
    })
        requires requires(Backend& backend) {
            backend.ApplyCy(control, target, control_noise, target_noise);
        }
    {
        backend_.ApplyCy(control, target, control_noise, target_noise);
    }

    /**
     * Print the backend's state representation to standard output.
     */
    void Print() const noexcept(noexcept(backend_.Print())) { backend_.Print(); }

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
