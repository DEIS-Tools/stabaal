#ifndef STABAAL_DETAIL_INVERSE_TABLEAU_HPP
#define STABAAL_DETAIL_INVERSE_TABLEAU_HPP

#include <cstddef>
#include <cstdint>
#include <memory>
#include <optional>

#include "pauli_string_view.hpp"
#include "stabaal/types.hpp"

namespace stabaal::detail {

/**
 * Simulate Clifford gates and Pauli measurements using an inverse stabilizer tableau.
 */
class InverseTableau {
   public:
    enum class Basis : std::uint8_t {
        kComputational,
        kHadamard,
    };

    /**
     * Create a tableau with num_qubits qubits initialized to |0>^num_qubits.
     * @param num_qubits Number of qubits.
     */
    explicit InverseTableau(QubitSize num_qubits);

    /**
     * Copy the state into independent storage.
     * @param other Source tableau.
     */
    InverseTableau(const InverseTableau& other);

    /**
     * Replace this state with an independent copy of other.
     * @param other Source tableau.
     * @return A reference to this tableau.
     */
    InverseTableau& operator=(const InverseTableau& other);

    /**
     * Take ownership of other's tableau storage.
     * @param other Tableau whose storage is transferred.
     */
    InverseTableau(InverseTableau&& other) noexcept = default;

    /**
     * Replace this tableau by taking ownership of other's storage.
     * @param other Tableau whose storage is transferred.
     * @return A reference to this tableau.
     */
    InverseTableau& operator=(InverseTableau&& other) noexcept = default;

    /**
     * Release the tableau storage.
     */
    ~InverseTableau() noexcept = default;

    /**
     * Return the number of active qubits.
     * @return Number of active qubits.
     */
    [[nodiscard]] QubitSize NumQubits() const noexcept;

    /**
     * Append a qubit initialized to |0>, growing storage if needed geometrically.
     */
    void Add();

    /**
     * Discards qubit q from the tableau.
     * @param q Qubit index.
     */
    void Discard(QubitSize q);

    /**
     * Measure qubit q in the X basis and update the state to match the outcome.
     * @param q Qubit index.
     * @return The outcome: false for eigenvalue +1, true for eigenvalue -1.
     */
    bool MeasureX(QubitSize q);

    /**
     * Measure qubit q in the Z basis and update the state to match the outcome.
     * @param q Qubit index.
     * @return The outcome: false for eigenvalue +1, true for eigenvalue -1.
     */
    bool MeasureZ(QubitSize q);

    /**
     * Checks wether measuring qubit q in the X basis will be a deterministic measurement
     * @param q Qubit index.
     * @return The outcome: true if the measurement will be deterministic, false otherwise
     */
    bool IsDeterministicMeasureX(QubitSize q);

    /**
     * Checks wether measuring qubit q in the Z basis will be a deterministic measurement
     * @param q Qubit index.
     * @return The outcome: true if the measurement will be deterministic, false otherwise
     */
    bool IsDeterministicMeasureZ(QubitSize q);

    /**
     * Forces a non-deterministic measurement in the X basis to have a certain outcome
     * @param q Qubit index.
     * @param measurement the desired outcome of the measurement.
     * @return The outcome of the measurement
     */
    bool ForceMeasureX(QubitSize q, bool measurement);

    /**
     * Forces a non-deterministic measurement in the Z basis to have a certain outcome
     * @param q Qubit index.
     * @param measurement the desired outcome of the measurement.
     * @return The outcome of the measurement
     */
    bool ForceMeasureZ(QubitSize q, bool measurement);

    /**
     * Apply a Hadamard gate to qubit q.
     * @param q Qubit index.
     */
    void ApplyH(QubitSize q) noexcept;

    /**
     * Apply an S phase gate to qubit q.
     * @param q Qubit index.
     */
    void ApplyS(QubitSize q) noexcept;

    /**
     * Apply an S† gate to qubit q.
     * @param q Qubit index.
     */
    void ApplySdg(QubitSize q) noexcept;

    /**
     * Apply the identity gate to qubit q.
     * @param q Qubit index.
     */
    void ApplyI(QubitSize q) noexcept;

    /**
     * Apply a Pauli X gate to qubit q.
     * @param q Qubit index.
     */
    void ApplyX(QubitSize q) noexcept;

    /**
     * Apply a Pauli Z gate to qubit q.
     * @param q Qubit index.
     */
    void ApplyZ(QubitSize q) noexcept;

    /**
     * Apply a Pauli Y gate to qubit q.
     * @param q Qubit index.
     */
    void ApplyY(QubitSize q) noexcept;

    /**
     * Apply a controlled X gate with control c and target t.
     * @param c Control qubit index.
     * @param t Target qubit index.
     */
    void ApplyCx(QubitSize c, QubitSize t) noexcept;

    /**
     * Apply a controlled Z gate with control c and target t.
     * @param c Control qubit index.
     * @param t Target qubit index.
     */
    void ApplyCz(QubitSize c, QubitSize t) noexcept;

    /**
     * Apply a controlled Y gate with control c and target t.
     * @param c Control qubit index.
     * @param t Target qubit index.
     */
    void ApplyCy(QubitSize c, QubitSize t) noexcept;

    /**
     * Print the active generators' X components, Z components, and signs to standard output.
     */
    void Print() const;

    /**
     * Return the bytes used by this object and its allocated tableau, including spare capacity.
     * Excludes allocator overhead and temporary measurement storage.
     * @return Memory usage in bytes.
     */
    [[nodiscard]] std::size_t MemoryUsage() const noexcept;

   private:
    std::unique_ptr<BitWord[]> words_;  // NOLINT(cppcoreguidelines-avoid-c-arrays)

    /**
     * Number of qubits in the tableau.
     */
    QubitSize num_qubits_;

    /**
     * Number of qubits the tableau can currently hold.
     */
    QubitSize capacity_;

    /**
     * Return the total word count for generator components and flags at the given capacity.
     * @param capacity Number of qubits the tableau can hold currently.
     * @return Total number of words required for components and flags.
     */
    [[nodiscard]] static std::size_t NumWords(QubitSize capacity) noexcept;

    /**
     * Return a writable view of packed bits in a row.
     * @param row Row index.
     * @return A view of the selected row, including spare capacity.
     */
    [[nodiscard]] BitWordSpan RowWords(std::size_t row) noexcept;

    /**
     * Return a read-only view of packed bits in a row.
     * @param row Row index.
     * @return A view of the selected row, including spare capacity.
     */
    [[nodiscard]] ConstBitWordSpan RowWords(std::size_t row) const noexcept;

    /**
     * Return a writable view of packed flags.
     * @param index Storage group: 0 or 1 for signs, 2 for swap flags.
     * @return A view of the selected packed flags.
     */
    [[nodiscard]] BitWordSpan Flags(std::size_t index) noexcept;

    /**
     * Return a read-only view of packed flags.
     * @param index Storage group: 0 or 1 for signs, 2 for swap flags.
     * @return A view of the selected packed flags.
     */
    [[nodiscard]] ConstBitWordSpan Flags(std::size_t index) const noexcept;

    /**
     * Return storage group 0 or 1 for the requested generator and basis, accounting for basis
     * swaps.
     * @param generator_index Generator number.
     * @param basis Generator basis.
     * @return Storage group index, either 0 or 1.
     */
    [[nodiscard]] std::size_t PhysicalBasis(QubitSize generator_index, Basis basis) const noexcept;

    /**
     * Return a writable view of the selected generator, accounting for basis swaps.
     * @param generator_index Generator number.
     * @param basis Basis to get generator in.
     * @return A writable view of the generator's X and Z components.
     */
    [[nodiscard]] PauliStringView Generator(QubitSize generator_index, Basis basis) noexcept;

    /**
     * Return a writable view of the selected generator, not accounting for basis swaps.
     * @param generator_index Generator number.
     * @param basis Basis to get generator in.
     * @return A writable view of the stored generator's components.
     */
    [[nodiscard]] PauliStringView PhysicalGenerator(QubitSize generator_index,
                                                    std::size_t basis) noexcept;

    /**
     * Return the selected generator's sign: false for +, true for -.
     * @param generator_index Generator number.
     * @param basis Generator basis.
     * @return The sign: false for positive, true for negative.
     */
    [[nodiscard]] bool Sign(QubitSize generator_index, Basis basis) const noexcept;

    /**
     * Toggle the selected generator's sign in the logical basis.
     * @param generator_index Generator number.
     * @param basis Generator basis.
     */
    void FlipSign(QubitSize generator_index, Basis basis) noexcept;

    /**
     * Copy generator components and flags from other into already allocated storage.
     * The destination capacity must accommodate other's active qubits.
     * Does not update the active qubit count or allocate storage.
     * @param other Source tableau.
     */
    void CopyRows(const InverseTableau& other) noexcept;

    /**
     * Measure qubit q in the selected basis, collapsing the state only if the outcome is random.
     * @param q Qubit index.
     * @param basis Generator basis.
     * @return The outcome: false for eigenvalue +1, true for eigenvalue -1.
     */
    bool Measure(QubitSize q, Basis basis);

    /**
     * Checks wether measuring in the selected basis will be a deterministic measurement
     * @param q Qubit index.
     * @param basis the basis to measure in.
     * @return the outcome: true if the measurement will be deterministic, false otherwise
     */
    [[nodiscard]] bool IsDeterministicMeasure(QubitSize q, Basis basis);

    /**
     * Forces a non-deterministic measurement in the selected basis to have a deterministic outcome
     * @param q Qubit index.
     * @param basis the basis to measure in.
     * @return the ourcome of the measurement
     */
    bool ForceMeasure(QubitSize q, bool measurement, Basis basis);

    /**
     * Apply measurement cleanup to every active generator and update their signs.
     * @param pivot Control position for the cleanup Cx gates.
     * @param targets Packed bits marking target positions, excluding the pivot.
     * @param apply_s_dagger Whether to apply S dagger to the pivot after the Cx gates.
     */
    void CleanUp(QubitSize pivot, ConstBitWordSpan targets, bool apply_s_dagger) noexcept;

    /**
     * Update all generators and signs for a random measurement after cleanup.
     * @param pivot Position selected during cleanup.
     * @param p Sampled measurement result.
     */
    void Collapse(QubitSize pivot, bool p) noexcept;

    /**
     * Clean up and collapse the tableau for a measurement with a random outcome.
     * The generator for qubit q in basis b must contain an X component.
     * @param q Qubit index.
     * @param b Measurement basis.
     * @param measurement_opt Force measurement to a specific outcome if not std::nullopt.
     * @return The outcome: false for eigenvalue +1, true for eigenvalue -1.
     */
    bool RandomMeasurement(QubitSize q, Basis b,
                           std::optional<bool> measurement_opt = std::nullopt);
};

}  // namespace stabaal::detail

#endif
