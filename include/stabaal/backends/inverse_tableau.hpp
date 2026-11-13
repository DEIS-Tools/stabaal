#ifndef STABAAL_BACKENDS_INVERSE_TABLEAU_HPP
#define STABAAL_BACKENDS_INVERSE_TABLEAU_HPP

#include <cstddef>
#include <cstdint>
#include <limits>
#include <memory>
#include <optional>
#include <span>

#include "stabaal/detail/pauli_string_view.hpp"
#include "stabaal/random.hpp"
#include "stabaal/types.hpp"

namespace stabaal::detail {

/**
 * Simulate pure and mixed stabilizer states using an inverse tableau and hidden purification slots.
 */
class InverseTableau {
   public:
    enum class Basis : std::uint8_t {
        kComputational,
        kHadamard,
    };
    struct GateInformation {
        std::string name;
        std::string description;
        u_int8_t arity;
    };

    /**
     * Create a tableau with num_qubits qubits initialized to |0>^num_qubits.
     * @param num_qubits Number of qubits.
     */
    explicit InverseTableau(QubitSize num_qubits);

    /**
     * Create a tableau with num_qubits qubits initialized to |0>^num_qubits.
     * @param num_qubits Number of qubits.
     * @param seed Measurement random seed, or std::nullopt for automatic seeding.
     */
    explicit InverseTableau(QubitSize num_qubits, std::optional<std::uint64_t> seed);

    /**
     * Copy the quantum and measurement RNG state into independent storage.
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
     * Allocate a qubit initialized to |0>, reusing a slot or growing storage geometrically.
     * @return The stable id of the allocated qubit.
     */
    QubitId Add();

    /**
     * Discard a qubit.
     * Storage contracts when half of the qubits are unused.
     * @param id Qubit id.
     * @throws std::out_of_range If id is not active.
     */
    void Discard(QubitId id);

    /**
     * Measure qubit q in the X basis and update the state to match the outcome.
     * @param id Qubit id.
     * @return The outcome: false for eigenvalue +1, true for eigenvalue -1.
     */
    bool MeasureX(QubitId id);

    /**
     * Measure qubit q in the Z basis and update the state to match the outcome.
     * @param id Qubit id.
     * @return The outcome: false for eigenvalue +1, true for eigenvalue -1.
     */
    bool MeasureZ(QubitId id);

    /**
     * Measure and discard qubit q in the X basis.
     * @param id Qubit id.
     * @return The outcome: false for eigenvalue +1, true for eigenvalue -1.
     */
    bool MeasureAndDiscardX(QubitId id);

    /**
     * Measure and discard qubit q in the Z basis.
     * @param id Qubit id.
     * @return The outcome: false for eigenvalue +1, true for eigenvalue -1.
     */
    bool MeasureAndDiscardZ(QubitId id);

    /**
     * Checks wether measuring qubit q in the X basis will be a deterministic measurement
     * @param id Qubit id.
     * @return The outcome: true if the measurement will be deterministic, false otherwise
     * @throws std::out_of_range If id is not active.
     */
    [[nodiscard]] bool IsDeterministicMeasureX(QubitId id) const;

    /**
     * Checks wether measuring qubit q in the Z basis will be a deterministic measurement
     * @param id Stable qubit id.
     * @return The outcome: true if the measurement will be deterministic, false otherwise
     * @throws std::out_of_range If id is not active.
     */
    [[nodiscard]] bool IsDeterministicMeasureZ(QubitId id) const;

    /**
     * Forces a non-deterministic measurement in the X basis to have a certain outcome
     * @param id Stable qubit id.
     * @param measurement the desired outcome of the measurement.
     * @return The outcome of the measurement
     */
    bool ForceMeasureX(QubitId id, bool measurement);

    /**
     * Forces a non-deterministic measurement in the Z basis to have a certain outcome
     * @param id Qubit id.
     * @param measurement the desired outcome of the measurement.
     * @return The outcome of the measurement
     */
    bool ForceMeasureZ(QubitId id, bool measurement);

    /**
     * Apply a Hadamard gate to qubit q.
     * @param id Qubit id.
     */
    void ApplyH(QubitId id);

    /**
     * Apply an S phase gate to qubit q.
     * @param id Qubit id.
     */
    void ApplyS(QubitId id);

    /**
     * Apply an S† gate to qubit q.
     * @param id Qubit id.
     */
    void ApplySdg(QubitId id);

    /**
     * Apply an √X gate to qubit q.
     * @param id Qubit id.
     */
    void ApplySx(QubitId id);

    /**
     * Apply the identity gate to qubit q.
     * @param id Qubit id.
     */
    void ApplyI(QubitId id);

    /**
     * Apply a Pauli X gate to qubit q.
     * @param id Qubit id.
     */
    void ApplyX(QubitId id);

    /**
     * Apply a Pauli Z gate to qubit q.
     * @param id Qubit id.
     */
    void ApplyZ(QubitId id);

    /**
     * Apply a Pauli Y gate to qubit q.
     * @param id Qubit id.
     */
    void ApplyY(QubitId id);

    /**
     * Apply a controlled X gate from control to target.
     * @param control Control qubit id.
     * @param target Target qubit id.
     */
    void ApplyCx(QubitId control, QubitId target);

    /**
     * Apply a controlled Z gate from control to target.
     * @param control Control qubit id.
     * @param target Target qubit id.
     */
    void ApplyCz(QubitId control, QubitId target);

    /**
     * Apply a controlled Y gate from control to target.
     * @param control Control qubit id.
     * @param target Target qubit id.
     */
    void ApplyCy(QubitId control, QubitId target);

    /**
     * Print stored generators: [a] active, [p] purification, [†] tombstone.
     */
    void Print() const;

    /**
     * Return the persistent bytes owned by this inverse tableau.
     * @return Persistent memory usage in bytes.
     */
    [[nodiscard]] std::size_t MemoryUsage() const noexcept;

    using UnaryGateFunction = void (stabaal::detail::InverseTableau::*)(stabaal::QubitId);

    using BinaryGateFunction = void (stabaal::detail::InverseTableau::*)(stabaal::QubitId,
                                                                         stabaal::QubitId);

    /**
     * Return a struct containing the name, arity and a description of a gate
     * @param gate Pointer to a unary gate function
     * @return information about the gate
     */
    [[nodiscard]] static GateInformation GetGateInformation(UnaryGateFunction gate);

    /**
     * Return a struct containing the name, arity and a description of a gate
     * @param gate Pointer to a binary gate function
     * @return Struct containing information about the gate
     */
    [[nodiscard]] static GateInformation GetGateInformation(BinaryGateFunction gate);

    /**
     * Return pointers to all functions that apply a unary gate
     * @return Array of function pointers to unary gate functions
     */
    [[nodiscard]] static std::array<UnaryGateFunction, 8> SupportedUnaryGates();

    /**
     * Return pointers to all functions that apply a binary gate
     * @return Array of function pointers to binary gate functions
     */
    [[nodiscard]] static std::array<BinaryGateFunction, 3> SupportedBinaryGates();

   private:
    std::unique_ptr<BitWord[]> words_;

    /**
     * Number of active qubits.
     */
    QubitSize num_qubits_;

    /** Stored slots, including tombstones and purification slots. */
    QubitSize num_slots_;

    /**
     * Number of qubits the tableau can currently hold.
     */
    QubitSize capacity_;

    /**
     * Sentinel value indicating that there is no free slots that can be reused.
     */
    static constexpr QubitSize kNoSlot{std::numeric_limits<QubitSize>::max()};

    /**
     * First free slot that can be reused.
     */
    QubitSize first_free_slot_{kNoSlot};

    /**
     * Number of free slots that can be reused.
     */
    QubitSize num_free_slots_{0};

    /**
     * Number of logical ID's created.
     */
    QubitSize num_ids_{0};

    /**
     * Number of ID's that fit in the mapping.
     */
    QubitSize id_capacity_{0};

    /**
     * First ID that is free to use.
     */
    QubitSize first_free_id_{kNoSlot};

    /**
     * ID generations and physical slot mapping when needed.
     */
    std::unique_ptr<std::uint32_t[]> generations_;

    /**
     * Highest ID generation issued, and remembered even after the last qubit is discarded.
     */
    std::uint64_t generation_epoch_{0};

    Sfc64 measurement_engine_;
    RandomBit<Sfc64::result_type> measurement_random_;

    /**
     * Validate a qubit ID and find its physical tableau slot.
     * @param id Qubit ID to resolve.
     * @return Physical tableau slot for the active qubit.
     * @throws std::out_of_range If the ID is not active.
     */
    [[nodiscard]] QubitSize ResolveQubit(QubitId id) const;

    /**
     * Check whether a discarded qubit must be retained for purification.
     * @param q Physical tableau slot.
     * @return True if X, Y, and Z measurements are all nondeterministic.
     */
    [[nodiscard]] bool RequiresPurification(QubitSize q) const noexcept;

    /**
     * Read an ID slot's generation.
     * @param id_slot Logical ID slot.
     * @return Current generation for the ID slot.
     */
    [[nodiscard]] std::uint32_t Generation(QubitSize id_slot) const noexcept;

    /**
     * Access an ID slot's stored generation for modification.
     * @param id_slot Logical ID slot.
     * @return Reference to the ID slot's generation.
     * @pre Generation storage has been allocated.
     */
    [[nodiscard]] std::uint32_t& Generation(QubitSize id_slot) noexcept;

    /**
     * Read an ID entry's physical slot or next free ID entry.
     * @param id_slot Logical ID slot.
     * @return Physical slot if live, or the next free ID slot if dead.
     */
    [[nodiscard]] QubitSize IdIndex(QubitSize id_slot) const noexcept;

    /**
     * Access an ID entry's index for modification.
     * @param id_slot Logical ID slot.
     * @return Reference to the physical slot or next free ID slot.
     */
    [[nodiscard]] QubitSize& IdIndex(QubitSize id_slot) noexcept;

    /**
     * Create the mapping from logical ID slots to physical slots.
     */
    void MapIds();

    /**
     * Ensure the ID mapping has room for another entry.
     */
    void ReserveId();

    /**
     * Find the next reusable physical slot and update the cached index.
     */
    void AdvanceFreeSlot() noexcept;

    /**
     * Append a physical slot initialized to |0>, growing storage if needed.
     * @return Index of the appended physical slot.
     */
    [[nodiscard]] QubitSize AppendSlot();

    /**
     * Reset a pure product qubit to |0> for reuse.
     * @param q Physical tableau slot.
     */
    void ResetProductSlot(QubitSize q) noexcept;

    /**
     * Apply a phase operation to a physical tableau slot.
     * @param q Physical tableau slot.
     * @param phase_rule Phase exponent used to update the generator sign.
     */
    void ApplyPhase(QubitSize q, PhaseExponent phase_rule) noexcept;

    /**
     * Remove reusable slots and shrink storage when the contraction threshold is met.
     */
    void Contract();

    /**
     * Remove a reset product slot while preserving the remaining tableau state.
     * @param q Physical tableau slot to remove.
     * @param targets Scratch space for packed target bits.
     * @param origins Original physical slot indices, updated when slots move.
     */
    void RemoveProductSlot(QubitSize q, BitWordSpan targets, std::span<QubitSize> origins) noexcept;

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
     * @param index Storage group: 0/1 signs, 2 swaps, 3 tombstones, 4 purification slots.
     * @return A view of the selected packed flags.
     */
    [[nodiscard]] BitWordSpan Flags(std::size_t index) noexcept;

    /**
     * Return a read-only view of packed flags.
     * @param index Storage group: 0/1 signs, 2 swaps, 3 tombstones, 4 purification slots.
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
     * The destination capacity must accommodate all of other's stored slots.
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
    [[nodiscard]] bool IsDeterministicMeasure(QubitSize q, Basis basis) const noexcept;

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

namespace stabaal {

using InverseTableau = detail::InverseTableau;

}  // namespace stabaal

#endif
