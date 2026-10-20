#ifndef STABAAL_SIMULATOR_HPP
#define STABAAL_SIMULATOR_HPP

#include <concepts>
#include <cstddef>

#include "detail/inverse_tableau.hpp"
#include "detail/types.hpp"  // IWYU pragma: export

namespace stabaal {

template <typename T>
concept SimulationBackend =
    std::constructible_from<T, QubitIdx> &&
    requires(T& backend, const T& const_backend, QubitIdx q, QubitIdx target) {
        { backend.Add() } -> std::same_as<void>;
        { backend.MeasureX(q) } -> std::same_as<bool>;
        { backend.MeasureZ(q) } -> std::same_as<bool>;
        { backend.ApplyH(q) } -> std::same_as<void>;
        { backend.ApplyS(q) } -> std::same_as<void>;
        { backend.ApplyX(q) } -> std::same_as<void>;
        { backend.ApplyZ(q) } -> std::same_as<void>;
        { backend.ApplyCx(q, target) } -> std::same_as<void>;
        { const_backend.Print() } -> std::same_as<void>;
        { const_backend.MemoryUsage() } noexcept -> std::same_as<std::size_t>;
    };

template <SimulationBackend Backend = detail::InverseTableau>
class Simulator {
   public:
    explicit Simulator(QubitIdx num_qubits) : backend_(num_qubits) {}

    void Add() { backend_.Add(); }

    bool MeasureX(QubitIdx q) { return backend_.MeasureX(q); }

    bool MeasureZ(QubitIdx q) { return backend_.MeasureZ(q); }

    void H(QubitIdx q) { backend_.ApplyH(q); }

    void S(QubitIdx q) { backend_.ApplyS(q); }

    void X(QubitIdx q) { backend_.ApplyX(q); }

    void Z(QubitIdx q) { backend_.ApplyZ(q); }

    void Cx(QubitIdx ctrl, QubitIdx trgt) { backend_.ApplyCx(ctrl, trgt); }

    void Print() const { backend_.Print(); }

    std::size_t MemoryUsage() const noexcept { return backend_.MemoryUsage(); }

   private:
    Backend backend_;
};
}  // namespace stabaal

#endif
