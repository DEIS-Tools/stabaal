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
        { backend.add() } -> std::same_as<void>;
        { backend.measureX(q) } -> std::same_as<bool>;
        { backend.measureZ(q) } -> std::same_as<bool>;
        { backend.applyH(q) } -> std::same_as<void>;
        { backend.applyS(q) } -> std::same_as<void>;
        { backend.applyX(q) } -> std::same_as<void>;
        { backend.applyZ(q) } -> std::same_as<void>;
        { backend.applyCx(q, target) } -> std::same_as<void>;
        { const_backend.print() } -> std::same_as<void>;
        { const_backend.memoryUsage() } noexcept -> std::same_as<std::size_t>;
    };

template <SimulationBackend Backend = detail::InverseTableau>
class Simulator {
   public:
    explicit Simulator(QubitIdx num_qubits) : backend_(num_qubits) {}

    void add() { backend_.add(); }

    bool measureX(QubitIdx q) { return backend_.measureX(q); }

    bool measureZ(QubitIdx q) { return backend_.measureZ(q); }

    void h(QubitIdx q) { backend_.applyH(q); }

    void s(QubitIdx q) { backend_.applyS(q); }

    void x(QubitIdx q) { backend_.applyX(q); }

    void z(QubitIdx q) { backend_.applyZ(q); }

    void cx(QubitIdx ctrl, QubitIdx trgt) { backend_.applyCx(ctrl, trgt); }

    void print() const { backend_.print(); }

    std::size_t memoryUsage() const noexcept { return backend_.memoryUsage(); }

   private:
    Backend backend_;
};
}  // namespace stabaal

#endif
