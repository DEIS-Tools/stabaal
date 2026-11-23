#ifndef STABAAL_BACKENDS_SIMULATION_BACKEND_HPP
#define STABAAL_BACKENDS_SIMULATION_BACKEND_HPP

#include <concepts>
#include <cstddef>

#include "stabaal/types.hpp"

namespace stabaal {

/**
 * Operations required of a Simulator backend.
 */
template <typename T>
concept SimulationBackend = std::constructible_from<T, QubitSize> &&
                            requires(T& backend, const T& const_backend, QubitId q, QubitId control,
                                     QubitId target, bool measurement) {
                                { const_backend.NumQubits() } noexcept -> std::same_as<QubitSize>;
                                { backend.Add() } -> std::same_as<QubitId>;
                                { backend.Discard(q) } -> std::same_as<void>;
                                { backend.MeasureX(q) } -> std::same_as<bool>;
                                { backend.MeasureZ(q) } -> std::same_as<bool>;
                                { backend.MeasureAndDiscardX(q) } -> std::same_as<bool>;
                                { backend.MeasureAndDiscardZ(q) } -> std::same_as<bool>;
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
                                { backend.ApplyCx(control, target) } -> std::same_as<void>;
                                { backend.ApplyCz(control, target) } -> std::same_as<void>;
                                { backend.ApplyCy(control, target) } -> std::same_as<void>;
                                { const_backend.Print() } -> std::same_as<void>;
                                {
                                    const_backend.MemoryUsage()
                                } noexcept -> std::same_as<std::size_t>;
                            };

}  // namespace stabaal

#endif
