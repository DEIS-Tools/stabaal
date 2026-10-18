#include "inverse_tableau.hpp"

#include <algorithm>
#include <cassert>
#include <cstddef>
#include <cstdint>
#include <iostream>
#include <random>
#include <stdexcept>
#include <string>

#include "basis_table.hpp"
#include "pauli_string.hpp"
#include "types.hpp"

using namespace stabaal;

InverseTableau::InverseTableau(const QubitIdx num_qubits)
    : computational_{num_qubits}, hadamard_{num_qubits}, num_qubits_{num_qubits} {
    for (std::size_t i{0}; i < num_qubits; ++i) {
        computational_.gens[i].SetZ(i);
        hadamard_.gens[i].SetX(i);
    }
}

std::mt19937& InverseTableau::RandomEngine() {
    static thread_local std::mt19937 engine{std::random_device{}()};
    return engine;
}

void InverseTableau::Add() {
    ++num_qubits_;
    computational_.Add();
    hadamard_.Add();

    const auto last{num_qubits_ - 1};
    computational_.gens[last].SetZ(last);
    hadamard_.gens[last].SetX(last);
}

void InverseTableau::CleanUp(const QubitIdx pivot, const PauliString& gen) noexcept {
    // cx between pivot and all columns with x in the measured row.
    for (std::size_t i{static_cast<std::size_t>(pivot) + 1}; i < num_qubits_; ++i) {
        if (!gen.GetX(i)) continue;

        for (std::size_t j{0}; j < num_qubits_; ++j) {
            computational_.signs.Set(
                j, computational_.signs.Get(j) !=
                       (computational_.gens[j].GetX(pivot) && computational_.gens[j].GetZ(i) &&
                        (computational_.gens[j].GetX(i) == computational_.gens[j].GetZ(pivot))));
            computational_.gens[j].SetX(
                i, computational_.gens[j].GetX(i) != computational_.gens[j].GetX(pivot));
            computational_.gens[j].SetZ(
                pivot, computational_.gens[j].GetZ(pivot) != computational_.gens[j].GetZ(i));

            hadamard_.signs.Set(j,
                                hadamard_.signs.Get(j) !=
                                    (hadamard_.gens[j].GetX(pivot) && hadamard_.gens[j].GetZ(i) &&
                                     (hadamard_.gens[j].GetX(i) == hadamard_.gens[j].GetZ(pivot))));
            hadamard_.gens[j].SetX(i, hadamard_.gens[j].GetX(i) != hadamard_.gens[j].GetX(pivot));
            hadamard_.gens[j].SetZ(pivot,
                                   hadamard_.gens[j].GetZ(pivot) != hadamard_.gens[j].GetZ(i));
        }
    }

    if (!gen.GetZ(pivot)) return;

    // s^dagger applied to pivot
    for (std::size_t i{0}; i < num_qubits_; ++i) {
        hadamard_.signs.Set(i, hadamard_.signs.Get(i) != (hadamard_.gens[i].GetX(pivot) &&
                                                          !hadamard_.gens[i].GetZ(pivot)));
        hadamard_.gens[i].SetZ(pivot,
                               hadamard_.gens[i].GetZ(pivot) != hadamard_.gens[i].GetX(pivot));

        computational_.signs.Set(
            i, computational_.signs.Get(i) !=
                   (computational_.gens[i].GetX(pivot) && !computational_.gens[i].GetZ(pivot)));
        computational_.gens[i].SetZ(
            pivot, computational_.gens[i].GetZ(pivot) != computational_.gens[i].GetX(pivot));
    }
}

void InverseTableau::Collapse(const QubitIdx pivot, const bool p) noexcept {
    // Apply X on the pivot when the measured outcome is 1 (p).
    for (std::size_t i{0}; i < num_qubits_; ++i) {
        computational_.gens[i].SwapXZ(pivot);
        computational_.signs.Set(
            i, computational_.signs.Get(i) != (computational_.gens[i].GetZ(pivot) &&
                                               (p != computational_.gens[i].GetX(pivot))));

        hadamard_.gens[i].SwapXZ(pivot);
        hadamard_.signs.Set(i, hadamard_.signs.Get(i) != (hadamard_.gens[i].GetZ(pivot) &&
                                                          (p != hadamard_.gens[i].GetX(pivot))));
    }
}

bool InverseTableau::RandomMeasurement(const QubitIdx q, const Basis b) {
    const auto& table{(b == Basis::kComputational) ? computational_ : hadamard_};
    const auto& gen{table.gens[q]};
    QubitIdx pivot{0};

    for (QubitIdx i{0}; i < num_qubits_; ++i) {
        if (gen.GetX(i)) {
            pivot = i;
            break;
        }
    }

    CleanUp(pivot, gen);
    const bool p{std::bernoulli_distribution{0.5}(RandomEngine())};
    Collapse(pivot, p);
    return table.signs.Get(q);
}

bool InverseTableau::MeasureZ(const QubitIdx q) {
    if (q >= num_qubits_) [[unlikely]] {
        throw std::out_of_range{"Qubit index out of range"};
    }

    if (computational_.gens[q].AnyX()) {
        return RandomMeasurement(q, Basis::kComputational);
    }

    return computational_.signs.Get(q);
}

bool InverseTableau::MeasureX(const QubitIdx q) {
    if (q >= num_qubits_) [[unlikely]] {
        throw std::out_of_range{"Qubit index out of range"};
    }

    if (hadamard_.gens[q].AnyX()) {
        return RandomMeasurement(q, Basis::kHadamard);
    }

    return hadamard_.signs.Get(q);
}

void InverseTableau::ApplyH(const QubitIdx q) {
    if (q >= num_qubits_) [[unlikely]] {
        throw std::out_of_range{"Qubit index out of range"};
    }

    std::swap(computational_.gens[q], hadamard_.gens[q]);
    const bool sign{computational_.signs.Get(q)};
    computational_.signs.Set(q, hadamard_.signs.Get(q));
    hadamard_.signs.Set(q, sign);
}

void InverseTableau::ApplyCx(const QubitIdx c, const QubitIdx t) {
    if (c >= num_qubits_ || t >= num_qubits_) [[unlikely]] {
        throw std::out_of_range{"Qubit index out of range"};
    }

    const auto hadamard_phase{hadamard_.gens[c].XorWithPhase(hadamard_.gens[t])};
    hadamard_.signs.Set(c,
                        Phase(hadamard_phase, hadamard_.signs.Get(c), hadamard_.signs.Get(t), 0));

    const auto computational_phase{computational_.gens[t].XorWithPhase(computational_.gens[c])};
    computational_.signs.Set(
        t, Phase(computational_phase, computational_.signs.Get(t), computational_.signs.Get(c), 0));
}

void InverseTableau::ApplyX(const QubitIdx q) {
    if (q >= num_qubits_) [[unlikely]] {
        throw std::out_of_range{"Qubit index out of range"};
    }

    computational_.signs.Set(q, !computational_.signs.Get(q));
}

void InverseTableau::ApplyS(const QubitIdx q) {
    if (q >= num_qubits_) [[unlikely]] {
        throw std::out_of_range{"Qubit index out of range"};
    }

    const auto exponent{hadamard_.gens[q].XorWithPhase(computational_.gens[q])};
    hadamard_.signs.Set(q, Phase(exponent, hadamard_.signs.Get(q), computational_.signs.Get(q), 3));
}

void InverseTableau::ApplyZ(const QubitIdx q) {
    if (q >= num_qubits_) [[unlikely]] {
        throw std::out_of_range{"Qubit index out of range"};
    }

    hadamard_.signs.Set(q, !hadamard_.signs.Get(q));
}

std::size_t InverseTableau::MemoryUsage() const noexcept {
    return sizeof(*this) + computational_.AllocatedBytes() + hadamard_.AllocatedBytes();
}

void InverseTableau::Print() const {
    for (std::size_t j{0}; j < num_qubits_; ++j) {
        std::cout << (j == 0 ? "\n    | " : "") << "x" << (j == num_qubits_ - 1 ? " | " : " ");
    }

    for (std::size_t j{0}; j < num_qubits_; ++j) {
        std::cout << "z" << (j == num_qubits_ - 1 ? " |" : " ");
    }

    std::cout << " r |";

    for (std::size_t i{0}; i < num_qubits_; ++i) {
        for (std::size_t j{0}; j < num_qubits_; ++j) {
            std::cout << (j == 0 ? "\nX" + std::to_string(i) + ": | " : "")
                      << hadamard_.gens[i].GetX(j) << (j == num_qubits_ - 1 ? " | " : " ");
        }

        for (std::size_t j{0}; j < num_qubits_; ++j) {
            std::cout << hadamard_.gens[i].GetZ(j) << (j == num_qubits_ - 1 ? " | " : " ");
        }

        std::cout << hadamard_.signs.Get(i) << " |";
    }

    for (std::size_t i{0}; i < num_qubits_; ++i) {
        for (std::size_t j{0}; j < num_qubits_; ++j) {
            std::cout << (j == 0 ? "\nZ" + std::to_string(i) + ": | " : "")
                      << computational_.gens[i].GetX(j) << (j == num_qubits_ - 1 ? " | " : " ");
        }

        for (std::size_t j{0}; j < num_qubits_; ++j) {
            std::cout << computational_.gens[i].GetZ(j) << (j == num_qubits_ - 1 ? " | " : " ");
        }

        std::cout << computational_.signs.Get(i) << " |";
    }

    std::cout << '\n';
}

bool InverseTableau::Phase(const std::int64_t exponent, const bool a_sign, const bool b_sign,
                           const std::uint8_t phase_rule) noexcept {
    assert(phase_rule < 4);
    const auto phase{static_cast<std::uint64_t>(exponent + phase_rule) & 3U};
    assert(phase % 2 == 0);
    return (phase == 2) != (a_sign != b_sign);
}
