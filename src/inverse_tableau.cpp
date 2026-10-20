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
        computational_.gens[i].setZ(i);
        hadamard_.gens[i].setX(i);
    }
}

std::mt19937 &InverseTableau::randomEngine() {
    static thread_local std::mt19937 engine{std::random_device{}()};
    return engine;
}

void InverseTableau::add() {
    ++num_qubits_;
    computational_.add();
    hadamard_.add();

    const auto last{num_qubits_ - 1};
    computational_.gens[last].setZ(last);
    hadamard_.gens[last].setX(last);
}

void InverseTableau::cleanUp(const QubitIdx pivot, const PauliString &gen) noexcept {
    // cx between pivot and all columns with x in the measured row.
    for (std::size_t i{static_cast<std::size_t>(pivot) + 1}; i < num_qubits_; ++i) {
        if (!gen.getX(i)) continue;

        for (std::size_t j{0}; j < num_qubits_; ++j) {
            computational_.signs.set(
                j, computational_.signs.get(j) !=
                       (computational_.gens[j].getX(pivot) && computational_.gens[j].getZ(i) &&
                        (computational_.gens[j].getX(i) == computational_.gens[j].getZ(pivot))));
            computational_.gens[j].setX(
                i, computational_.gens[j].getX(i) != computational_.gens[j].getX(pivot));
            computational_.gens[j].setZ(
                pivot, computational_.gens[j].getZ(pivot) != computational_.gens[j].getZ(i));

            hadamard_.signs.set(j,
                                hadamard_.signs.get(j) !=
                                    (hadamard_.gens[j].getX(pivot) && hadamard_.gens[j].getZ(i) &&
                                     (hadamard_.gens[j].getX(i) == hadamard_.gens[j].getZ(pivot))));
            hadamard_.gens[j].setX(i, hadamard_.gens[j].getX(i) != hadamard_.gens[j].getX(pivot));
            hadamard_.gens[j].setZ(pivot,
                                   hadamard_.gens[j].getZ(pivot) != hadamard_.gens[j].getZ(i));
        }
    }

    if (!gen.getZ(pivot)) return;

    // s^dagger applied to pivot
    for (std::size_t i{0}; i < num_qubits_; ++i) {
        hadamard_.signs.set(i, hadamard_.signs.get(i) != (hadamard_.gens[i].getX(pivot) &&
                                                          !hadamard_.gens[i].getZ(pivot)));
        hadamard_.gens[i].setZ(pivot,
                               hadamard_.gens[i].getZ(pivot) != hadamard_.gens[i].getX(pivot));

        computational_.signs.set(
            i, computational_.signs.get(i) !=
                   (computational_.gens[i].getX(pivot) && !computational_.gens[i].getZ(pivot)));
        computational_.gens[i].setZ(
            pivot, computational_.gens[i].getZ(pivot) != computational_.gens[i].getX(pivot));
    }
}

void InverseTableau::collapse(const QubitIdx pivot, const bool p) noexcept {
    // Apply X on the pivot when the measured outcome is 1 (p).
    for (std::size_t i{0}; i < num_qubits_; ++i) {
        computational_.gens[i].swapXZ(pivot);
        computational_.signs.set(
            i, computational_.signs.get(i) != (computational_.gens[i].getZ(pivot) &&
                                               (p != computational_.gens[i].getX(pivot))));

        hadamard_.gens[i].swapXZ(pivot);
        hadamard_.signs.set(i, hadamard_.signs.get(i) != (hadamard_.gens[i].getZ(pivot) &&
                                                          (p != hadamard_.gens[i].getX(pivot))));
    }
}

bool InverseTableau::randomMeasurement(const QubitIdx q, const Basis b) {
    const auto &table{(b == Basis::kComputational) ? computational_ : hadamard_};
    const auto &gen{table.gens[q]};
    QubitIdx pivot{0};

    for (QubitIdx i{0}; i < num_qubits_; ++i) {
        if (gen.getX(i)) {
            pivot = i;
            break;
        }
    }

    cleanUp(pivot, gen);
    const bool p{std::bernoulli_distribution{0.5}(randomEngine())};
    collapse(pivot, p);
    return table.signs.get(q);
}

bool InverseTableau::measureZ(const QubitIdx q) {
    if (q >= num_qubits_) [[unlikely]] {
        throw std::out_of_range{"Qubit index out of range"};
    }

    if (computational_.gens[q].anyX()) {
        return randomMeasurement(q, Basis::kComputational);
    }

    return computational_.signs.get(q);
}

bool InverseTableau::measureX(const QubitIdx q) {
    if (q >= num_qubits_) [[unlikely]] {
        throw std::out_of_range{"Qubit index out of range"};
    }

    if (hadamard_.gens[q].anyX()) {
        return randomMeasurement(q, Basis::kHadamard);
    }

    return hadamard_.signs.get(q);
}

void InverseTableau::applyH(const QubitIdx q) {
    if (q >= num_qubits_) [[unlikely]] {
        throw std::out_of_range{"Qubit index out of range"};
    }

    std::swap(computational_.gens[q], hadamard_.gens[q]);
    const bool sign{computational_.signs.get(q)};
    computational_.signs.set(q, hadamard_.signs.get(q));
    hadamard_.signs.set(q, sign);
}

void InverseTableau::applyCx(const QubitIdx c, const QubitIdx t) {
    if (c >= num_qubits_ || t >= num_qubits_) [[unlikely]] {
        throw std::out_of_range{"Qubit index out of range"};
    }

    const auto hadamard_phase{hadamard_.gens[c].xorWithPhase(hadamard_.gens[t])};
    hadamard_.signs.set(c,
                        phase(hadamard_phase, hadamard_.signs.get(c), hadamard_.signs.get(t), 0));

    const auto computational_phase{computational_.gens[t].xorWithPhase(computational_.gens[c])};
    computational_.signs.set(
        t, phase(computational_phase, computational_.signs.get(t), computational_.signs.get(c), 0));
}

void InverseTableau::applyX(const QubitIdx q) {
    if (q >= num_qubits_) [[unlikely]] {
        throw std::out_of_range{"Qubit index out of range"};
    }

    computational_.signs.set(q, !computational_.signs.get(q));
}

void InverseTableau::applyS(const QubitIdx q) {
    if (q >= num_qubits_) [[unlikely]] {
        throw std::out_of_range{"Qubit index out of range"};
    }

    const auto exponent{hadamard_.gens[q].xorWithPhase(computational_.gens[q])};
    hadamard_.signs.set(q, phase(exponent, hadamard_.signs.get(q), computational_.signs.get(q), 3));
}

void InverseTableau::applyZ(const QubitIdx q) {
    if (q >= num_qubits_) [[unlikely]] {
        throw std::out_of_range{"Qubit index out of range"};
    }

    hadamard_.signs.set(q, !hadamard_.signs.get(q));
}

std::size_t InverseTableau::memoryUsage() const noexcept {
    auto bytes{sizeof(*this) + computational_.allocatedBytes() + hadamard_.allocatedBytes() +
               scratch_.capacity() * sizeof(PauliString)};
    for (const auto &row : scratch_) bytes += row.allocatedBytes();
    return bytes;
}

void InverseTableau::print() const {
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
                      << hadamard_.gens[i].getX(j) << (j == num_qubits_ - 1 ? " | " : " ");
        }

        for (std::size_t j{0}; j < num_qubits_; ++j) {
            std::cout << hadamard_.gens[i].getZ(j) << (j == num_qubits_ - 1 ? " | " : " ");
        }

        std::cout << hadamard_.signs.get(i) << " |";
    }

    for (std::size_t i{0}; i < num_qubits_; ++i) {
        for (std::size_t j{0}; j < num_qubits_; ++j) {
            std::cout << (j == 0 ? "\nZ" + std::to_string(i) + ": | " : "")
                      << computational_.gens[i].getX(j) << (j == num_qubits_ - 1 ? " | " : " ");
        }

        for (std::size_t j{0}; j < num_qubits_; ++j) {
            std::cout << computational_.gens[i].getZ(j) << (j == num_qubits_ - 1 ? " | " : " ");
        }

        std::cout << computational_.signs.get(i) << " |";
    }

    std::cout << '\n';
}

bool InverseTableau::phase(const std::int64_t exponent, const bool a_sign, const bool b_sign,
                           const std::uint8_t phase_rule) noexcept {
    assert(phase_rule < 4);
    const auto phase{static_cast<std::uint64_t>(exponent + phase_rule) & 3U};
    assert(phase % 2 == 0);
    return (phase == 2) != (a_sign != b_sign);
}
