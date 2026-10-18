#include "tableau.hpp"

#include <algorithm>
#include <cassert>
#include <cstddef>
#include <cstdint>
#include <iostream>
#include <random>
#include <string>

#include "basis_table.hpp"
#include "pauli_string.hpp"
#include "types.hpp"

using namespace stabaal;

Tableau::Tableau(QubitIdx num_qubits)
    : computational_{num_qubits}, hadamard_{num_qubits}, num_qubits_{num_qubits} {
    for (std::size_t i{0}; i < num_qubits; ++i) {
        computational_.gens[i].SetZ(i);
        hadamard_.gens[i].SetX(i);
    }
}

std::mt19937& Tableau::RandomEngine() {
    static thread_local std::mt19937 engine{std::random_device{}()};
    return engine;
}

void Tableau::Add() {
    ++num_qubits_;
    computational_.Add();
    hadamard_.Add();

    const auto last{num_qubits_ - 1};
    computational_.gens[last].SetZ(last);
    hadamard_.gens[last].SetX(last);
}

void Tableau::CleanUp(QubitIdx pivot, const PauliString& gen) noexcept {
    // cx between pivot and all columns with x in the measured row.
    for (std::size_t i{static_cast<std::size_t>(pivot) + 1}; i < num_qubits_; ++i) {
        if (!gen.GetX(i)) continue;

        for (std::size_t j{0}; j < num_qubits_; ++j) {
            computational_.signs[j] =
                computational_.signs[j] !=
                (computational_.gens[j].GetX(pivot) && computational_.gens[j].GetZ(i) &&
                 (computational_.gens[j].GetX(i) == computational_.gens[j].GetZ(pivot)));
            computational_.gens[j].SetX(
                i, computational_.gens[j].GetX(i) != computational_.gens[j].GetX(pivot));
            computational_.gens[j].SetZ(
                pivot, computational_.gens[j].GetZ(pivot) != computational_.gens[j].GetZ(i));

            hadamard_.signs[j] = hadamard_.signs[j] !=
                                 (hadamard_.gens[j].GetX(pivot) && hadamard_.gens[j].GetZ(i) &&
                                  (hadamard_.gens[j].GetX(i) == hadamard_.gens[j].GetZ(pivot)));
            hadamard_.gens[j].SetX(i, hadamard_.gens[j].GetX(i) != hadamard_.gens[j].GetX(pivot));
            hadamard_.gens[j].SetZ(pivot,
                                   hadamard_.gens[j].GetZ(pivot) != hadamard_.gens[j].GetZ(i));
        }
    }
    if (!gen.GetZ(pivot)) return;

    // s^dagger applied to pivot
    for (std::size_t i{0}; i < num_qubits_; ++i) {
        hadamard_.signs[i] =
            hadamard_.signs[i] != (hadamard_.gens[i].GetX(pivot) && !hadamard_.gens[i].GetZ(pivot));
        hadamard_.gens[i].SetZ(pivot,
                               hadamard_.gens[i].GetZ(pivot) != hadamard_.gens[i].GetX(pivot));

        computational_.signs[i] = computational_.signs[i] != (computational_.gens[i].GetX(pivot) &&
                                                              !computational_.gens[i].GetZ(pivot));
        computational_.gens[i].SetZ(
            pivot, computational_.gens[i].GetZ(pivot) != computational_.gens[i].GetX(pivot));
    }
}

void Tableau::Collapse(QubitIdx pivot, bool p) noexcept {
    // Apply X on the pivot when the measured outcome is 1 (p).
    for (std::size_t i{0}; i < num_qubits_; ++i) {
        computational_.gens[i].SwapXZ(pivot);
        computational_.signs[i] =
            computational_.signs[i] !=
            (computational_.gens[i].GetZ(pivot) && (p != computational_.gens[i].GetX(pivot)));

        hadamard_.gens[i].SwapXZ(pivot);
        hadamard_.signs[i] = hadamard_.signs[i] != (hadamard_.gens[i].GetZ(pivot) &&
                                                    (p != hadamard_.gens[i].GetX(pivot)));
    }
}

bool Tableau::RandomMeasurement(QubitIdx q, Basis b) {
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
    return table.signs[q];
}

bool Tableau::MeasureZ(QubitIdx q) {
    if (computational_.gens[q].AnyX()) {
        return RandomMeasurement(q, Basis::kComputational);
    }
    return computational_.signs[q];
}

bool Tableau::MeasureX(QubitIdx q) {
    if (hadamard_.gens[q].AnyX()) {
        return RandomMeasurement(q, Basis::kHadamard);
    }
    return hadamard_.signs[q];
}

void Tableau::ApplyH(QubitIdx q) {
    std::swap(computational_.gens[q], hadamard_.gens[q]);
    const bool sign{computational_.signs[q]};
    computational_.signs[q] = hadamard_.signs[q];
    hadamard_.signs[q] = sign;
}

void Tableau::ApplyCx(QubitIdx c, QubitIdx t) {
    hadamard_.signs[c] =
        Phase(hadamard_.gens[c], hadamard_.signs[c], hadamard_.gens[t], hadamard_.signs[t], 0);
    hadamard_.gens[c].XorWith(hadamard_.gens[t]);

    computational_.signs[t] = Phase(computational_.gens[c], computational_.signs[c],
                                    computational_.gens[t], computational_.signs[t], 0);
    computational_.gens[t].XorWith(computational_.gens[c]);
}

void Tableau::ApplyX(QubitIdx q) noexcept { computational_.signs[q] = !computational_.signs[q]; }

void Tableau::ApplyS(QubitIdx q) {
    hadamard_.signs[q] = Phase(hadamard_.gens[q], hadamard_.signs[q], computational_.gens[q],
                               computational_.signs[q], 3);
    hadamard_.gens[q].XorWith(computational_.gens[q]);
}

void Tableau::ApplyZ(QubitIdx q) noexcept { hadamard_.signs[q] = !hadamard_.signs[q]; }

void Tableau::Print() const {
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
        std::cout << hadamard_.signs[i] << " |";
    }

    for (std::size_t i{0}; i < num_qubits_; ++i) {
        for (std::size_t j{0}; j < num_qubits_; ++j) {
            std::cout << (j == 0 ? "\nZ" + std::to_string(i) + ": | " : "")
                      << computational_.gens[i].GetX(j) << (j == num_qubits_ - 1 ? " | " : " ");
        }
        for (std::size_t j{0}; j < num_qubits_; ++j) {
            std::cout << computational_.gens[i].GetZ(j) << (j == num_qubits_ - 1 ? " | " : " ");
        }

        std::cout << computational_.signs[i] << " |";
    }
    std::cout << '\n';
}

bool Tableau::Phase(const PauliString& a, bool a_sign, const PauliString& b, bool b_sign,
                    std::uint8_t phase_rule) {
    assert(phase_rule < 4);  // NOLINT(misc-static-assert)
    const auto val{
        static_cast<std::uint64_t>(a.PhaseExponent(b) + (a_sign ? 2 : 0) + (b_sign ? 2 : 0) +
                                   phase_rule) &
            3U,
    };
    assert(val % 2 == 0);  // NOLINT(misc-static-assert)
    return val == 2;
}
