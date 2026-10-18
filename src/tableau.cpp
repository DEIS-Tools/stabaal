#include "tableau.hpp"

#include <algorithm>
#include <cassert>
#include <cstddef>
#include <cstdint>
#include <iostream>
#include <random>
#include <string>

#include "pauli_string.hpp"

using namespace stabaal;

Tableau::Tableau(int num_qubits)
    : computational_(num_qubits, PauliString(num_qubits)),
      hadamard_(num_qubits, PauliString(num_qubits)),
      num_qubits_(num_qubits) {
    for (auto i = 0; i < num_qubits; ++i) {
        computational_[i].z_[i] = true;
        hadamard_[i].x_[i] = true;
    }
}

std::mt19937& Tableau::RandomEngine() {
    static thread_local std::mt19937 engine{std::random_device{}()};
    return engine;
}

void Tableau::Add() {
    ++num_qubits_;

    for (auto& gen : computational_) {
        gen.x_.push_back(false);
        gen.z_.push_back(false);
    }
    for (auto& gen : hadamard_) {
        gen.x_.push_back(false);
        gen.z_.push_back(false);
    }
    computational_.emplace_back(num_qubits_);
    hadamard_.emplace_back(num_qubits_);

    const auto last = num_qubits_ - 1;
    computational_[last].z_[last] = true;
    hadamard_[last].x_[last] = true;
}

void Tableau::CleanUp(int pivot, const PauliString& gen) noexcept {
    // cx between pivot and all columns with x in the measured row.
    for (auto i = pivot + 1; i < num_qubits_; ++i) {
        if (!gen.x_[i]) continue;

        for (auto j = 0; j < num_qubits_; ++j) {
            computational_[j].sign_ = computational_[j].sign_ !=
                                      (computational_[j].x_[pivot] && computational_[j].z_[i] &&
                                       (computational_[j].x_[i] == computational_[j].z_[pivot]));
            computational_[j].x_[i] ^= computational_[j].x_[pivot];
            computational_[j].z_[pivot] ^= computational_[j].z_[i];

            hadamard_[j].sign_ =
                hadamard_[j].sign_ != (hadamard_[j].x_[pivot] && hadamard_[j].z_[i] &&
                                       (hadamard_[j].x_[i] == hadamard_[j].z_[pivot]));
            hadamard_[j].x_[i] ^= hadamard_[j].x_[pivot];
            hadamard_[j].z_[pivot] ^= hadamard_[j].z_[i];
        }
    }
    if (!gen.z_[pivot]) return;

    // s^dagger applied to pivot
    for (int i = 0; i < num_qubits_; ++i) {
        hadamard_[i].sign_ =
            hadamard_[i].sign_ != (hadamard_[i].x_[pivot] && !hadamard_[i].z_[pivot]);
        hadamard_[i].z_[pivot] ^= hadamard_[i].x_[pivot];

        computational_[i].sign_ = computational_[i].sign_ !=
                                  (computational_[i].x_[pivot] && !computational_[i].z_[pivot]);
        computational_[i].z_[pivot] ^= computational_[i].x_[pivot];
    }
}

void Tableau::Collapse(int pivot, bool p) noexcept {
    // Apply X on the pivot when the measured outcome is 1 (p).
    for (int i = 0; i < num_qubits_; ++i) {
        SwapBit(computational_[i].x_[pivot], computational_[i].z_[pivot]);
        computational_[i].sign_ = computational_[i].sign_ != (computational_[i].z_[pivot] &&
                                                              (p != computational_[i].x_[pivot]));

        SwapBit(hadamard_[i].x_[pivot], hadamard_[i].z_[pivot]);
        hadamard_[i].sign_ =
            hadamard_[i].sign_ != (hadamard_[i].z_[pivot] && (p != hadamard_[i].x_[pivot]));
    }
}

bool Tableau::RandomMeasurement(int q, Basis b) {
    const auto& gen = (b == Basis::kComputational) ? computational_[q] : hadamard_[q];
    auto pivot = 0;

    for (auto i = 0; i < num_qubits_; ++i) {
        if (gen.x_[i]) {
            pivot = i;
            break;
        }
    }
    CleanUp(pivot, gen);
    const bool p = std::bernoulli_distribution{0.5}(RandomEngine());
    Collapse(pivot, p);
    return gen.sign_;
}

bool Tableau::MeasureZ(int q) {
    if (computational_[q].x_.any()) {
        return RandomMeasurement(q, Basis::kComputational);
    }
    return computational_[q].sign_;
}

bool Tableau::MeasureX(int q) {
    if (hadamard_[q].x_.any()) {
        return RandomMeasurement(q, Basis::kHadamard);
    }
    return hadamard_[q].sign_;
}

void Tableau::ApplyH(int q) { std::swap(computational_[q], hadamard_[q]); }

void Tableau::ApplyCx(int c, int t) {
    hadamard_[c].sign_ = Phase(hadamard_[c], hadamard_[t], 0);
    hadamard_[c].x_ ^= hadamard_[t].x_;
    hadamard_[c].z_ ^= hadamard_[t].z_;

    computational_[t].sign_ = Phase(computational_[c], computational_[t], 0);
    computational_[t].x_ ^= computational_[c].x_;
    computational_[t].z_ ^= computational_[c].z_;
}

void Tableau::ApplyX(int q) noexcept { computational_[q].sign_ = !computational_[q].sign_; }

void Tableau::ApplyS(int q) {
    hadamard_[q].sign_ = Phase(hadamard_[q], computational_[q], 3);
    hadamard_[q].x_ ^= computational_[q].x_;
    hadamard_[q].z_ ^= computational_[q].z_;
}

void Tableau::ApplyZ(int q) noexcept { hadamard_[q].sign_ = !hadamard_[q].sign_; }

void Tableau::Print() const {
    for (size_t j = 0; j < num_qubits_; ++j) {
        std::cout << (j == 0 ? "\n    | " : "") << "x" << (j == num_qubits_ - 1 ? " | " : " ");
    }
    for (size_t j = 0; j < num_qubits_; ++j) {
        std::cout << "z" << (j == num_qubits_ - 1 ? " |" : " ");
    }
    std::cout << " r |";

    for (size_t i = 0; i < num_qubits_; ++i) {
        for (size_t j = 0; j < num_qubits_; ++j) {
            std::cout << (j == 0 ? "\nX" + std::to_string(i) + ": | " : "") << hadamard_[i].x_[j]
                      << (j == num_qubits_ - 1 ? " | " : " ");
        }
        for (size_t j = 0; j < num_qubits_; ++j) {
            std::cout << hadamard_[i].z_[j] << (j == num_qubits_ - 1 ? " | " : " ");
        }
        std::cout << hadamard_[i].sign_ << " |";
    }

    for (size_t i = 0; i < num_qubits_; ++i) {
        for (size_t j = 0; j < num_qubits_; ++j) {
            std::cout << (j == 0 ? "\nZ" + std::to_string(i) + ": | " : "")
                      << computational_[i].x_[j] << (j == num_qubits_ - 1 ? " | " : " ");
        }
        for (size_t j = 0; j < num_qubits_; ++j) {
            std::cout << computational_[i].z_[j] << (j == num_qubits_ - 1 ? " | " : " ");
        }

        std::cout << computational_[i].sign_ << " |";
    }
    std::cout << '\n';
}

void Tableau::SwapBit(BitRef a, BitRef b) noexcept {
    const bool temp = a;
    a = b;
    b = temp;
}

std::int64_t Tableau::E(const PauliString& a, const PauliString& b) {
    return static_cast<std::int64_t>((a.z_ & a.x_).count()) +
           static_cast<std::int64_t>((b.z_ & b.x_).count()) +
           2 * static_cast<std::int64_t>((a.z_ & b.x_).count()) -
           static_cast<std::int64_t>(((a.x_ ^ b.x_) & (a.z_ ^ b.z_)).count());
}

bool Tableau::Phase(const PauliString& a, const PauliString& b, int phase_rule) {
    const auto val =
        static_cast<std::uint64_t>(E(a, b) + (a.sign_ ? 2 : 0) + (b.sign_ ? 2 : 0) + phase_rule) &
        3U;
    assert(val % 2 == 0);  // NOLINT(misc-static-assert)
    return val == 2;
}
