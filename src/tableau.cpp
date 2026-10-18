#include "tableau.hpp"

#include <algorithm>
#include <cstddef>
#include <iostream>
#include <random>
#include <stdexcept>
#include <string>

#include "generator.hpp"

namespace stabaal {

Tableau::Tableau(int qubitCount)
    : computational_(qubitCount, Generator(qubitCount)),
      hadamard_(qubitCount, Generator(qubitCount)),
      qubits_(qubitCount) {
    for (auto i = 0; i < qubits_; ++i) {
        computational_[i].z[i] = true;
        hadamard_[i].x[i] = true;
    }
}

std::mt19937& Tableau::randomEngine_() {
    static thread_local std::mt19937 engine{std::random_device{}()};
    return engine;
}

void Tableau::add() {
    ++qubits_;

    for (auto& gen : computational_) {
        gen.x.push_back(false);
        gen.z.push_back(false);
    }
    for (auto& gen : hadamard_) {
        gen.x.push_back(false);
        gen.z.push_back(false);
    }
    computational_.emplace_back(qubits_);
    hadamard_.emplace_back(qubits_);

    const auto last = qubits_ - 1;
    computational_[last].z[last] = true;
    hadamard_[last].x[last] = true;
}

void Tableau::cleanUp_(int pivot, const Generator& gen) {
    // cx between pivot and all columns with x in the measured row.
    for (auto i = pivot + 1; i < qubits_; ++i) {
        if (!gen.x[i]) continue;

        for (auto j = 0; j < qubits_; ++j) {
            computational_[j].r =
                computational_[j].r != (computational_[j].x[pivot] && computational_[j].z[i] &&
                                        (computational_[j].x[i] == computational_[j].z[pivot]));
            computational_[j].x[i] ^= computational_[j].x[pivot];
            computational_[j].z[pivot] ^= computational_[j].z[i];

            hadamard_[j].r = hadamard_[j].r != (hadamard_[j].x[pivot] && hadamard_[j].z[i] &&
                                                (hadamard_[j].x[i] == hadamard_[j].z[pivot]));
            hadamard_[j].x[i] ^= hadamard_[j].x[pivot];
            hadamard_[j].z[pivot] ^= hadamard_[j].z[i];
        }
    }
    if (!gen.z[pivot]) return;

    // s^dagger applied to pivot
    for (int i = 0; i < qubits_; ++i) {
        hadamard_[i].r = hadamard_[i].r != (hadamard_[i].x[pivot] && !hadamard_[i].z[pivot]);
        hadamard_[i].z[pivot] ^= hadamard_[i].x[pivot];

        computational_[i].r =
            computational_[i].r != (computational_[i].x[pivot] && !computational_[i].z[pivot]);
        computational_[i].z[pivot] ^= computational_[i].x[pivot];
    }
}

void Tableau::collapse_(int pivot, bool p) {
    // Apply X on the pivot when the measured outcome is 1 (p).
    for (int i = 0; i < qubits_; ++i) {
        swapBit_(computational_[i].x[pivot], computational_[i].z[pivot]);
        computational_[i].r = computational_[i].r !=
                              (computational_[i].z[pivot] && (p != computational_[i].x[pivot]));

        swapBit_(hadamard_[i].x[pivot], hadamard_[i].z[pivot]);
        hadamard_[i].r = hadamard_[i].r != (hadamard_[i].z[pivot] && (p != hadamard_[i].x[pivot]));
    }
}

bool Tableau::randomMeasurement_(int q, Basis b) {
    const auto& gen = (b == Basis::COMPUTATIONAL) ? computational_[q] : hadamard_[q];
    auto pivot = 0;

    for (auto i = 0; i < qubits_; ++i) {
        if (gen.x[i]) {
            pivot = i;
            break;
        }
    }
    cleanUp_(pivot, gen);
    const bool p = std::bernoulli_distribution{0.5}(randomEngine_());
    collapse_(pivot, p);
    return gen.r;
}

bool Tableau::measureZ(int q) {
    if (computational_[q].x.any()) {
        return randomMeasurement_(q, Basis::COMPUTATIONAL);
    }
    return computational_[q].r;
}

bool Tableau::measureX(int q) {
    if (hadamard_[q].x.any()) {
        return randomMeasurement_(q, Basis::HADAMARD);
    }
    return hadamard_[q].r;
}

void Tableau::applyH(int q) { std::swap(computational_[q], hadamard_[q]); }

void Tableau::applyCx(int c, int t) {
    hadamard_[c].r = phase_(hadamard_[c], hadamard_[t], 0);
    hadamard_[c].x ^= hadamard_[t].x;
    hadamard_[c].z ^= hadamard_[t].z;

    computational_[t].r = phase_(computational_[c], computational_[t], 0);
    computational_[t].x ^= computational_[c].x;
    computational_[t].z ^= computational_[c].z;
}

void Tableau::applyX(int q) { computational_[q].r = !computational_[q].r; }

void Tableau::applyS(int q) {
    hadamard_[q].r = phase_(hadamard_[q], computational_[q], 3);
    hadamard_[q].x ^= computational_[q].x;
    hadamard_[q].z ^= computational_[q].z;
}

void Tableau::applyZ(int q) { hadamard_[q].r = !hadamard_[q].r; }

void Tableau::print() const {
    for (size_t j = 0; j < qubits_; ++j) {
        std::cout << (j == 0 ? "\n    | " : "") << "x" << (j == qubits_ - 1 ? " | " : " ");
    }
    for (size_t j = 0; j < qubits_; ++j) {
        std::cout << "z" << (j == qubits_ - 1 ? " |" : " ");
    }
    std::cout << " r |";

    for (size_t i = 0; i < qubits_; ++i) {
        for (size_t j = 0; j < qubits_; ++j) {
            std::cout << (j == 0 ? "\nX" + std::to_string(i) + ": | " : "") << hadamard_[i].x[j]
                      << (j == qubits_ - 1 ? " | " : " ");
        }
        for (size_t j = 0; j < qubits_; ++j) {
            std::cout << hadamard_[i].z[j] << (j == qubits_ - 1 ? " | " : " ");
        }
        std::cout << hadamard_[i].r << " |";
    }

    for (size_t i = 0; i < qubits_; ++i) {
        for (size_t j = 0; j < qubits_; ++j) {
            std::cout << (j == 0 ? "\nZ" + std::to_string(i) + ": | " : "")
                      << computational_[i].x[j] << (j == qubits_ - 1 ? " | " : " ");
        }
        for (size_t j = 0; j < qubits_; ++j) {
            std::cout << computational_[i].z[j] << (j == qubits_ - 1 ? " | " : " ");
        }

        std::cout << computational_[i].r << " |";
    }
    std::cout << '\n';
}

void Tableau::swapBit_(BitRef a, BitRef b) {
    const bool temp = a;
    a = b;
    b = temp;
}

int64_t Tableau::e_(const Generator& a, const Generator& b) {
    return static_cast<int64_t>((a.z & a.x).count()) + static_cast<int64_t>((b.z & b.x).count()) +
           2 * static_cast<int64_t>((a.z & b.x).count()) -
           static_cast<int64_t>(((a.x ^ b.x) & (a.z ^ b.z)).count());
}

bool Tableau::phase_(const Generator& a, const Generator& b, int phaseRule) {
    const auto val =
        static_cast<uint64_t>(e_(a, b) + (a.r ? 2 : 0) + (b.r ? 2 : 0) + phaseRule) & 3U;
    if (val % 2 == 1) {
        throw std::runtime_error(
            "Final phase came out to uneven numbers, meaning there has been a mistake");
    }

    return val == 2;
}

}  // namespace stabaal
