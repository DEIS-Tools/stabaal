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
    : _computational(qubitCount, Generator(qubitCount)),
      _hadamard(qubitCount, Generator(qubitCount)),
      _qubits(qubitCount) {
    for (auto i = 0; i < _qubits; ++i) {
        _computational[i].z[i] = true;
        _hadamard[i].x[i] = true;
    }
}

std::mt19937 &Tableau::_randomEngine() {
    static thread_local std::mt19937 engine{std::random_device{}()};
    return engine;
}

void Tableau::add() {
    ++_qubits;

    for (auto &gen : _computational) {
        gen.x.push_back(false);
        gen.z.push_back(false);
    }
    for (auto &gen : _hadamard) {
        gen.x.push_back(false);
        gen.z.push_back(false);
    }
    _computational.emplace_back(_qubits);
    _hadamard.emplace_back(_qubits);

    const auto last = _qubits - 1;
    _computational[last].z[last] = true;
    _hadamard[last].x[last] = true;
}

void Tableau::_cleanUp(int pivot, const Generator &gen) {
    // cx between pivot and all columns with x in the measured row.
    for (auto i = pivot + 1; i < _qubits; ++i) {
        if (!gen.x[i]) continue;

        for (auto j = 0; j < _qubits; ++j) {
            _computational[j].r =
                _computational[j].r != (_computational[j].x[pivot] && _computational[j].z[i] &&
                                        (_computational[j].x[i] == _computational[j].z[pivot]));
            _computational[j].x[i] ^= _computational[j].x[pivot];
            _computational[j].z[pivot] ^= _computational[j].z[i];

            _hadamard[j].r = _hadamard[j].r != (_hadamard[j].x[pivot] && _hadamard[j].z[i] &&
                                                (_hadamard[j].x[i] == _hadamard[j].z[pivot]));
            _hadamard[j].x[i] ^= _hadamard[j].x[pivot];
            _hadamard[j].z[pivot] ^= _hadamard[j].z[i];
        }
    }
    if (!gen.z[pivot]) return;

    // s^dagger applied to pivot
    for (int i = 0; i < _qubits; ++i) {
        _hadamard[i].r = _hadamard[i].r != (_hadamard[i].x[pivot] && !_hadamard[i].z[pivot]);
        _hadamard[i].z[pivot] ^= _hadamard[i].x[pivot];

        _computational[i].r =
            _computational[i].r != (_computational[i].x[pivot] && !_computational[i].z[pivot]);
        _computational[i].z[pivot] ^= _computational[i].x[pivot];
    }
}

void Tableau::_collapse(int pivot, bool p) {
    // Apply X on the pivot when the measured outcome is 1 (p).
    for (int i = 0; i < _qubits; ++i) {
        _swapBit(_computational[i].x[pivot], _computational[i].z[pivot]);
        _computational[i].r = _computational[i].r !=
                              (_computational[i].z[pivot] && (p != _computational[i].x[pivot]));

        _swapBit(_hadamard[i].x[pivot], _hadamard[i].z[pivot]);
        _hadamard[i].r = _hadamard[i].r != (_hadamard[i].z[pivot] && (p != _hadamard[i].x[pivot]));
    }
}

bool Tableau::_randomMeasurement(int q, Basis b) {
    const auto &gen = (b == Basis::COMPUTATIONAL) ? _computational[q] : _hadamard[q];
    auto pivot = 0;

    for (auto i = 0; i < _qubits; ++i) {
        if (gen.x[i]) {
            pivot = i;
            break;
        }
    }
    _cleanUp(pivot, gen);
    const bool p = std::bernoulli_distribution{0.5}(_randomEngine());
    _collapse(pivot, p);
    return gen.r;
}

bool Tableau::measureZ(int q) {
    if (_computational[q].x.any()) {
        return _randomMeasurement(q, Basis::COMPUTATIONAL);
    }
    return _computational[q].r;
}

bool Tableau::measureX(int q) {
    if (_hadamard[q].x.any()) {
        return _randomMeasurement(q, Basis::HADAMARD);
    }
    return _hadamard[q].r;
}

void Tableau::applyH(int q) { std::swap(_computational[q], _hadamard[q]); }

void Tableau::applyCx(int c, int t) {
    _hadamard[c].r = _phase(_hadamard[c], _hadamard[t], 0);
    _hadamard[c].x ^= _hadamard[t].x;
    _hadamard[c].z ^= _hadamard[t].z;

    _computational[t].r = _phase(_computational[c], _computational[t], 0);
    _computational[t].x ^= _computational[c].x;
    _computational[t].z ^= _computational[c].z;
}

void Tableau::applyX(int q) { _computational[q].r = !_computational[q].r; }

void Tableau::applyS(int q) {
    _hadamard[q].r = _phase(_hadamard[q], _computational[q], 3);
    _hadamard[q].x ^= _computational[q].x;
    _hadamard[q].z ^= _computational[q].z;
}

void Tableau::applyZ(int q) { _hadamard[q].r = !_hadamard[q].r; }

void Tableau::print() const {
    for (size_t j = 0; j < _qubits; ++j) {
        std::cout << (j == 0 ? "\n    | " : "") << "x" << (j == _qubits - 1 ? " | " : " ");
    }
    for (size_t j = 0; j < _qubits; ++j) {
        std::cout << "z" << (j == _qubits - 1 ? " |" : " ");
    }
    std::cout << " r |";

    for (size_t i = 0; i < _qubits; ++i) {
        for (size_t j = 0; j < _qubits; ++j) {
            std::cout << (j == 0 ? "\nX" + std::to_string(i) + ": | " : "") << _hadamard[i].x[j]
                      << (j == _qubits - 1 ? " | " : " ");
        }
        for (size_t j = 0; j < _qubits; ++j) {
            std::cout << _hadamard[i].z[j] << (j == _qubits - 1 ? " | " : " ");
        }
        std::cout << _hadamard[i].r << " |";
    }

    for (size_t i = 0; i < _qubits; ++i) {
        for (size_t j = 0; j < _qubits; ++j) {
            std::cout << (j == 0 ? "\nZ" + std::to_string(i) + ": | " : "")
                      << _computational[i].x[j] << (j == _qubits - 1 ? " | " : " ");
        }
        for (size_t j = 0; j < _qubits; ++j) {
            std::cout << _computational[i].z[j] << (j == _qubits - 1 ? " | " : " ");
        }

        std::cout << _computational[i].r << " |";
    }
    std::cout << '\n';
}

void Tableau::_swapBit(BitRef a, BitRef b) {
    const bool temp = a;
    a = b;
    b = temp;
}

int64_t Tableau::_e(const Generator &a, const Generator &b) {
    return static_cast<int64_t>((a.z & a.x).count()) + static_cast<int64_t>((b.z & b.x).count()) +
           2 * static_cast<int64_t>((a.z & b.x).count()) -
           static_cast<int64_t>(((a.x ^ b.x) & (a.z ^ b.z)).count());
}

bool Tableau::_phase(const Generator &a, const Generator &b, int phaseRule) {
    const auto val =
        static_cast<uint64_t>(_e(a, b) + (a.r ? 2 : 0) + (b.r ? 2 : 0) + phaseRule) & 3U;
    if (val % 2 == 1) {
        throw std::runtime_error(
            "Final phase came out to uneven numbers, meaning there has been a mistake");
    }

    return val == 2;
}

}  // namespace stabaal
