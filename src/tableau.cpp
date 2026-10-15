#include "tableau.h"

void Tableau::add() {
    ++qubits;
    
    for (auto gen : computational) {
        gen.x.append(false);
        gen.z.append(false);
    }
    for (auto gen : hadamard) {
        gen.x.append(false);
        gen.z.append(false);
    }
    computational.emplace_back(qubits);
    hadamard.emplace_back(qubits);

    auto last = qubits - 1;
    computational[last].z[last] = true;
    hadamard[last].x[last]      = true;
}

void Tableau::clean_up(int pivot, Generator gen) {
    for (auto i = pivot+1; i < qubits; ++i) {
        if (!gen.x[i]) continue;

        for (auto j = 0; j < qubits; ++j) {
            computational[j].r        ^= computational[j].x[pivot] * computational[j].z[i] * (computational[j].x[i] ^ computational[j].z[pivot] ^ 1);
            computational[j].x[i]     ^= computational[j].x[pivot];
            computational[j].z[pivot] ^= computational[j].z[i];

            hadamard[j].r        ^= hadamard[j].x[pivot] * hadamard[j].z[i] * (hadamard[j].x[i] ^ hadamard[j].z[pivot] ^ 1);
            hadamard[j].x[i]     ^= hadamard[j].x[pivot];
            hadamard[j].z[pivot] ^= hadamard[j].z[i];
        }
    }
    if (!gen.z[pivot]) return;
    
    // s^dagger applied to pivot
    gen.r        ^= gen.x[pivot] * !gen.z[pivot];
    gen.z[pivot] ^= gen.x[pivot];
}

void Tableau::collapse(int pivot, bool p) {
    // The terinary operators are x applied on pivot if measured 1 (p).
    for (int i = 0; i < qubits; ++i) {
        swap_bit(computational[i].x[pivot], computational[i].z[pivot]);
        computational[i].r ^= ((p) ? computational[i].z[pivot] : false) ^ (computational[i].z[pivot] & computational[i].x[pivot]);

        swap_bit(hadamard[i].x[pivot], hadamard[i].z[pivot]);
        hadamard[i].r ^= ((p) ? hadamard[i].z[pivot] : false) ^ (hadamard[i].z[pivot] & hadamard[i].x[pivot]);
    }
}

bool Tableau::random_measurement(int q, Basis b) {
    auto gen = (b == COMPUTATIONAL) ? computational[q] : hadamard[q];
    auto pivot = 0;

    for (auto i = 0; i < qubits; ++i) {
        if (gen.x[i]) { pivot = i; break; }
    }
    clean_up(pivot, gen);
    auto p = rand() % 2;
    collapse(pivot, p);
    return p;
}

bool Tableau::measure_z(int q) {
    if (computational[q].x.any()) {
        return random_measurement(q, COMPUTATIONAL);
    } 
    return computational[q].r;
}

bool Tableau::measure_x(int q) {
    if (hadamard[q].x.any()) {
        return random_measurement(q, HADAMARD);
    } 
    return hadamard[q].r;
}

void Tableau::apply_h(int q) {
    auto temp        = computational[q];
    computational[q] = hadamard[q];
    hadamard[q]      = temp;
}

void Tableau::apply_cx(int c, int t) {
    hadamard[c].r  = phase(hadamard[c], hadamard[t], 0);
    hadamard[c].x ^= hadamard[t].x;
    hadamard[c].z ^= hadamard[t].z;

    computational[c].r  = phase(computational[c], computational[t], 0);
    computational[t].x ^= computational[c].x;
    computational[t].z ^= computational[c].z;
}

void Tableau::apply_x(int q) {
    computational[q].r = !computational[q].r;
}

void Tableau::apply_s(int q) {
    hadamard[q].r  = phase(hadamard[q], computational[q], 3);
    hadamard[q].x ^= computational[q].x;
    hadamard[q].z ^= computational[q].z;
}

void Tableau::apply_z(int q) {
    hadamard[q].r = !hadamard[q].r;
}

void Tableau::print() {
    for (size_t j = 0; j < qubits; ++j)
    std::cout << (j == 0 ? "\n    | " : "") << "x" << (j == qubits - 1 ? " | " : " ");
    for (size_t j = 0; j < qubits; ++j)
    std::cout << "z" << (j == qubits - 1 ? " |" : " ");
    std::cout << " r |";
    
    for (size_t i = 0; i < qubits; ++i) {
        for (size_t j = 0; j < qubits; ++j)
            std::cout << (j == 0 ? "\nX" + std::to_string(i) + ": | " : "") << hadamard[i].x[j] << (j == qubits - 1 ? " | " : " ");
        for (size_t j = 0; j < qubits; ++j)
            std::cout << hadamard[i].z[j] << (j == qubits - 1 ? " | " : " ");
        std::cout << hadamard[i].r << " |";
    }

    for (size_t i = 0; i < qubits; ++i) {
        for (size_t j = 0; j < qubits; ++j)
            std::cout << (j == 0 ? "\nZ" + std::to_string(i) + ": | " : "") << computational[i].x[j] << (j == qubits - 1 ? " | " : " ");
        for (size_t j = 0; j < qubits; ++j)
            std::cout << computational[i].z[j] << (j == qubits - 1 ? " | " : " ");
        std::cout << computational[i].r << " |";
    }
    std::cout << '\n';
}

void Tableau::swap_bit(bit_ref a, bit_ref b) {
    bool temp = a;
    a = b;
    b = temp;
}

int Tableau::e(Generator a, Generator b) {
    return (a.z & a.x).count() + (b.z & b.x).count() + 2 * (a.z & b.x).count() - ((a.x ^ b.x) & (a.z ^ b.z)).count();
}

bool Tableau::phase(Generator a, Generator b, int phase_rule) {
    int val = ((e(a, b) + 2 * a.r + 2 * b.r + phase_rule) % 4 + 4);
    if (val % 2 == 1) 
        throw std::runtime_error("Final phase came out to uneven numbers, meaning there has been a mistake");
    return val == 2;
}
