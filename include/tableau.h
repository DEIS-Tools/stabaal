#pragma once
#ifndef TABLEAU_H

#include <iostream>
#include <iomanip>
#include <sstream>
#include <stdlib.h>
#include <time.h>

#include "generator.h"
    
using bit_ref = boost::dynamic_bitset<>::reference; 
using basis_table = std::vector<Generator>;

class Tableau {
   public: 
    enum Basis {
        COMPUTATIONAL,
        HADAMARD,
        IMAGINARY
    };

    Tableau(int qubits) 
        : computational(qubits, Generator(qubits)), 
          hadamard(qubits, Generator(qubits)), 
          qubits(qubits) {
     
        for (auto i = 0; i < qubits; ++i) {
            computational[i].z[i] = true; 
            hadamard[i].x[i]      = true; 
        }

        srand(time(NULL));
    }
    
    void add();
    
    bool measure_z(int q);
    bool measure_x(int q);
    
    void apply_h(int q);
    void apply_cx(int c, int t);
    void apply_x(int q);
    void apply_s(int q);
    void apply_z(int q);

    void print();
    
   private:
    basis_table computational;
    basis_table hadamard;    
    uint64_t qubits;

    void clean_up(int pivot, Generator gen);
    void collapse(int pivot, bool p);
    bool random_measurement(int q, Basis b);
    
    int e(Generator a, Generator b);
    bool phase(Generator a, Generator b, int phase_rule);

    void swap_bit(bit_ref a, bit_ref b);
};
#endif