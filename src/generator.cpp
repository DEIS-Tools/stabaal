#include "generator.hpp"

using namespace stabaal;

Generator::Generator() = default;
Generator::Generator(int qubits) : x(qubits), z(qubits) {}
