#include "generator.hpp"

namespace stabaal {

Generator::Generator() = default;

Generator::Generator(int qubits) : x(qubits), z(qubits) {}

}  // namespace stabaal
