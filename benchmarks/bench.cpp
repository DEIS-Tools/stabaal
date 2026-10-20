#include <algorithm>
#include <charconv>
#include <chrono>
#include <cstddef>
#include <cstdint>
#include <exception>
#include <iomanip>
#include <iostream>
#include <memory>
#include <random>
#include <span>
#include <stdexcept>
#include <string_view>
#include <system_error>
#include <vector>

#include "inverse_tableau.hpp"
#include "types.hpp"

namespace {

template <typename T>
T parse(const std::string_view text) {
    T value{};
    const auto *const last{std::to_address(text.end())};
    const auto [end, error]{std::from_chars(std::to_address(text.begin()), last, value)};
    if (error != std::errc{} || end != last || value == 0) {
        throw std::invalid_argument{"Arguments must be positive integers"};
    }

    return value;
}

enum class GateKind : std::uint8_t { kH, kX, kS, kZ, kCx };

struct Gate {
    GateKind kind;
    stabaal::QubitIdx q;
    stabaal::QubitIdx target;
};

// Runs the circuit with the hand-written gates, or with the ones generated from gate matrices.
template <bool Generated>
void run(stabaal::InverseTableau &tableau, const std::vector<Gate> &circuit) {
    for (const auto &gate : circuit) {
        switch (gate.kind) {
            case GateKind::kH:
                Generated ? tableau.h(gate.q) : tableau.applyH(gate.q);
                break;
            case GateKind::kX:
                Generated ? tableau.x(gate.q) : tableau.applyX(gate.q);
                break;
            case GateKind::kS:
                Generated ? tableau.s(gate.q) : tableau.applyS(gate.q);
                break;
            case GateKind::kZ:
                Generated ? tableau.z(gate.q) : tableau.applyZ(gate.q);
                break;
            case GateKind::kCx:
                Generated ? tableau.cx(gate.q, gate.target) : tableau.applyCx(gate.q, gate.target);
                break;
        }
    }
}

double median(std::vector<double> &times) {
    std::sort(times.begin(), times.end());
    return (times[(times.size() - 1) / 2] + times[times.size() / 2]) / 2;
}

}  // namespace

int main(int argc, char **argv) {
    const std::span<char *const> args{argv, static_cast<std::size_t>(argc)};
    if (argc > 4 || (argc > 1 && std::string_view{args[1]} == "--help")) {
        std::cout << "Usage: stabaal_bench [qubits=256] [gates=1000000] [repeats=5]\n";
        return argc > 4 ? 1 : 0;
    }

    try {
        const auto qubits{argc > 1 ? parse<stabaal::QubitIdx>(args[1]) : stabaal::QubitIdx{256}};
        const auto gates{argc > 2 ? parse<std::size_t>(args[2]) : std::size_t{1000000}};
        const auto repeats{argc > 3 ? parse<std::uint32_t>(args[3]) : std::uint32_t{5}};
        std::uniform_int_distribution<stabaal::QubitIdx> qubit{0, qubits - 1};
        std::vector<Gate> circuit{};
        circuit.reserve(gates);

        std::vector<double> times{};
        std::vector<double> generated_times{};
        times.reserve(repeats);
        generated_times.reserve(repeats);
        std::uint32_t checksum{0};
        std::size_t memory_bytes{0};
        for (std::uint32_t repeat{0}; repeat < repeats; ++repeat) {
            std::mt19937 rng{1234U + repeat};
            circuit.clear();
            for (std::size_t i{0}; i < gates; ++i) {
                const auto q{qubit(rng)};
                auto target{qubit(rng)};
                if (target == q) target = (q + 1) % qubits;
                circuit.push_back(
                    {static_cast<GateKind>(rng() % (qubits == 1 ? 4 : 5)), q, target});
            }

            stabaal::InverseTableau tableau{qubits};
            auto start{std::chrono::steady_clock::now()};
            run<false>(tableau, circuit);
            auto stop{std::chrono::steady_clock::now()};
            times.push_back(std::chrono::duration<double>{stop - start}.count());
            if (repeat == 0) memory_bytes = tableau.memoryUsage();
            if (tableau.measureZ(repeat % qubits)) ++checksum;

            stabaal::InverseTableau generated{qubits};
            start = std::chrono::steady_clock::now();
            run<true>(generated, circuit);
            stop = std::chrono::steady_clock::now();
            generated_times.push_back(std::chrono::duration<double>{stop - start}.count());
            if (generated.measureZ(repeat % qubits)) ++checksum;
        }

        const auto report{[gates](const char *name, std::vector<double> &samples) {
            const auto seconds{median(samples)};
            std::cout << name << ": median=" << seconds << " s, "
                      << seconds * 1e9 / static_cast<double>(gates) << " ns/gate, "
                      << static_cast<double>(gates) / seconds / 1e6 << " million gates/s\n";
        }};

        std::cout << "qubits=" << qubits << " gates=" << gates << " repeats=" << repeats
                  << std::fixed << std::setprecision(3) << "\ntableau_memory=" << memory_bytes
                  << " bytes (" << static_cast<double>(memory_bytes) / 1024.0
                  << " KiB; excluding allocator overhead)\n";
        report("hand-written", times);
        report("generated   ", generated_times);
        std::cout << "checksum=" << checksum << '\n';
    } catch (const std::exception &error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
