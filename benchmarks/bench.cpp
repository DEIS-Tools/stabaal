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

#include "tableau.hpp"
#include "types.hpp"

namespace {

template <typename T>
T Parse(std::string_view text) {
    T value{};
    const auto* const last{std::to_address(text.end())};
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

}  // namespace

int main(int argc, char** argv) {
    const std::span<char* const> args{argv, static_cast<std::size_t>(argc)};
    if (argc > 4 || (argc > 1 && std::string_view{args[1]} == "--help")) {
        std::cout << "Usage: stabaal_bench [qubits=256] [gates=1000000] [repeats=5]\n";
        return argc > 4 ? 1 : 0;
    }

    try {
        const auto qubits{argc > 1 ? Parse<stabaal::QubitIdx>(args[1]) : stabaal::QubitIdx{256}};
        const auto gates{argc > 2 ? Parse<std::size_t>(args[2]) : std::size_t{1000000}};
        const auto repeats{argc > 3 ? Parse<std::uint32_t>(args[3]) : std::uint32_t{5}};
        std::mt19937 rng{1234};
        std::uniform_int_distribution<stabaal::QubitIdx> qubit{0, qubits - 1};
        std::vector<Gate> circuit{};
        circuit.reserve(gates);
        for (std::size_t i{0}; i < gates; ++i) {
            const auto q{qubit(rng)};
            auto target{qubit(rng)};
            if (target == q) target = (q + 1) % qubits;
            circuit.push_back({static_cast<GateKind>(rng() % (qubits == 1 ? 4 : 5)), q, target});
        }

        std::vector<double> times{};
        times.reserve(repeats);
        std::uint32_t checksum{0};
        for (std::uint32_t repeat{0}; repeat < repeats; ++repeat) {
            stabaal::Tableau tableau{qubits};
            const auto start{std::chrono::steady_clock::now()};
            for (const auto& gate : circuit) {
                switch (gate.kind) {
                    case GateKind::kH:
                        tableau.ApplyH(gate.q);
                        break;
                    case GateKind::kX:
                        tableau.ApplyX(gate.q);
                        break;
                    case GateKind::kS:
                        tableau.ApplyS(gate.q);
                        break;
                    case GateKind::kZ:
                        tableau.ApplyZ(gate.q);
                        break;
                    case GateKind::kCx:
                        tableau.ApplyCx(gate.q, gate.target);
                        break;
                }
            }

            const auto stop{std::chrono::steady_clock::now()};
            times.push_back(std::chrono::duration<double>{stop - start}.count());
            if (tableau.MeasureZ(repeat % qubits)) ++checksum;
        }

        std::sort(times.begin(), times.end());
        const auto median{(times[(times.size() - 1) / 2] + times[times.size() / 2]) / 2};
        std::cout << "qubits=" << qubits << " gates=" << gates << " repeats=" << repeats
                  << std::fixed << std::setprecision(3) << "\nmedian=" << median << " s, "
                  << median * 1e9 / static_cast<double>(gates) << " ns/gate, "
                  << static_cast<double>(gates) / median / 1e6 << " million gates/s\n"
                  << "checksum=" << checksum << '\n';
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
