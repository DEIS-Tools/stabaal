#include <algorithm>
#include <array>
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

#include "stabaal.hpp"

namespace {

using BenchSimulator = stabaal::Simulator<>;

template <typename T>
[[nodiscard]] T Parse(const std::string_view text) {
    T value{};
    const auto* const last{std::to_address(text.end())};
    const auto [end, error]{std::from_chars(std::to_address(text.begin()), last, value)};
    if (error != std::errc{} || end != last || value == 0) {
        throw std::invalid_argument{"Arguments must be positive integers"};
    }

    return value;
}

void PrintTiming(const std::string_view name, std::vector<double>& samples,
                 const std::size_t operation_count) {
    std::sort(samples.begin(), samples.end());
    const auto median{
        (samples[(samples.size() - 1) / 2] + samples[samples.size() / 2]) / 2,
    };
    const auto operations{static_cast<double>(operation_count)};
    auto duration{median * 1e9 / operations};
    constexpr std::array<std::string_view, 4> units{"ns", "µs", "ms", "s"};
    std::size_t unit{0};
    while (duration >= 1000 && unit + 1 < units.size()) {
        duration /= 1000;
        ++unit;
    }

    std::cout << name << ": " << duration << ' ' << units[unit] << "/op, "
              << operations / median / 1e6 << " million ops/s\n";
}

enum class GateKind : std::uint8_t { kH, kX, kS, kSdg, kZ, kY, kCx, kCz, kCy };

struct Gate {
    GateKind kind;
    stabaal::QubitSize q;
    stabaal::QubitSize target;
};

}  // namespace

int main(int argc, char** argv) {
    const std::span<char* const> args{argv, static_cast<std::size_t>(argc)};
    const bool help{argc == 2 && std::string_view{args[1]} == "--help"};
    if (argc != 4 || help) {
        std::cout << "Usage: stabaal_bench <qubits> <gates> <repeats>\n";
        return help ? 0 : 1;
    }

    try {
        const auto qubits{Parse<stabaal::QubitSize>(args[1])};
        const auto gates{Parse<std::size_t>(args[2])};
        const auto repeats{Parse<std::uint32_t>(args[3])};
        std::uniform_int_distribution<stabaal::QubitSize> qubit{0, qubits - 1};
        std::vector<Gate> circuit{};
        circuit.reserve(gates);

        std::array<std::vector<double>, 3> times{};
        auto& [gate_times, measure_z_times, copy_times]{times};
        for (auto& samples : times) samples.reserve(repeats);
        constexpr std::array<std::string_view, 9> gate_names{"H", "X",  "S",  "Sdg", "Z",
                                                             "Y", "Cx", "Cz", "Cy"};
        std::array<std::vector<double>, gate_names.size()> individual_gate_times{};
        for (auto& samples : individual_gate_times) samples.reserve(repeats);
        std::size_t memory_bytes{0};
        [[maybe_unused]] volatile bool outcome{false};
        for (std::uint32_t repeat{0}; repeat < repeats; ++repeat) {
            std::mt19937 rng{1234U + repeat};
            circuit.clear();
            for (std::size_t i{0}; i < gates; ++i) {
                const auto q{qubit(rng)};
                auto target{qubit(rng)};
                if (target == q) target = (q + 1) % qubits;
                circuit.push_back(
                    {static_cast<GateKind>(rng() % (qubits == 1 ? 6 : 9)), q, target});
            }

            BenchSimulator simulator{qubits};
            const auto start{std::chrono::steady_clock::now()};
            for (const auto& gate : circuit) {
                switch (gate.kind) {
                    case GateKind::kH:
                        simulator.H(gate.q);
                        break;
                    case GateKind::kX:
                        simulator.X(gate.q);
                        break;
                    case GateKind::kS:
                        simulator.S(gate.q);
                        break;
                    case GateKind::kSdg:
                        simulator.Sdg(gate.q);
                        break;
                    case GateKind::kZ:
                        simulator.Z(gate.q);
                        break;
                    case GateKind::kY:
                        simulator.Y(gate.q);
                        break;
                    case GateKind::kCx:
                        simulator.Cx(gate.q, gate.target);
                        break;
                    case GateKind::kCz:
                        simulator.Cz(gate.q, gate.target);
                        break;
                    case GateKind::kCy:
                        simulator.Cy(gate.q, gate.target);
                        break;
                }
            }

            const auto stop{std::chrono::steady_clock::now()};
            gate_times.push_back(std::chrono::duration<double>{stop - start}.count());
            if (repeat == 0) memory_bytes = simulator.MemoryUsage();
            const auto copy_start{std::chrono::steady_clock::now()};
            auto measured{simulator};
            copy_times.push_back(
                std::chrono::duration<double>{std::chrono::steady_clock::now() - copy_start}
                    .count());
            const auto measurement_start{std::chrono::steady_clock::now()};
            for (stabaal::QubitSize q{0}; q < qubits; ++q) {
                outcome = measured.MeasureZ(q);
            }
            const auto measurement_stop{std::chrono::steady_clock::now()};
            measure_z_times.push_back(
                std::chrono::duration<double>{measurement_stop - measurement_start}.count());
            outcome = simulator.MeasureZ(0);

            const auto time_gate = [&](auto apply, std::vector<double>& samples) {
                auto individual{simulator};
                const auto gate_start{std::chrono::steady_clock::now()};
                for (const auto& gate : circuit) apply(individual, gate);
                const auto gate_stop{std::chrono::steady_clock::now()};
                samples.push_back(std::chrono::duration<double>{gate_stop - gate_start}.count());
                outcome = individual.MeasureZ(0);
            };
            time_gate([](auto& sim, const Gate& gate) { sim.H(gate.q); }, individual_gate_times[0]);
            time_gate([](auto& sim, const Gate& gate) { sim.X(gate.q); }, individual_gate_times[1]);
            time_gate([](auto& sim, const Gate& gate) { sim.S(gate.q); }, individual_gate_times[2]);
            time_gate([](auto& sim, const Gate& gate) { sim.S(gate.q); }, individual_gate_times[3]);
            time_gate([](auto& sim, const Gate& gate) { sim.Z(gate.q); }, individual_gate_times[4]);
            time_gate([](auto& sim, const Gate& gate) { sim.Y(gate.q); }, individual_gate_times[5]);
            if (qubits > 1) {
                time_gate([](auto& sim, const Gate& gate) { sim.Cx(gate.q, gate.target); },
                          individual_gate_times[6]);
                time_gate([](auto& sim, const Gate& gate) { sim.Cz(gate.q, gate.target); },
                          individual_gate_times[7]);
                time_gate([](auto& sim, const Gate& gate) { sim.Cy(gate.q, gate.target); },
                          individual_gate_times[8]);
            }
        }

        auto memory{static_cast<double>(memory_bytes)};
        constexpr std::array<std::string_view, 4> memory_units{"B", "KiB", "MiB", "GiB"};
        std::size_t memory_unit{0};
        while (memory >= 1024 && memory_unit + 1 < memory_units.size()) {
            memory /= 1024;
            ++memory_unit;
        }

        std::cout << "qubits=" << qubits << " gates=" << gates << " repeats=" << repeats
                  << std::fixed << std::setprecision(3) << "\nsimulator_memory=" << memory << ' '
                  << memory_units[memory_unit] << " (excluding allocator overhead)\n";
        PrintTiming("gates", gate_times, gates);
        for (std::size_t i{0}; i < gate_names.size(); ++i) {
            if (!individual_gate_times[i].empty()) {
                PrintTiming(gate_names[i], individual_gate_times[i], gates);
            }
        }

        PrintTiming("MeasureZ", measure_z_times, qubits);
        PrintTiming("copy", copy_times, 1);
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
