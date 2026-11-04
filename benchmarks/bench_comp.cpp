#include <algorithm>
#include <array>
#include <atomic>
#include <charconv>
#include <chrono>
#include <cstddef>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <exception>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <memory>
#include <mutex>
#include <random>
#include <span>
#include <sstream>
#include <stdexcept>
#include <string>
#include <string_view>
#include <system_error>
#include <thread>
#include <vector>
#include <utility>
#include <numeric>
#include <variant>

#include "stabaal.hpp"
#include "stim.h"

namespace {

using BenchSimulator = stabaal::Simulator<>;
using StimSimulator = stim::TableauSimulator<stim::MAX_BITWORD_WIDTH>;
using Distribution = std::uniform_int_distribution<stabaal::QubitSize>;

enum class GateKind : std::uint8_t { kI, kH, kX, kS, kSdg, kSx, kZ, kY, kCx, kCz, kCy, kMz, kMx, kMy };
constexpr auto kStimGates = std::to_array<stim::GateType>({
    stim::GateType::I,        // I
    stim::GateType::H,       // kH
    stim::GateType::X,       // kX
    stim::GateType::S,       // kS
    stim::GateType::S_DAG,   // kSdg
    stim::GateType::SQRT_X,  // kSx
    stim::GateType::Z,       // kZ
    stim::GateType::Y,       // kY
    stim::GateType::CX,      // kCx
    stim::GateType::CZ,      // kCz
    stim::GateType::CY,      // kCy
    stim::GateType::M,       // kMz
    stim::GateType::MX,      // kMx
    stim::GateType::MY,      // kMy
});

struct Gate {
    GateKind kind;
    stabaal::QubitSize q;
    stabaal::QubitSize target;
};

using Circuit = std::vector<Gate>;
using State = std::variant<BenchSimulator, StimSimulator>;

struct Result {
    double run_time;
    State state;
};

template <std::size_t W>
std::size_t simd_bits_heap(const stim::simd_bits<W>& bits) {
    return bits.num_simd_words * sizeof(stim::bitword<W>);
}

template <std::size_t W>
std::size_t stim_tableau_memory(const stim::Tableau<W>& t) {
    std::size_t bytes{sizeof(t)};
    for (const auto* half : {&t.xs, &t.zs}) {
        bytes += simd_bits_heap(half->xt.data) + simd_bits_heap(half->zt.data) + simd_bits_heap(half->signs);
    }
    return bytes;
}

std::string format_bytes(std::size_t bytes) {
    constexpr std::array<std::string_view, 4> units{"B", "KiB", "MiB", "GiB"};
    auto value{static_cast<double>(bytes)};
    std::size_t unit{0};
    while (value >= 1024 && unit + 1 < units.size()) {
        value /= 1024;
        ++unit;
    }
    std::ostringstream out;
    out << std::fixed << std::setprecision(2) << value << ' ' << units[unit] << " (" << bytes << " B)";
    return out.str();
}

}  // namespace

Result bench_stim(Circuit circuit, stabaal::QubitSize qubits) {
    StimSimulator simulator{std::mt19937_64{0}, qubits};
    const auto start{std::chrono::steady_clock::now()};
    for (const auto& gate : circuit) {
        const std::array targets{stim::GateTarget::qubit(gate.q), stim::GateTarget::qubit(gate.target)};
        const bool two_qubit{gate.kind >= GateKind::kCx && gate.kind <= GateKind::kCy};
        if (gate.kind >= GateKind::kMz) {
            simulator.sign_bias = (gate.target & 1) ? -1 : 1;
        }
        simulator.do_gate(stim::CircuitInstruction{
            kStimGates[static_cast<std::size_t>(gate.kind)],
            {},
            {targets.data(), targets.data() + (two_qubit ? 2 : 1)},
            ""});
    }

    const auto stop{std::chrono::steady_clock::now()};
    return {std::chrono::duration<double>{stop - start}.count(), simulator};
}

Result bench_stabaal(Circuit circuit, stabaal::QubitSize qubits) {
    BenchSimulator simulator{qubits};
    const auto start{std::chrono::steady_clock::now()};
    for (const auto& gate : circuit) {
        switch (gate.kind) {
            case GateKind::kI:
                simulator.I(gate.q);
                break;
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
            case GateKind::kSx:
                simulator.Sx(gate.q);
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
            case GateKind::kMz:
                (simulator.IsDeterministicMeasureZ(gate.q)) ?
                    simulator.MeasureZ(gate.q) :
                    simulator.ForceMeasureZ(gate.q, (gate.target & 1));
                break;
            case GateKind::kMx:
                (simulator.IsDeterministicMeasureX(gate.q)) ?
                    simulator.MeasureX(gate.q) :
                    simulator.ForceMeasureX(gate.q, (gate.target & 1));
                break;
            case GateKind::kMy:
                (simulator.IsDeterministicMeasureY(gate.q)) ?
                    simulator.MeasureY(gate.q) :
                    simulator.ForceMeasureY(gate.q, (gate.target & 1));
                break;
        }
    }

    const auto stop{std::chrono::steady_clock::now()};
    return {std::chrono::duration<double>{stop - start}.count(), simulator};
}

Circuit generate_circuit(size_t gates, std::mt19937 rng, stabaal::QubitSize qubits, Distribution qubit) {
    Circuit out{};
    for (std::size_t i{0}; i < gates; ++i) {
        const auto q{qubit(rng)};
        auto target{qubit(rng)};
        if (target == q) target = (q + 1) % qubits;
        out.emplace_back(static_cast<GateKind>(rng() % (qubits == 1 ? 8 : 11)), q, target);
    }
    return out;
}

// Writes how a single-qubit check came out to `why`, e.g. "X: stabaal random, stim +1".
void describe_check(std::ostream& why, char basis, bool stabaal_det, int stim_peek) {
    why << basis << ": stabaal " << (stabaal_det ? "deterministic" : "random") << ", stim "
        << (stim_peek == 0 ? "random" : (stim_peek > 0 ? "+1" : "-1"));
}

bool state_equivivalence(BenchSimulator stabaal, StimSimulator stim, stabaal::QubitSize qubits,
                         std::ostream& why) {
    std::vector<int> not_checked(qubits);
    std::iota(std::begin(not_checked), std::end(not_checked), 0);
    std::vector<int> collapsed;
    const auto fail = [&](int q) {
        if (!collapsed.empty()) {
            why << "after forcing MZ = 1 on qubit(s)";
            for (const auto c : collapsed) why << ' ' << c;
            why << ", ";
        }
        why << "qubit " << q << " differs:\n  ";
        describe_check(why, 'X', stabaal.IsDeterministicMeasureX(q), stim.peek_x(q));
        why << "\n  ";
        describe_check(why, 'Y', stabaal.IsDeterministicMeasureY(q), stim.peek_y(q));
        why << "\n  ";
        describe_check(why, 'Z', stabaal.IsDeterministicMeasureZ(q), stim.peek_z(q));
        why << '\n';
        return false;
    };

    while (!not_checked.empty()) {
        auto qubit = std::begin(not_checked);

        while (qubit != std::end(not_checked)) {
            const auto q{*qubit};
            auto det_x = stabaal.IsDeterministicMeasureX(q);
            auto det_y = stabaal.IsDeterministicMeasureY(q);
            auto det_z = stabaal.IsDeterministicMeasureZ(q);
            if (det_x != (stim.peek_x(q) != 0) ||
                det_y != (stim.peek_y(q) != 0) ||
                det_z != (stim.peek_z(q) != 0))
                return fail(q);

            // Deterministic measurements leave the state unchanged, so checking the sign is safe.
            if (det_x) {
                if (stabaal.MeasureX(q) != (stim.peek_x(q) == -1)) {
                    why << "stabaal measured X = " << !(stim.peek_x(q) == -1) << ", ";
                    return fail(q);
                }
                qubit = not_checked.erase(qubit);

            } else if (det_y) {
                if (stabaal.MeasureY(q) != (stim.peek_y(q) == -1)) {
                    why << "stabaal measured Y = " << !(stim.peek_y(q) == -1) << ", ";
                    return fail(q);
                }
                qubit = not_checked.erase(qubit);

            } else if (det_z) {
                if (stabaal.MeasureZ(q) != (stim.peek_z(q) == -1)) {
                    why << "stabaal measured Z = " << !(stim.peek_z(q) == -1) << ", ";
                    return fail(q);
                }
                qubit = not_checked.erase(qubit);

            } else {
                ++qubit;
            }
        }

        if (not_checked.empty()) break;

        auto rand_qubit = not_checked[not_checked.size()-1];
        not_checked.pop_back();
        collapsed.push_back(rand_qubit);

        stabaal.ForceMeasureZ(rand_qubit, true);
        stim.sign_bias = -1;
        const std::array targets{stim::GateTarget::qubit(rand_qubit)};
        stim.do_gate(stim::CircuitInstruction{
            kStimGates[static_cast<std::size_t>(GateKind::kMz)],
            {},
            {targets.data(), targets.data() + 1},
            ""
        });
    }
    return true;
}

bool is_measurement(GateKind kind) {
    return kind >= GateKind::kMz && kind <= GateKind::kMy;
}

bool is_two_qubit(GateKind kind) {
    return kind >= GateKind::kCx && kind <= GateKind::kCy;
}


std::vector<Gate> gate_options(stabaal::QubitSize qubits) {
    std::vector<Gate> options;
    for (stabaal::QubitSize q{0}; q < qubits; ++q) {
        for (auto kind : {GateKind::kH, GateKind::kX, GateKind::kS, GateKind::kSdg, GateKind::kSx,
                          GateKind::kZ, GateKind::kY}) {
            options.push_back({kind, q, 0});
        }
        for (stabaal::QubitSize t{0}; t < qubits; ++t) {
            if (t == q) continue;
            for (auto kind : {GateKind::kCx, GateKind::kCz, GateKind::kCy}) options.push_back({kind, q, t});
        }
        // for (auto kind : {GateKind::kMz, GateKind::kMx, GateKind::kMy}) {
        //     options.push_back({kind, q, 0});
        //     options.push_back({kind, q, 1});
        // }
    }
    return options;
}

// Lowest qubit the gate touches.
stabaal::QubitSize lowest_qubit(const Gate& gate) {
    return is_two_qubit(gate.kind) ? std::min(gate.q, gate.target) : gate.q;
}

bool shares_qubit(const Gate& a, const Gate& b) {
    const auto touches = [](const Gate& gate, stabaal::QubitSize q) {
        return gate.q == q || (is_two_qubit(gate.kind) && gate.target == q);
    };
    return touches(b, a.q) || (is_two_qubit(a.kind) && touches(b, a.target));
}

bool is_canonical(const Circuit& circuit) {
    for (std::size_t i{1}; i < circuit.size(); ++i) {
        const auto& a{circuit[i - 1]};
        const auto& b{circuit[i]};
        if (is_measurement(a.kind) || is_measurement(b.kind) || shares_qubit(a, b)) continue;
        if (lowest_qubit(a) > lowest_qubit(b)) return false;
    }
    return true;
}

// Steps the counter of option indices to the next circuit. Returns false after the last one.
bool next(std::vector<std::size_t>& digits, std::size_t num_options) {
    for (std::size_t i{digits.size()}; i > 0; --i) {
        if (++digits[i - 1] < num_options) return true;
        digits[i - 1] = 0;
    }
    return false;
}

constexpr auto kGateNames = std::to_array<std::string_view>(
    {"I", "H", "X", "S", "S_DAG", "SQRT_X", "Z", "Y", "CX", "CZ", "CY", "M", "MX", "MY"});
constexpr std::string_view kFailureLog{"equivalence_failures.log"};
constexpr std::string_view kRunLog{"equivalence_run.log"};
constexpr std::size_t kMaxFailures{100};
constexpr std::chrono::minutes kTimeout{5};

std::mutex output_mutex;
std::size_t failures{0};

// Stim-style circuit text, e.g. "H 0\nCX 0 1\nMX 1 (forced 1)".
void print_circuit(std::ostream& out, const Circuit& circuit) {
    for (const auto& gate : circuit) {
        out << "  " << kGateNames[static_cast<std::size_t>(gate.kind)] << ' ' << gate.q;
        if (is_two_qubit(gate.kind)) out << ' ' << gate.target;
        if (is_measurement(gate.kind)) out << " (forced " << (gate.target & 1) << ')';
        out << '\n';
    }
}

// Checks one circuit. A mismatch is appended to kFailureLog and the run keeps going.
// Returns false once kMaxFailures mismatches have been logged, to stop the run.
bool do_job(const Circuit& circuit, stabaal::QubitSize qubits) {
    auto stabaal_result = bench_stabaal(circuit, qubits);
    auto stim_result = bench_stim(circuit, qubits);

    auto stabaal_state = std::get<BenchSimulator>(stabaal_result.state);
    auto stim_state = std::get<StimSimulator>(stim_result.state);

    std::ostringstream why;
    if (state_equivivalence(stabaal_state, stim_state, qubits, why)) return true;

    // Only one thread writes at a time, so reports don't get mixed together.
    const std::scoped_lock lock{output_mutex};
    if (failures >= kMaxFailures) return false;
    ++failures;

    std::ofstream log{std::string{kFailureLog}, std::ios::app};
    log << "=== Mismatch on " << qubits << " qubit(s) ===\ncircuit:\n";
    print_circuit(log, circuit);
    log << why.str() << "stabaal tableau:";
    {
        // Print() writes to std::cout, so point std::cout at the log while it runs.
        auto* const old{std::cout.rdbuf(log.rdbuf())};
        stabaal_state.Print();
        std::cout.rdbuf(old);
    }
    log << "stim stabilizers:\n";
    for (const auto& s : stim_state.canonical_stabilizers()) log << "  " << s << '\n';
    log << '\n';

    std::cout << "Mismatch " << failures << " written to " << kFailureLog << '\n';
    return failures < kMaxFailures;
}

// hardware_concurrency() counts the whole node, not the CPUs Slurm allocated to the job.
unsigned worker_count() {
    if (const char* const cpus{std::getenv("SLURM_CPUS_PER_TASK")}) {
        return static_cast<unsigned>(std::max(1, std::atoi(cpus)));
    }
    return std::max(1U, std::thread::hardware_concurrency());
}

bool exhaustive_equivivalence(stabaal::QubitSize start_qubits, stabaal::QubitSize end_qubits,
                              std::size_t start_length, std::size_t end_length,
                              std::chrono::minutes timeout) {
    const auto deadline{std::chrono::steady_clock::now() + timeout};
    constexpr std::size_t kBatchSize{256};
    const unsigned workers{worker_count()};
    std::atomic<bool> stop{false};

    for (stabaal::QubitSize qubits{start_qubits}; qubits <= end_qubits && !stop; ++qubits) {
        const auto options{gate_options(qubits)};
        for (std::size_t length{start_length}; length <= end_length && !stop; ++length) {
            const auto start{std::chrono::steady_clock::now()};
            std::vector<std::size_t> digits(length, 0);
            bool done{false};
            std::atomic<std::size_t> tested{0};
            std::size_t skipped{0};
            std::mutex generator_mutex;

            auto worker = [&] {
                std::vector<Circuit> batch;
                Circuit circuit(length);
                while (!stop) {
                    batch.clear();
                    {
                        const std::scoped_lock lock{generator_mutex};
                        while (!done && batch.size() < kBatchSize) {
                            // Checked per circuit: long runs of skipped circuits may never fill a batch.
                            if (std::chrono::steady_clock::now() >= deadline) {
                                stop = true;
                                break;
                            }
                            for (std::size_t i{0}; i < length; ++i) circuit[i] = options[digits[i]];
                            if (is_canonical(circuit)) batch.push_back(circuit);
                            else ++skipped;
                            done = !next(digits, options.size());
                        }
                    }
                    if (batch.empty()) return;

                    for (const auto& job : batch) {
                        if (!do_job(job, qubits)) {
                            stop = true;
                            return;
                        }
                    }
                    tested += batch.size();
                }
            };

            {
                // jthreads join automatically when pool goes out of scope.
                std::vector<std::jthread> pool;
                for (unsigned i{0}; i < workers; ++i) pool.emplace_back(worker);
            }

            const std::chrono::duration<double> elapsed{std::chrono::steady_clock::now() - start};
            std::cout << qubits << " qubits, " << length << " gates: " << tested << " tested, "
                      << skipped << " skipped (" << elapsed.count() << " s)\n";
        }
    }
    const std::scoped_lock lock{output_mutex};
    return failures == 0;
}

int main() {
    std::remove(std::string{kFailureLog}.c_str());
    std::ofstream run_log{std::string{kRunLog}};
    auto* const terminal{std::cout.rdbuf(run_log.rdbuf())};
    std::cout << std::unitbuf;  // Flush every write, since a run usually ends by being killed.
    bool error{false};
    size_t i = 1; 
    do {
        error = !exhaustive_equivivalence(4, 4, i, i, kTimeout);
        ++i;
    } while (!error);

    std::cout << failures << " mismatch(es), see " << kFailureLog << '\n';
    std::cout.rdbuf(terminal);  // run_log dies before std::cout, which must not keep using it.
    return EXIT_FAILURE;
    
    //  try {
    //     const stabaal::QubitSize qubits{1U << 8};
    //     const auto gates{100000};
    //     const auto repeats{100};
    //     Distribution qubit{0, qubits - 1};
    //     Circuit circuit{};
    //     circuit.reserve(gates);
    //     [[maybe_unused]] volatile bool outcome{false};

    //     for (std::uint32_t repeat{0}; repeat < repeats; ++repeat) {
    //         std::mt19937 rng{1234U + repeat};
    //         circuit.clear();     
    //         circuit = generate_circuit(gates, rng, qubits, qubit);
    //         auto stabaal_result = bench_stabaal(circuit, qubits);
    //         auto stim_result = bench_stim(circuit, qubits);
    //         std::cout << "stabaal time: " << stabaal_result.run_time << "\n";
    //         std::cout << "stim    time: " << stim_result.run_time    << "\n\n";

    //         if (repeat == 0) {
    //             const auto& stim_sim{std::get<StimSimulator>(stim_result.state)};
    //             std::cout << "stabaal tableau mem: "
    //                       << format_bytes(std::get<BenchSimulator>(stabaal_result.state).MemoryUsage()) << "\n";
    //             std::cout << "stim    tableau mem: " << format_bytes(stim_tableau_memory(stim_sim.inv_state)) << "\n";
    //             std::cout << "stim    record  mem: "
    //                       << format_bytes(stim_sim.measurement_record.storage.capacity() / 8) << "\n\n";
    //         }
    //     }

    // } catch (const std::exception& error) {
    //     std::cerr << error.what() << '\n';
    //     return 1;
    // }
}
