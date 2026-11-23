#include <pybind11/cast.h>
#include <pybind11/iostream.h>
#include <pybind11/pybind11.h>
#include <pybind11/pytypes.h>
#include <pybind11/stl.h>

#include "stabaal/backends/pauli_noise_backend.hpp"
#include "stabaal/noise/pauli_noise.hpp"
#include "stabaal/simulator.hpp"
#include "stabaal/types.hpp"

namespace py = pybind11;

using InverseSim = stabaal::Simulator<>;

template <typename Sim>
auto BindSimulator(py::module_& m, const char* name, const char* description) {
    auto simulator =
        py::class_<Sim>(m, name, description)
            .def(py::init<stabaal::QubitSize>(), py::arg("num_qubits"),
                 "Create a simulator with all qubits initialized to |0>.")
            .def("add", &Sim::Add, "Allocate a qubit initialized to |0> and return its id.")
            .def("num_qubits", &Sim::NumQubits, "Return the number of active qubits.")
            .def("discard", &Sim::Discard, py::arg("q"), "Discard qubit q.")
            .def("h", py::overload_cast<stabaal::QubitId>(&Sim::H), py::arg("q"),
                 "Apply a Hadamard gate to qubit q.")
            .def("s", py::overload_cast<stabaal::QubitId>(&Sim::S), py::arg("q"),
                 "Apply an S gate to qubit q.")
            .def("sdg", py::overload_cast<stabaal::QubitId>(&Sim::Sdg), py::arg("q"),
                 "Apply an S† gate to qubit q.")
            .def("sx", py::overload_cast<stabaal::QubitId>(&Sim::Sx), py::arg("q"),
                 "Apply an √X gate to qubit q.")
            .def("i", py::overload_cast<stabaal::QubitId>(&Sim::I), py::arg("q"),
                 "Apply the identity gate to qubit q.")
            .def("x", py::overload_cast<stabaal::QubitId>(&Sim::X), py::arg("q"),
                 "Apply a X gate to qubit q.")
            .def("z", py::overload_cast<stabaal::QubitId>(&Sim::Z), py::arg("q"),
                 "Apply a Z gate to qubit q.")
            .def("y", py::overload_cast<stabaal::QubitId>(&Sim::Y), py::arg("q"),
                 "Apply a Y gate to qubit q.")
            .def("cx", py::overload_cast<stabaal::QubitId, stabaal::QubitId>(&Sim::Cx),
                 py::arg("control"), py::arg("target"), "Apply a CX gate from control to target.")
            .def("cz", py::overload_cast<stabaal::QubitId, stabaal::QubitId>(&Sim::Cz),
                 py::arg("control"), py::arg("target"), "Apply a CZ gate from control to target.")
            .def("cy", py::overload_cast<stabaal::QubitId, stabaal::QubitId>(&Sim::Cy),
                 py::arg("control"), py::arg("target"), "Apply a CY gate from control to target.")
            .def("measure_x", &Sim::MeasureX, py::arg("q"), "Measure qubit q in Hadamard basis.")
            .def("measure_z", &Sim::MeasureZ, py::arg("q"),
                 "Measure qubit q in computational basis.")
            .def("measure_and_discard_x", &Sim::MeasureAndDiscardX, py::arg("q"),
                 "Measure qubit q in the X basis, discard it, and return the boolean outcome.")
            .def("measure_and_discard_z", &Sim::MeasureAndDiscardZ, py::arg("q"),
                 "Measure qubit q in the Z basis, discard it, and return the boolean outcome.")
            .def("is_deterministic_measure_x", &Sim::IsDeterministicMeasureX, py::arg("q"),
                 "Return whether an X measurement is deterministic.")
            .def("is_deterministic_measure_z", &Sim::IsDeterministicMeasureZ, py::arg("q"),
                 "Return whether a Z measurement is deterministic.")
            .def("force_measure_x", &Sim::ForceMeasureX, py::arg("q"), py::arg("measurement"),
                 "Force a non-deterministic X measurement outcome.")
            .def("force_measure_z", &Sim::ForceMeasureZ, py::arg("q"), py::arg("measurement"),
                 "Force a non-deterministic Z measurement outcome.")
            .def("print", &Sim::Print, "Print the backend's tableau to stdout.")
            .def(
                "__repr__",
                [](const Sim& sim) {
                    auto output = py::module_::import("io").attr("StringIO")();
                    {
                        py::scoped_ostream_redirect redirect(std::cout, output);
                        sim.Print();
                    }
                    return output.attr("getvalue")().cast<std::string>();
                },
                "Print the backend's tableau to stdout")
            .def("memory_usage", &Sim::MemoryUsage,
                 "Return the simulator's memory usage in bytes.");

    if constexpr (std::constructible_from<Sim, stabaal::QubitSize, std::optional<std::uint64_t>,
                                          std::optional<std::uint64_t>>) {
        simulator.def(py::init<stabaal::QubitSize, std::optional<std::uint64_t>,
                               std::optional<std::uint64_t>>(),
                      py::arg("num_qubits"), py::arg("seed") = py::none(),
                      py::arg("noise_seed") = py::none(),
                      "Create a simulator with all qubits initialized to |0>. "
                      "Use seed for measurement randomness and noise_seed for independent noise "
                      "randomness. Each seed defaults to None for automatic seeding.");
        simulator.def(
            "h",
            [](Sim& sim, stabaal::QubitId q, const stabaal::PauliNoise& noise) { sim.H(q, noise); },
            py::arg("q"), py::arg("noise"),
            "Apply a Hadamard gate to qubit q, then sample and apply a Pauli error.");
        simulator.def(
            "s",
            [](Sim& sim, stabaal::QubitId q, const stabaal::PauliNoise& noise) { sim.S(q, noise); },
            py::arg("q"), py::arg("noise"),
            "Apply an S gate to qubit q, then sample and apply a Pauli error.");
        simulator.def(
            "sdg",
            [](Sim& sim, stabaal::QubitId q, const stabaal::PauliNoise& noise) {
                sim.Sdg(q, noise);
            },
            py::arg("q"), py::arg("noise"),
            "Apply an S† gate to qubit q, then sample and apply a Pauli error.");
        simulator.def(
            "sx",
            [](Sim& sim, stabaal::QubitId q, const stabaal::PauliNoise& noise) {
                sim.Sx(q, noise);
            },
            py::arg("q"), py::arg("noise"),
            "Apply a √X gate to qubit q, then sample and apply a Pauli error.");
        simulator.def(
            "i",
            [](Sim& sim, stabaal::QubitId q, const stabaal::PauliNoise& noise) { sim.I(q, noise); },
            py::arg("q"), py::arg("noise"),
            "Apply the identity gate to qubit q, then sample and apply a Pauli error.");
        simulator.def(
            "x",
            [](Sim& sim, stabaal::QubitId q, const stabaal::PauliNoise& noise) { sim.X(q, noise); },
            py::arg("q"), py::arg("noise"),
            "Apply a X gate to qubit q, then sample and apply a Pauli error.");
        simulator.def(
            "z",
            [](Sim& sim, stabaal::QubitId q, const stabaal::PauliNoise& noise) { sim.Z(q, noise); },
            py::arg("q"), py::arg("noise"),
            "Apply a Z gate to qubit q, then sample and apply a Pauli error.");
        simulator.def(
            "y",
            [](Sim& sim, stabaal::QubitId q, const stabaal::PauliNoise& noise) { sim.Y(q, noise); },
            py::arg("q"), py::arg("noise"),
            "Apply a Y gate to qubit q, then sample and apply a Pauli error.");
        simulator.def(
            "cx",
            [](Sim& sim, stabaal::QubitId control, stabaal::QubitId target,
               const stabaal::PauliNoise& control_noise, const stabaal::PauliNoise& target_noise) {
                sim.Cx(control, target, control_noise, target_noise);
            },
            py::arg("control"), py::arg("target"), py::arg("control_noise"),
            py::arg("target_noise"),
            "Apply a CX gate from control to target, then independently sample and apply Pauli "
            "errors to control and target, in that order.");
        simulator.def(
            "cz",
            [](Sim& sim, stabaal::QubitId control, stabaal::QubitId target,
               const stabaal::PauliNoise& control_noise, const stabaal::PauliNoise& target_noise) {
                sim.Cz(control, target, control_noise, target_noise);
            },
            py::arg("control"), py::arg("target"), py::arg("control_noise"),
            py::arg("target_noise"),
            "Apply a CZ gate from control to target, then independently sample and apply Pauli "
            "errors to control and target, in that order.");
        simulator.def(
            "cy",
            [](Sim& sim, stabaal::QubitId control, stabaal::QubitId target,
               const stabaal::PauliNoise& control_noise, const stabaal::PauliNoise& target_noise) {
                sim.Cy(control, target, control_noise, target_noise);
            },
            py::arg("control"), py::arg("target"), py::arg("control_noise"),
            py::arg("target_noise"),
            "Apply a CY gate from control to target, then independently sample and apply Pauli "
            "errors to control and target, in that order.");
    } else {
        simulator.def(py::init<stabaal::QubitSize, std::optional<std::uint64_t>>(),
                      py::arg("num_qubits"), py::arg("seed") = py::none(),
                      "Create a simulator with all qubits initialized to |0>. "
                      "Use seed for measurement randomness, or None for automatic seeding.");
    }

    return simulator;
}

PYBIND11_MODULE(stabaal, m) {
    m.doc() = "Dynamic Stabilizer Circuit Simulation";
    py::class_<stabaal::PauliNoise>(m, "PauliNoise")
        .def(py::init<double, double, double, double>(), py::arg("wi"), py::arg("wx"),
             py::arg("wy"), py::arg("wz"),
             "Create a Pauli noise channel from relative I, X, Y, Z weights, normalized to "
             "probabilities. Weights must be finite, nonnegative, and not all zero.");
    auto inverse = BindSimulator<InverseSim>(
        m, "InverseTableau", "Stabilizer simulator using an inverse tableau backend.");
    BindSimulator<stabaal::Simulator<stabaal::PauliNoiseBackend<stabaal::InverseTableau>>>(
        m, "NoisyInverseTableau",
        "Stabilizer simulator using an inverse tableau backend with per-gate Pauli noise.");

    m.def(
        "Simulator",
        [](stabaal::QubitSize num_qubits, const py::object& backend,
           std::optional<std::uint64_t> seed, std::optional<std::uint64_t> noise_seed) {
            if (noise_seed)
                return backend(num_qubits, py::arg("seed") = seed,
                               py::arg("noise_seed") = noise_seed);
            if (seed) return backend(num_qubits, py::arg("seed") = seed);
            return backend(num_qubits);
        },
        py::arg("num_qubits"), py::arg_v("backend", inverse, "InverseTableau"),
        py::arg("seed") = py::none(), py::arg("noise_seed") = py::none(),
        "Create a simulator using the selected backend constructor. Defaults to InverseTableau. "
        "Use seed for measurement randomness and noise_seed for independent noise randomness "
        "on noise-capable backends. Each seed defaults to None for automatic seeding.");
}
