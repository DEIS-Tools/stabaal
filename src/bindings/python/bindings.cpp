#include <pybind11/cast.h>
#include <pybind11/iostream.h>
#include <pybind11/pybind11.h>
#include <pybind11/pytypes.h>

#include "stabaal/simulator.hpp"
#include "stabaal/types.hpp"

namespace py = pybind11;

using InverseSim = stabaal::Simulator<>;

PYBIND11_MODULE(stabaal, m) {
    m.doc() = "Dynamic Stabilizer Circuit Simulation";

    auto inverse =
        py::class_<InverseSim>(m, "InverseTableau",
                               "Stabilizer simulator using an inverse tableau backend.")
            .def(py::init<stabaal::QubitSize>(), py::arg("num_qubits"),
                 "Create a simulator with all qubits initialized to |0>.")
            .def("add", &InverseSim::Add, "Append a qubit initialized to |0>.")
            .def("num_qubits", &InverseSim::NumQubits, "Return the number of active qubits.")
            .def("discard", &InverseSim::Discard, py::arg("q"), "Discard qubit q from the tableau.")
            .def("h", &InverseSim::H, py::arg("q"), "Apply a Hadamard gate to qubit q.")
            .def("s", &InverseSim::S, py::arg("q"), "Apply an S gate to qubit q.")
            .def("sdg", &InverseSim::Sdg, py::arg("q"), "Apply an S dagger gate to qubit q.")
            .def("i", &InverseSim::I, py::arg("q"), "Apply the identity gate to qubit q.")
            .def("x", &InverseSim::X, py::arg("q"), "Apply a X gate to qubit q.")
            .def("z", &InverseSim::Z, py::arg("q"), "Apply a Z gate to qubit q.")
            .def("y", &InverseSim::Y, py::arg("q"), "Apply a Y gate to qubit q.")
            .def("cx", &InverseSim::Cx, py::arg("control"), py::arg("target"),
                 "Apply a CX gate from control to target.")
            .def("cz", &InverseSim::Cz, py::arg("control"), py::arg("target"),
                 "Apply a CZ gate from control to target.")
            .def("cy", &InverseSim::Cy, py::arg("control"), py::arg("target"),
                 "Apply a CY gate from control to target.")
            .def("measure_x", &InverseSim::MeasureX, py::arg("q"),
                 "Measure qubit q in Hadamard basis.")
            .def("measure_z", &InverseSim::MeasureZ, py::arg("q"),
                 "Measure qubit q in computational basis.")
            .def("is_deterministic_measure_x", &InverseSim::IsDeterministicMeasureX,
                 py::arg("q"), "Return whether an X measurement is deterministic.")
            .def("is_deterministic_measure_z", &InverseSim::IsDeterministicMeasureZ,
                 py::arg("q"), "Return whether a Z measurement is deterministic.")
            .def("force_measure_x", &InverseSim::ForceMeasureX, py::arg("q"),
                 py::arg("measurement"), "Force a non-deterministic X measurement outcome.")
            .def("force_measure_z", &InverseSim::ForceMeasureZ, py::arg("q"),
                 py::arg("measurement"), "Force a non-deterministic Z measurement outcome.")
            .def("print", &InverseSim::Print, "Print the backend's tableau to stdout.")
            .def(
                "__repr__",
                [](const InverseSim& sim) {
                    auto output = py::module_::import("io").attr("StringIO")();
                    {
                        py::scoped_ostream_redirect redirect(std::cout, output);
                        sim.Print();
                    }
                    return output.attr("getvalue")().cast<std::string>();
                },
                "Print the backend's tableau to stdout")
            .def("memory_usage", &InverseSim::MemoryUsage,
                 "Return the simulator's memory usage in bytes.");

    m.def(
        "Simulator", [](stabaal::QubitSize num_qubits) { return InverseSim(num_qubits); },
        py::arg("num_qubits"), "Create a simulator using the default inverse tableau backend.");

    m.def(
        "Simulator",
        [](stabaal::QubitSize num_qubits, const py::object& backend) {
            return backend(num_qubits);
        },
        py::arg("num_qubits"), py::arg_v("backend", inverse, "InverseTableau"),
        "Create a simulator using the selected backend constructor. Defaults to InverseTableau.");
}
