#include <pybind11/cast.h>
#include <pybind11/pybind11.h>
#include <pybind11/pytypes.h>

#include "stabaal/simulator.hpp"

namespace py = pybind11;

using InverseSim = stabaal::Simulator<>;

PYBIND11_MODULE(stabaal, m) {
    m.doc() = "Dynamic Stabilizer Circuit Simulation";

    auto inverse = py::class_<InverseSim>(m, "InverseTableau",
                                          "Stabilizer simulator using an inverse tableau backend.")
                       .def(py::init<stabaal::QubitIdx>(), py::arg("num_qubits"),
                            "Create a simulator with all qubits initialized to |0>.")
                       .def("add", &InverseSim::add, "Append a qubit initialized to |0>.")
                       .def("h", &InverseSim::h, py::arg("q"), "Apply a Hadamard gate to qubit q.")
                       .def("s", &InverseSim::s, py::arg("q"), "Apply an S gate to qubit q.")
                       .def("x", &InverseSim::x, py::arg("q"), "Apply a X gate to qubit q.")
                       .def("z", &InverseSim::z, py::arg("q"), "Apply a Z gate to qubit q.")
                       .def("cx", &InverseSim::cx, py::arg("control"), py::arg("target"),
                            "Apply a CX gate from control to target.")
                       .def("measure_x", &InverseSim::measureX, py::arg("q"),
                            "Measure qubit q in Hadamard basis.")
                       .def("measure_z", &InverseSim::measureZ, py::arg("q"),
                            "Measure qubit q in computational basis.")
                       .def("print", &InverseSim::print, "Print the backend's tableau to stdout.")
                       .def("memory_usage", &InverseSim::memoryUsage,
                            "Return the simulator's memory usage in bytes.");

    m.def(
        "Simulator",
        [](stabaal::QubitIdx num_qubits, const py::object& backend) { return backend(num_qubits); },
        py::arg("num_qubits"), py::arg("backend") = inverse,
        "Create a simulator using the selected backend constructor. Defaults to InverseTableau.");
}
