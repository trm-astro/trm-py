#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>
#include "py_lcurve.h"
#include "trm/lcurve.h"

namespace py = pybind11;

PYBIND11_MODULE(_cpp_lcurve, m) {
    m.doc() = "Lcurve Library";

    // Initialize the submodule `lcurve`
    init_lcurve(m);
}

void init_lcurve(py::module_ &m) {
    // Bind the Lcurve::run_lroche function
    m.def("run_lroche",
          &Lcurve::run_lroche,
          "Run the lroche function with model and data files",
          py::arg("model_file"), py::arg("data_file"));

    // Add other bindings as needed
    // For example, you can bind the Lcurve::Data class, Fobj class, etc.
    // m.def("some_other_function", &Lcurve::some_other_function);
}