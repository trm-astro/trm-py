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
    m.def("lroche",
          &Lcurve::run_lroche,
          "Run the lroche function with model and data files",
          py::arg("model_file"), py::arg("data_file"));
    
    m.def("lprofile",
          &Lcurve::run_lprofile,
          "Run the lprofile function with a model file",
          py::arg("model_file"));
    
    m.def("levmarq",
          &Lcurve::run_levmarq,
          "Run the levmarq function with model and data files",
          py::arg("model_file"), py::arg("data_file"));

}