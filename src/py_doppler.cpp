#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include "doppler.h"
#include "py_doppler.h"

namespace py = pybind11;

void init_doppler(py::module &);

PYBIND11_MODULE(_cpp_doppler, m) {
    m.doc() = "Doppler Library";
    // Initialize the submodule `doppler`
    init_doppler(m);
}

void init_doppler(py::module_ &m) {
    m.def("comdat", &doppler_comdat, py::arg("map"), py::arg("data"),
          "Computes the data equivalent to a Doppler map using a Data object as the template.");

    m.def("comdef", &doppler_comdef, py::arg("map"),
          "Computes the default image for a Doppler map.");

    m.def("datcom", &doppler_datcom, py::arg("data"), py::arg("map"),
          "Computes the transpose operation to comdat for testing.");

    m.def("memit", &doppler_memit, py::arg("map"), py::arg("data"), py::arg("niter"), py::arg("caim"),
          py::arg("tlim") = 1.e-4, py::arg("rmax") = 0.2,
          "Carries out MEM iterations on a map given data.");
}

