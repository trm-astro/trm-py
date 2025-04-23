#ifndef PY_DOPPLER_H
#define PY_DOPPLER_H

#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>

// Define the alias "py" for the namespace "pybind11"
namespace py = pybind11;

void init_doppler(py::module_ &m);


#endif // PY_DOPPLER_H