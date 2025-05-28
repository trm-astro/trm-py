#ifndef PY_LCURVE_H
#define PY_LCURVE_H

#include <pybind11/pybind11.h>

// Define the alias "py" for the namespace "pybind11"
namespace py = pybind11;

void init_lcurve(py::module_ &m);

#endif // PY_LCURVE_H