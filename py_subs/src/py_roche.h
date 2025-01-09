
#ifndef PY_ROCHE_H
#define PY_ROCHE_H

#include <pybind11/pybind11.h>

// Define the alias "py" for the namespace "pybind11"
namespace py = pybind11;

void init_roche(py::module_ &m);

#endif // PY_ROCHE_H