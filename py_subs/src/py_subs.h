#ifndef PY_SUBS_H
#define PY_SUBS_H

#include <pybind11/pybind11.h>

// Define the alias "py" for the namespace "pybind11"
namespace py = pybind11;

void init_subs(py::module_ &m);

#endif // PY_SUBS_H