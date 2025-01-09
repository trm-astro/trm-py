#include <pybind11/pybind11.h>
#include "py_subs.h"
#include "py_roche.h"

// Define the alias "py" for the namespace "pybind11"
namespace py = pybind11;

void init_subs(py::module &);
// void init_binary(py::module &);
// void init_colly(py::module &);
// void init_lcurve(py::module &);
void init_roche(py::module &);


PYBIND11_MODULE(py_subs, m) {
    m.doc() = "A package with multiple features (submodules)";

    // Create the submodule `subs`
    py::module_ cpp_subs = m.def_submodule("cpp_subs", "interface to cpp-subs");
    init_subs(cpp_subs);

    // // Initialize the submodule `binary`
    // py::module_ binary = m.def_submodule("binary", "interface to cpp-binary");
    // init_binary(m);

    // // Initialize the submodule `colly`
    // py::module_ colly = m.def_submodule("colly", "interface to cpp-colly");
    // init_colly(m);

    // // Initialize the submodule `lcurve`
    // py::module_ lcurve = m.def_submodule("lcurve", "interface to cpp-lcurve");
    // init_lcurve(m);

    // Initialize the submodule `roche`
    py::module_ roche = m.def_submodule("roche", "interface to cpp-roche");
    init_roche(roche);

}