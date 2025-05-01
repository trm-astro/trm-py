#include <pybind11/pybind11.h>
#include "py_subs.h"
#include "trm/subs.h"

// Define the alias "py" for the namespace "pybind11"
namespace py = pybind11;

void init_subs(py::module &);

PYBIND11_MODULE(_cpp_subs, m) {
    m.doc() = "Subs Library";

    // Create the submodule `subs`
    //py::module_ cpp_subs = m.def_submodule("cpp_subs", "interface to cpp-subs");
    init_subs(m);
}

void init_subs(py::module_ &m) {
    // Include any further python bindings here
    // Overloads should be made in the C++ and then exposed to Python here

    // Overload resolution for the same name 'voigt'
    m.def("voigt", 
          [](double a, double v, double eps = 1.e-8) {
              return Subs::voigt(a, v, eps);
          }, 
          "The Voigt function from cpp-subs",
          py::arg("a"), py::arg("v"), py::arg("eps") = 1.e-8
    );
    m.def("voigt", 
            [](double x, const std::vector<double>& sigma, double eps) { // Lambda for vectorized version
                double* result; // To store results
                Subs::voigt(x, sigma.data(), result, sigma.size(), eps); // Call the voigt function
                return result; // Return the result
            }, 
            "The Voigt function from cpp-subs"
    );  // Documentation for array input
    
    // Overload resolution for the same name 'gammaq'
    m.def("gammq",
        [](double a, double x) {
            return Subs::gammq(a, x);
        },
        "The incomplete gamma function from cpp-subs",
        py::arg("a"), py::arg("x")
    );
    m.def("gammq", 
        [](double a, std::vector<double>& x) { // Lambda for vectorized version
            double* result; // To store results
            Subs::gammq(a, x.data(), result, x.size()); // Call the vectorized gammaq
            return result; // Return the result
        }, 
        "The incomplete gamma function from cpp-subs");
}