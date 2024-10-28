#include <pybind11/pybind11.h>
#include "cpp-subs/include/trm/subs.h"

// Define the alias "py" for the namespace "pybind11"
namespace py = pybind11;

void init_subs(py::module_ &m) {
    // Include any further python bindings here
    // Overloads should be made in the C++ and then exposed to Python here

    // Overload resolution for the same name 'voigt'
    m.def("voigt", &voigt, "The Voigt function from cpp-subs",
        py::arg("eps") = 1.e-8);
    m.def("voigt", 
            [](double x, const std::vector<double>& sigma, double eps) { // Lambda for vectorized version
                std::vector<double> result; // To store results
                voigt(x, sigma.data(), result, sigma.size(), eps); // Call the voigt function
                return result; // Return the result
            }, 
            "The Voigt function from cpp-subs");  // Documentation for array input
    
    // Overload resolution for the same name 'gammaq'
    m.def("gammaq", &gammaq, "The incomplete gamma function from cpp-subs");
    m.def("gammaq", 
            [](double a, const std::vector<double>& x) { // Lambda for vectorized version
                std::vector<double> result; // To store results
                gammaq(a, x.data(), result, x.size()); // Call the vectorized gammaq
                return result; // Return the result
            }, 
            "The incomplete gamma function from cpp-subs");  // Documentation for array input
}