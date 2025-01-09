#include <pybind11/pybind11.h>
#include "py_subs.h"
#include "py_roche.h"
#include "trm/subs.h"
#include "trm/roche.h"


// Define the alias "py" for the namespace "pybind11"
namespace py = pybind11;

void init_roche(py::module_ &m) {
    // Functions from croche.pxd

    // The STAR enum from roche.h
    py::enum_<Roche::STAR>(m, "STAR")
        .value("PRIMARY", Roche::STAR::PRIMARY)
        .value("SECONDARY", Roche::STAR::SECONDARY)
        .export_values();

    m.def("face",
            [](double q, Roche::STAR star, double spin, 
                const Subs::Vec3& dirn, double rref, double pref, double acc) {
                Subs::Vec3 pvec, dvec;
                double r, g;
                Roche::face(q, star, spin, dirn, rref, pref, acc, pvec, dvec, r, g);
                return std::make_tuple(pvec, dvec, r, g);
            },
            "Computes the position of a point on the Roche distorted surface",
            py::arg("q"), py::arg("star"), py::arg("spin"), py::arg("dirn"), py::arg("rref"), py::arg("pref"), py::arg("acc")
    );

}