#include <pybind11/pybind11.h>
#include "py_subs.h"
#include "py_roche.h"
#include "trm/subs.h"
#include "trm/roche.h"


// Define the alias "py" for the namespace "pybind11"
namespace py = pybind11;

void init_roche(py::module_ &m) {
    // Functions from croche.pxd

    // TODO: Add a converter to convert Subs::Vec3 to a tuple/list and vice versa
    py::class_<Subs::Vec3>(m, "Vec3")
        .def(py::init<>())
        .def(py::init<double, double, double>())
        .def("x", (const double& (Vec3::*)() const) &Vec3::x)
        .def("y", (const double& (Vec3::*)() const) &Vec3::y)
        .def("z", (const double& (Vec3::*)() const) &Vec3::z)
        .def("set", (void (Vec3::*)(double, double, double)) &Vec3::set)
        .def("set", (void (Vec3::*)(double*)) &Vec3::set)
        .def("get", &Vec3::get)
        .def("__repr__", [](const Vec3 &v) {
            return "<Vec3 x=" + std::to_string(v.x()) + 
                   " y=" + std::to_string(v.y()) + 
                   " z=" + std::to_string(v.z()) + ">";
        })
        .def(py::init([](py::list l) {
            if (l.size() != 3) {
                throw std::runtime_error("List must have exactly 3 elements");
            }
            return new Vec3(l[0].cast<double>(), l[1].cast<double>(), l[2].cast<double>());
        }))
        .def(py::init([](py::tuple t) {
            if (t.size() != 3) {
                throw std::runtime_error("Tuple must have exactly 3 elements");
            }
            return new Vec3(t[0].cast<double>(), t[1].cast<double>(), t[2].cast<double>());
        }));

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
            "p,d,r,g = face(q, spin, dirn, rref, pref, star=2, acc=1.e-5), returns position and direction of element of specific Roche potential.
            q      -- mass ratio = M2/M1
            spin   -- ratio spin/orbital frequency
            dirn   -- direction (a Vec3) to take from centre of mass of star in question.
            rref   -- reference radius greater than any radius of potential in question.
            pref   -- the potential to aim for.
            star   -- 1 or 2 for primary or secondary star. acc    -- accuracy in terms of separation of location.

            Returns p = position, d = direction perpendicular to face, r = radius from centre of mass, g = gravity.
            ",
            py::arg("q"), py::arg("star"), py::arg("spin"), py::arg("dirn"), py::arg("rref"), py::arg("pref"), py::arg("acc")
    );

    m.def("bspot",
            [](double q, double rad, double acc = 1.e-7) {
                Subs::Vec3 r, v, rs, vs;
                Roche::strinit(q, r, v);
                try{
                    Roche::stradv(q, r, v, rad, acc, 1.e-2);
                    return std::make_tuple(r.x(), r.y(), v.x(), v.y());
                }
                catch(const Roche::Roche_Error& err){
                    throw std::runtime_error("roche.stradv: never achieved desired radius");
                }

            },
            "returns position and stream velocity on stream at radius rad

            (x,y,vx,vy) = bspot(q, rad, acc=1.e-7)

            Args:
                q (float): mass ratio = M2/M1
                rad (float): radius to aim for
                acc (float[optional]): computationa accuracy parameter",
            py::arg("q"), py::arg("rad"), py::arg("acc") = 1.e-7
    );    

    m.def("ref_sphere",
            [](double q, double spin, double ffac, int star = 2) {
                double rref, pref;
                Roche::ref_sphere(q, star == 1 ? Roche::STAR::PRIMARY : Roche::STAR::SECONDARY, spin, ffac, &rref, &pref);
                return std::make_tuple(rref, pref);
            },
            "    
            (rref,pref) = ref_sphere(q, spin, ffac, star=2), returns reference radius and potential needed for face.
            q      -- mass ratio = M2/M1
            spin   -- ratio spin/orbital frequency
            ffac   -- linear filling factor of star in question, defined as the radius of the star along the line of
                    centres towards its companion divided by the Roche lobe radius in that direction. For spin = 1
                    the latter is simply the distance to the L1 point, but otherwise you need to use modified L1
                    radii as returned by xl11 or xl12.
            star   -- 1 or 2 for primary or secondary star.",
            py::arg("q"), py::arg("spin"), py::arg("ffac"), py::arg("star") = 2
    );

    m.def("findi",
            [](double q, double deltaphi, double acc = 1.e-4, double di = 1.e-5) {
                // do assertion checks
                if(q <= 0.){
                    throw std::runtime_error("roche.findi: q must be > 0");
                }
                if(deltaphi <= 0. || deltaphi > 0.25){
                    throw std::runtime_error("roche.findi: deltaphi out of range 0 to 0.25");
                }
                if(acc <= 0. || acc > 0.1){
                    throw std::runtime_error("roche.findi: acc <= 0 or acc > 0.1");
                }
                if(di <= 0. || di > 10.){
                    throw std::runtime_error("roche.findi: di <= 0 or di > 10.");
                }

                double ilo = 65., ihi = 90.;
                double phi = deltaphi/2.;
                Subs::Vec3 earth1 = Roche::set_earth(ilo, phi);
                Subs::Vec3 earth2 = Roche::set_earth(ihi, phi);
                Subs::Vec3 r;
                bool elo = Roche::fblink(q, Roche::SECONDARY, 1.0, 1.0, acc, earth1, r);
                bool ehi = Roche::fblink(q, Roche::SECONDARY, 1.0, 1.0, acc, earth2, r);
                double iangle;
                if(elo && ehi){
                    iangle = -2.;
                }else if(!elo && !ehi){
                    iangle = -1.;
                }else{
                    while(ihi - ilo > di){
                        iangle = (ilo+ihi)/2.;
                        if(Roche::fblink(q, Roche::SECONDARY, 1.0, 1.0, acc, Roche::set_earth(iangle, phi), r))
                            ihi = iangle;
                        else
                            ilo = iangle;
                    }
                    iangle = (ilo+ihi)/2.;
                }
                return iangle; 
            },
        "findi(q, deltaphi, acc=1.e-4, di=1.e-5), 
        computes inclination for a given mass ratio and phase width",
        py::arg("q"), py::arg("deltaphi"), py::arg("acc") = 1.e-4, py::arg("di") = 1.e-5
    );

    m.def("findq",
            [](double i, double deltaphi, acc=1.e-4, dq=1.0e-5, qlo=0.001, qhi=2.){
                // do assertion checks
                if(i <= 0. || i > 90.){
                    throw std::runtime_error("roche.findq: i out of range 0 to 90");
                }
                if(deltaphi <= 0. || deltaphi > 0.25){
                    throw std::runtime_error("roche.findq: deltaphi out of range 0 to 0.25");
                }
                if(acc <= 0. || acc > 0.1){
                    throw std::runtime_error("roche.findq: acc <= 0 or acc > 0.1");
                }
                if(dq <= 0. || dq > 0.1){
                    throw std::runtime_error("roche.findq: di <= 0 or di > 0.1");
                }

                double phi = deltaphi/2.;
                Subs::Vec3 r;
                Subs::Vec3 earth = Roche::set_earth(i, phi);
                bool elo = Roche::fblink(qlo, Roche::SECONDARY, 1.0, 1.0, acc, earth, r);
                bool ehi = Roche::fblink(qhi, Roche::SECONDARY, 1.0, 1.0, acc, earth, r);
                double q;
                if(elo && ehi){
                    q = -2.;
                }else if(!elo && !ehi){
                    q = -1.;
                }else{
                    while(qhi - qlo > dq){
                        q = (qlo+qhi)/2.;
                        if(Roche::fblink(q, Roche::SECONDARY, 1.0, 1.0, acc, earth, r))
                            qhi = q;
                        else
                            qlo = q;
                    }
                    q = (qlo+qhi)/2.;
                }
                return q
            },
            "findi(q, deltaphi, acc=1.e-4, di=1.e-5),
            computes mass ratio q for a given phase width and inclination",
            py::arg("i"), py::arg("deltaphi"), py::arg("acc") = 1.e-4, py::arg("dq") = 1.e-5, py::arg("qlo") = 0.005, py::arg("qhi") = 2.
    );

    m.def("findphi",
            [](double q, double i, double delta=1.e-6){
                // do assertion checks
                if(q <= 0.){
                    throw std::runtime_error("roche.findphi: q <= 0");
                }
                if(i <= 0. || i > 90.){
                    throw std::runtime_error("roche.findphi: i out of range 0 to 90");
                }
                if(delta <= 0. || delta > 0.001){
                    throw std::runtime_error("roche.findphi: delta <= 0 or delta > 0.001");
                }
                Subs::Vec3 r(0,0,0);
                double ingress, egress;
                if(!Roche::ingress_egress(q, Roche::SECONDARY, 1.0, 1.0, i, delta, r, ingress, egress)){
                    throw std::runtime_error("roche.findphi: the centre of mass of the white dwarf is not eclipsed");
                }
                return egress - ingress;
            },
            "findphi(q, i, delta=1.e-6), computes deltaphi for a given mass ratio and inclination",
            py::arg("q"), py::arg("i"), py::arg("delta") = 1.e-6
    );

    m.def("fblink",
        [](double q, double iangle, double phi, Subs::Vec3 r, double ffac=1., double acc=1.0e-4, Roche::STAR star=Roche::SECONDARY, int spin = 1){
            // do assertion checks
            if(q <= 0.){
                throw std::runtime_error("roche.fblink: q <= 0");
            }
            if(iangle <= 0. || iangle > 90){
                throw std::runtime_error("roche.fblink: iangle <= 0 or > 90");
            }
            if(ffac < 0. || ffac > 1.0){
                throw std::runtime_error("roche.fblink: ffac < 2");
            }
            if(acc <= 0. || acc > 0.1){
                throw std::runtime_error("roche.fblink: acc <= 0 or acc > 0.1");
            }
            if(star < 1 || star > 2){
                throw std::runtime_error("roche.fblink: star must be either 1 or 2");
            }

            int eclipse;
            if(Roche::fblink(q, star == 1 ? Roche::PRIMARY : Roche::SECONDARY, spin, ffac, acc, Roche::set_earth(iangle, phi), r))
                eclipse = 1;
            else
                eclipse = 0;
            return eclipse;
        }
        "fblink(q, i, phi, r, ffac=1., acc=1.e-4, star=2, spin=1), computes whether a point is eclipsed or not",
        py::arg("q"), py::arg("iangle"), py::arg("phi"), py::arg("r"), py::arg("ffac") = 1., py::arg("acc") = 1.e-4, py::arg("star") = 2, py::arg("spin") = 1
    );

    m.def("ieng",
        [](double q, double iangle, double x, double y, double z=0., double ffac=1., double delta=1.0e-7, double star=2, double spin=1){
            // do assertation checks
            if(q <= 0.){
                throw std::runtime_error("roche.ieng: q <= 0");
            }
            if(iangle <= 0. || iangle > 90){
                throw std::runtime_error("roche.ieng: iangle <= 0 or > 90");
            }
            if(ffac < 0. || ffac > 1.0){
                throw std::runtime_error("roche.ieng: ffac < 2");
            }
            if(delta <= 0.){
                throw std::runtime_error("roche.ieng: delta <= 0");
            }
            if(star < 1 || star > 2){
                throw std::runtime_error("roche.ieng: star must be either 1 or 2");
            }
            Subs::Vec3 r;
            r.set(x, y, z);
            double ingress, egress;
            if(!Roche::ingress_egress(q, star == 1 ? Roche::PRIMARY : Roche::SECONDARY, spin, ffac, iangle, delta, r, ingress, egress)){
                throw std::runtime_error("roche.ieng: point is not eclipsed");
            }
            return std::make_tuple(ingress, egress);
        }
        "ieng(q, iangle, x, y, z=0., ffac=1., delta=1.0e-7, star=2, spin=1), computes ingress and egress phases of a point",
        py::arg("q"), py::arg("iangle"), py::arg("x"), py::arg("y"), py::arg("z") = 0., py::arg("ffac") = 1., py::arg("delta") = 1.0e-7, py::arg("star") = 2, py::arg("spin") = 1
    );

    m.def("lobe1",
        [](double q, int n=200){
            // do assertion checks
            if(q <= 0.){
                throw std::runtime_error("roche.lobe1: q <= 0");
            }
            if(n < 2){
                throw std::runtime_error("roche.lobe1: n < 2");
            }
            double* x[n]
            double* y[n]

            Roche::lobe1(q, x, y, n);
            return std::make_tuple(x, y);
        }
        "lobe1(q, n=200), returns tuple of x, y arrays representing the primary star's Roche lobe",
        py::arg("q"), py::arg("n") = 200
    );

    m.def("lobe2",
        [](double q, int n=200){
            // do assertion checks
            if(q <= 0.){
                throw std::runtime_error("roche.lobe2: q <= 0");
            }
            if(n < 2){
                throw std::runtime_error("roche.lobe2: n < 2");
            }
            double* x[n]
            double* y[n]

            Roche::lobe2(q, x, y, n);
            return std::make_tuple(x, y);
        }
        "lobe2(q, n=200), returns tuple of x, y arrays representing the secondary star's Roche lobe",
        py::arg("q"), py::arg("n") = 200
    );
    
    m.def("rpot",
        [](double q, const Subs::Vec3& r){
            // do assertion checks
            if(q <= 0.){
                throw std::runtime_error("roche.rpot: q <= 0");
            }
            double rp = Roche::rpot(q, r);

            return rp;
        }
        "rpot(q, r), computes Roche potential at a given point",
        py::arg("q"), py::arg("r")
        );
    
    m.def("rpot1",
        [](double q, double spin, const Subs::Vec3& r){
            // do assertion checks
            if(q <= 0.){
                throw std::runtime_error("roche.rpot1: q <= 0");
            }
            double rp = Roche::rpot1(q, spin, r);

            return rp;
        }
        "rpot1(q, spin, r), computes asynchronous Roche potential, star 1",
        py::arg("q"), py::arg("spin"), py::arg("r")
        );
    
    m.def("rpot2",
        [](double q, double spin, const Subs::Vec3& r){
            // do assertion checks
            if(q <= 0.){
                throw std::runtime_error("roche.rpot2: q <= 0");
            }
            double rp = Roche::rpot2(q, spin, r);

            return rp;
        }
        "rpot2(q, spin, r), computes asynchronous Roche potential, star 2",
        py::arg("q"), py::arg("spin"), py::arg("r")
    );

    m.def("drpot",
        [](double q, const Subs::Vec3& r){
            // do assertion checks
            if(q <= 0.){
                throw std::runtime_error("roche.drpot: q <= 0");
            }
            Subs::Vec3 drp = Roche::drpot(q, r);

            return std::make_tuple(drp.x(), drp.y(), drp.z());
        }
        "drpot(q, r), computes derivative of Roche potential at a given point",
        py::arg("q"), py::arg("r")
    );

    m.def("drpot1",
        [](double q, double spin, const Subs::Vec3& r){
            // do assertion checks
            if(q <= 0.){
                throw std::runtime_error("roche.drpot1: q <= 0");
            }
            Subs::Vec3 drp = Roche::drpot1(q, spin, r);

            return std::make_tuple(drp.x(), drp.y(), drp.z());
        }
        "drpot1(q, spin, r), computes derivative of asynchronous Roche potential, star 1",
        py::arg("q"), py::arg("spin"), py::arg("r")
    );

    m.def("drpot2",
        [](double q, double spin, const Subs::Vec3& r){
            // do assertion checks
            if(q <= 0.){
                throw std::runtime_error("roche.drpot2: q <= 0");
            }
            Subs::Vec3 drp = Roche::drpot2(q, spin, r);

            return std::make_tuple(drp.x(), drp.y(), drp.z());
        }
        "drpot2(q, spin, r), computes derivative of asynchronous Roche potential, star 2",
        py::arg("q"), py::arg("spin"), py::arg("r")
    );

    m.def("shadow",
        [](double q, double iangle, double phi, int n=200, double dist=5., double acc=1.e-4){
            // do assertion checks
            if(q <= 0.){
                throw std::runtime_error("roche.shadow: q <= 0");
            }
            if(iangle <= 0. || iangle > 90){
                throw std::runtime_error("roche.shadow: iangle <= 0 or > 90");
            }
            if(n < 2){
                throw std::runtime_error("roche.shadow: n < 2");
            }
            if(dist <= 0.){
                throw std::runtime_error("roche.shadow: dist <= 0");
            }
            if(acc <= 0. || acc > 0.1){
                throw std::runtime_error("roche.shadow: acc <= 0 or acc > 0.1");
            }
            double* x[n];
            double* y[n];
            bool* s[n];
            Roche::roche_shadow(q, iangle, phi, x, y, n, dist, acc);
            return std::make_tuple(x, y, s);
        }
    )


}