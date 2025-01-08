//
// Python/C++ interface file for roche routines.
//

#include <Python.h>
#include "numpy/arrayobject.h"
#include "trm/roche.h"

// Helper functions

static int rconv(PyObject* obj, void* p) {

    // Converter function. Assumes that 'obj' points to a
    // subs.Vec3 from Python and uses it to set the value of
    // a Subs::Vec3 pointed to by p (must be set up by the
    // calling routine). For use in combination with
    // PyArg_ParseTuple and "O&".

    PyObject* px = PyObject_GetAttrString(obj, "x");
    PyObject* py = PyObject_GetAttrString(obj, "y");
    PyObject* pz = PyObject_GetAttrString(obj, "z");

    int status = 0;

    if(px && py && pz){
        // Next functions can convert integers to doubles
        double x = PyFloat_AsDouble(px);
        double y = PyFloat_AsDouble(py);
        double z = PyFloat_AsDouble(pz);
        if(!PyErr_Occurred()){
            ((Subs::Vec3*)p)->set(x,y,z);
            status = 1;
        }
    }

    Py_XDECREF(px);
    Py_XDECREF(py);
    Py_XDECREF(pz);

    return status;
}
    
//----------------------------------------------------------------------------------------
// bspot returns position of bright-spot

static PyObject* 
roche_bspot(PyObject *self, PyObject *args)
{
    double q, rad, acc=1.e-7;
    if(!PyArg_ParseTuple(args, "dd|dd", &q, &rad, &acc))
        return NULL;
    if(q <= 0){
        PyErr_SetString(PyExc_ValueError, "roche.stradv: q <= 0");
        return NULL;
    }
    if(rad <= 0){
        PyErr_SetString(PyExc_ValueError, "roche.stradv: rad <= 0");
        return NULL;
    }
    if(acc <= 0){
        PyErr_SetString(PyExc_ValueError, "roche.stradv: acc <= 0");
        return NULL;
    }

    Subs::Vec3 r, v, rs, vs;
    Roche::strinit(q,r,v);
    try{
        Roche::stradv(q, r, v, rad, acc, 1.e-2);
        return Py_BuildValue("dddd", r.x(), r.y(), v.x(), v.y());
    }
    catch(const Roche::Roche_Error& err){
        PyErr_SetString(PyExc_ValueError, "roche.stradv: never achieved desired radius");
        return NULL;
    }
};

//----------------------------------------------------------------------------------------
// Computes position of a point of given potential in Roche geometry

static PyObject* 
roche_face(PyObject *self, PyObject *args)
{
    
    double q, spin, rref, pref;
    Subs::Vec3 dirn;
    double acc=1.e-9;
    int star = 2;
    if(!PyArg_ParseTuple(args, "ddO&dd|id", &q, &spin, rconv, (void*)&dirn, &rref, &pref, &star, &acc))
        return NULL;
    if(q <= 0.){
        PyErr_SetString(PyExc_ValueError, "roche.face: q <= 0");
        return NULL;
    }
    if(rref <= 0.){
        PyErr_SetString(PyExc_ValueError, "roche.face: rref must be > 0");
        return NULL;
    }
    if(acc <= 0. || acc > 0.1){
        PyErr_SetString(PyExc_ValueError, "roche.face: acc <= 0 or acc > 0.1");
        return NULL;
    }
    if(star < 1 || star > 2){
        PyErr_SetString(PyExc_ValueError, "roche.face: star must be either 1 or 2");
        return NULL;
    }

    // Compute Roche lobe
    Subs::Vec3 pvec, dvec;
    double r, g;
    try{
        Roche::face(q, star == 1 ? Roche::PRIMARY : Roche::SECONDARY, spin, dirn, rref, pref, acc, pvec, dvec, r, g);
        return Py_BuildValue("(ddd)(ddd)dd", pvec.x(), pvec.y(), pvec.z(), dvec.x(), dvec.y(), dvec.z(), r, g);
    }
    catch(const Roche::Roche_Error& err){
        PyErr_SetString(PyExc_ValueError, ("roche.face: " + err).c_str());
        return NULL;
    }
};

//----------------------------------------------------------------------------------------
// Computes the reference potential

static PyObject* 
roche_ref_sphere(PyObject *self, PyObject *args)
{
    
    double q, spin, ffac;
    int star = 2;
    if(!PyArg_ParseTuple(args, "ddd|i", &q, &spin, &ffac, &star))
        return NULL;
    if(q <= 0.){
        PyErr_SetString(PyExc_ValueError, "roche.ref_sphere: q <= 0");
        return NULL;
    }
    if(ffac <= 0. || ffac > 1.){
        PyErr_SetString(PyExc_ValueError, "roche.ref_sphere: require 0 < ffac <= 1");
        return NULL;
    }
    if(star < 1 || star > 2){
        PyErr_SetString(PyExc_ValueError, "roche.ref_sphere: star must be either 1 or 2");
        return NULL;
    }

    // Compute Roche lobe
    double rref, pref;
    Roche::ref_sphere(q, star == 1 ? Roche::PRIMARY : Roche::SECONDARY, spin, ffac, rref, pref);
    return Py_BuildValue("dd", rref, pref);
};


//----------------------------------------------------------------------------------------
// Finds inclination corresponding to a given mass ratio and white dwarf eclipse phase width
// returns -1 if there is no valid value

static PyObject* 
roche_findi(PyObject *self, PyObject *args)
{
    
    double q, pwidth, acc=1.e-4, di=1.e-5;
    if(!PyArg_ParseTuple(args, "dd|dddd", &q, &pwidth, &acc, &di))
        return NULL;
    if(q <= 0.){
        PyErr_SetString(PyExc_ValueError, "roche.findi: q must be > 0");
        return NULL;
    }
    if(pwidth <= 0. || pwidth > 0.25){
        PyErr_SetString(PyExc_ValueError, "roche.findi: pwidth out of range 0 to 0.25");
        return NULL;
    }
    if(acc <= 0. || acc > 0.1){
        PyErr_SetString(PyExc_ValueError, "roche.findi: acc <= 0 or acc > 0.1");
        return NULL;
    }
    if(di <= 0. || di > 10.){
        PyErr_SetString(PyExc_ValueError, "roche.findi: di <= 0 or di > 10.");
        return NULL;
    }

    double ilo = 65., ihi = 90.;
    double phi = pwidth/2.;
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
    return Py_BuildValue("d", iangle);
};

//----------------------------------------------------------------------------------------
// Finds mass ratio corresponding to a given angle and white dwarf eclipse phase width

static PyObject* 
roche_findq(PyObject *self, PyObject *args)
{
    
    double iangle, pwidth, acc=1.e-4, dq=1.e-5, qlo=0.001, qhi=2.0;
    if(!PyArg_ParseTuple(args, "dd|dddd", &iangle, &pwidth, &acc, &dq, &qlo, &qhi))
        return NULL;
    if(iangle <= 0. || iangle > 90.){
        PyErr_SetString(PyExc_ValueError, "roche.findq: iangle out of range 0 to 90");
        return NULL;
    }
    if(pwidth <= 0. || pwidth > 0.25){
        PyErr_SetString(PyExc_ValueError, "roche.findq: pwidth out of range 0 to 0.25");
        return NULL;
    }
    if(acc <= 0. || acc > 0.1){
        PyErr_SetString(PyExc_ValueError, "roche.findq: acc <= 0 or acc > 0.1");
        return NULL;
    }
    if(dq <= 0. || dq > 0.1){
        PyErr_SetString(PyExc_ValueError, "roche.findq: dq <= 0 or dq > 0.1");
        return NULL;
    }

    Subs::Vec3 r;
    double phi = pwidth/2.;
    Subs::Vec3 earth = Roche::set_earth(iangle, phi);
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
    return Py_BuildValue("d", q);
};

//----------------------------------------------------------------------------------------
// Finds mass ratio corresponding to a given angle and white dwarf eclipse phase width

static PyObject* 
roche_findphi(PyObject *self, PyObject *args)
{
    
    double q, iangle, delta=1.e-6;
    if(!PyArg_ParseTuple(args, "dd|dddd", &q, &iangle, &delta))
        return NULL;
    if(q <= 0.){
        PyErr_SetString(PyExc_ValueError, "roche.findphi: q <= 0");
        return NULL;
    }
    if(iangle <= 0. || iangle > 90.){
        PyErr_SetString(PyExc_ValueError, "roche.findphi: iangle out of range 0 to 90");
        return NULL;
    }
    if(delta <= 0. || delta > 0.001){
        PyErr_SetString(PyExc_ValueError, "roche.findphi: delta <= 0 or delta > 0.001");
        return NULL;
    }

    Subs::Vec3 r(0,0,0);
    double ingress, egress;
    if(!Roche::ingress_egress(q, Roche::SECONDARY, 1.0, 1.0, iangle, delta, r, ingress, egress)){
        PyErr_SetString(PyExc_ValueError, "roche.findphi: the centre of mass of the white dwarf is not eclipsed");
        return NULL;
    }
    return Py_BuildValue("d", egress-ingress);
};

//----------------------------------------------------------------------------------------
// Computes whether a given point is eclipsed

static PyObject* 
roche_fblink(PyObject *self, PyObject *args)
{
    
    double q, iangle, phi;
    Subs::Vec3 r;
    double ffac=1., acc=1.e-4, spin = 1.0;
    int star = 2;
    if(!PyArg_ParseTuple(args, "dddO&|ddid", &q, &iangle, &phi,
                         rconv, (void*)&r, &ffac, &acc, &star, &spin))
        return NULL;
    if(q <= 0.){
        PyErr_SetString(PyExc_ValueError, "roche.fblink: q <= 0");
        return NULL;
    }
    if(iangle <= 0. || iangle > 90.){
        PyErr_SetString(PyExc_ValueError, "roche.fblink: iangle out of range 0 to 90");
        return NULL;
    }
    if(ffac <= 0. || ffac > 1.){
        PyErr_SetString(PyExc_ValueError, "roche.fblink: ffac out of range 0 to 1");
        return NULL;
    }
    if(acc <= 0. || acc > 0.1){
        PyErr_SetString(PyExc_ValueError, "roche.fblink: acc <= 0 or acc > 0.1");
        return NULL;
    }
    if(star < 1 || star > 2){
        PyErr_SetString(PyExc_ValueError, "roche.fblink: star must be either 1 or 2");
        return NULL;
    }

    // Compute Roche lobe
    int eclipse;
    if(Roche::fblink(q, star == 1 ? Roche::PRIMARY : Roche::SECONDARY, spin, ffac, acc, Roche::set_earth(iangle, phi), r))
        eclipse = 1;
    else
        eclipse = 0;
    return Py_BuildValue("i", eclipse);
};

//----------------------------------------------------------------------------------------
// Computes inegress and egress phases of a point

static PyObject* 
roche_ineg(PyObject *self, PyObject *args)
{
    
    double q, iangle, x, y, z=0.;
    double ffac=1., delta=1.e-7, spin = 1.0;
    int star = 2;
    if(!PyArg_ParseTuple(args, "dddd|dddid", &q, &iangle, &x, &y, &z, &ffac, &delta, &star, &spin))
        return NULL;
    if(q <= 0.){
        PyErr_SetString(PyExc_ValueError, "roche.ineg: q <= 0");
        return NULL;
    }
    if(iangle <= 0. || iangle > 90.){
        PyErr_SetString(PyExc_ValueError, "roche.ineg: iangle out of range 0 to 90");
        return NULL;
    }
    if(ffac <= 0. || ffac > 1.){
        PyErr_SetString(PyExc_ValueError, "roche.ineg: ffac out of range 0 to 1");
        return NULL;
    }
    if(delta <= 0.){
        PyErr_SetString(PyExc_ValueError, "roche.ineg: delta <= 0");
        return NULL;
    }
    if(star < 1 || star > 2){
        PyErr_SetString(PyExc_ValueError, "roche.ineg: star must be either 1 or 2");
        return NULL;
    }

    // Compute Roche lobe
    Subs::Vec3 r(x,y,z);
    double ingress, egress;
    if(!Roche::ingress_egress(q, star == 1 ? Roche::PRIMARY : Roche::SECONDARY, spin, ffac, iangle, delta, r, ingress, egress)){
        PyErr_SetString(PyExc_ValueError, "roche.ineg: point is not eclipsed");
        return NULL;
    }
    return Py_BuildValue("dd", ingress, egress);
};

//----------------------------------------------------------------------------------------
// Returns tuple of x, y arrays representing the secondary star's Roche lobe

static PyObject* 
roche_lobe1(PyObject *self, PyObject *args)
{

    double q;
    int n = 200;
    if(!PyArg_ParseTuple(args, "d|i:roche.lobe1", &q, &n))
        return NULL;
    if(q <= 0.){
        PyErr_SetString(PyExc_ValueError, "roche.lobe1: q <= 0");
        return NULL;
    }
    if(n < 2){
        PyErr_SetString(PyExc_ValueError, "roche.lobe1: n < 2");
        return NULL;
    }

    // Create output array containing both x and y arrays.
    npy_intp dim[1] = {n};

    PyArrayObject* xo = NULL;
    xo = (PyArrayObject*) PyArray_SimpleNew(1, dim, PyArray_FLOAT);
    if(xo == NULL) return NULL;

    PyArrayObject* yo = NULL;
    yo = (PyArrayObject*) PyArray_SimpleNew(1, dim, PyArray_FLOAT);
    if(yo == NULL){
        Py_DECREF(xo);
        return NULL;
    }

    float* x = (float*)xo->data;
    float* y = (float*)yo->data;

    // Compute Roche lobe
    Roche::lobe1(q, x, y, n);

    return Py_BuildValue("NN", xo, yo);

};

//----------------------------------------------------------------------------------------
// Returns tuple of x, y arrays representing the secondary star's Roche lobe

static PyObject* 
roche_lobe2(PyObject *self, PyObject *args)
{

    double q;
    int n = 200;
    if(!PyArg_ParseTuple(args, "d|i:roche.lobe2", &q, &n))
        return NULL;
    if(q <= 0.){
        PyErr_SetString(PyExc_ValueError, "roche.lobe2: q <= 0");
        return NULL;
    }
    if(n < 2){
        PyErr_SetString(PyExc_ValueError, "roche.lobe2: n < 2 in");
        return NULL;
    }

    // Create output array containing both x and y arrays.
    npy_intp dim[1] = {n};

    PyArrayObject* xo = NULL;
    xo = (PyArrayObject*) PyArray_SimpleNew(1, dim, PyArray_FLOAT);
    if(xo == NULL) return NULL;

    PyArrayObject* yo = NULL;
    yo = (PyArrayObject*) PyArray_SimpleNew(1, dim, PyArray_FLOAT);
    if(yo == NULL){
        Py_DECREF(xo);
        return NULL;
    }

    float* x = (float*)xo->data;
    float* y = (float*)yo->data;

    // Compute Roche lobe
    Roche::lobe2(q, x, y, n);

    return Py_BuildValue("NN", xo, yo);

};

//----------------------------------------------------------------------------------------
// Computes Roche potential at a given point

static PyObject*
roche_rpot(PyObject *self, PyObject *args)
{

    double q;
    Subs::Vec3 p;
    if(!PyArg_ParseTuple(args, "dO&", &q, rconv, (void*)&p))
        return NULL;

    if(q <= 0.){
        PyErr_SetString(PyExc_ValueError, "roche.rpot: q <= 0");
        return NULL;
    }

    // Compute Roche lobe
    double rp = Roche::rpot(q, p);
    return Py_BuildValue("d", rp);
};

//----------------------------------------------------------------------------------------
// Computes Roche potential for star 1 asynchronous

static PyObject*
roche_rpot1(PyObject *self, PyObject *args)
{

    double q, spin;
    Subs::Vec3 p;
    if(!PyArg_ParseTuple(args, "ddO&", &q, &spin, rconv, (void*)&p))
        return NULL;
    if(q <= 0.){
        PyErr_SetString(PyExc_ValueError, "roche.rpot1: q <= 0");
        return NULL;
    }

    double rp = Roche::rpot1(q, spin, p);
    return Py_BuildValue("d", rp);
};

//----------------------------------------------------------------------------------------
// Computes Roche potential for star 2 asynchronous

static PyObject*
roche_rpot2(PyObject *self, PyObject *args)
{

    double q, spin;
    Subs::Vec3 p;
    if(!PyArg_ParseTuple(args, "ddO&", &q, &spin, rconv, (void*)&p))
        return NULL;
    if(q <= 0.){
        PyErr_SetString(PyExc_ValueError, "roche.rpot2: q <= 0");
        return NULL;
    }

    double rp = Roche::rpot2(q, spin, p);
    return Py_BuildValue("d", rp);
};

//----------------------------------------------------------------------------------------
// Computes derivative of Roche potential at a given point

static PyObject*
roche_drpot(PyObject *self, PyObject *args)
{

    double q;
    Subs::Vec3 p;
    if(!PyArg_ParseTuple(args, "dO&", &q, rconv, (void*)&p))
        return NULL;

    if(q <= 0.){
        PyErr_SetString(PyExc_ValueError, "roche.drpot: q <= 0");
        return NULL;
    }

    // Compute Roche lobe
    Subs::Vec3 drp = Roche::drpot(q, p);
    return Py_BuildValue("ddd", drp.x(), drp.y(), drp.z());
};

//----------------------------------------------------------------------------------------
// Computes derivatives of Roche potential for star 1 asynchronous

static PyObject*
roche_drpot1(PyObject *self, PyObject *args)
{

    double q, spin;
    Subs::Vec3 p;
    if(!PyArg_ParseTuple(args, "ddO&", &q, &spin, rconv, (void*)&p))
        return NULL;
    if(q <= 0.){
        PyErr_SetString(PyExc_ValueError, "roche.drpot1: q <= 0");
        return NULL;
    }

    Subs::Vec3 drp = Roche::drpot1(q, spin, p);
    return Py_BuildValue("ddd", drp.x(), drp.y(), drp.z());
};

//----------------------------------------------------------------------------------------
// Returns tuple of x, y arrays representing the "shadow" of the secondary star's Roche lobe

static PyObject*
roche_shadow(PyObject *self, PyObject *args)
{

    double q, iangle, phi, dist=5., acc=1.e-4;
    int n = 200;
    if(!PyArg_ParseTuple(args, "ddd|idd:roche.shadow", &q, &iangle, &phi, &n, &dist, &acc))
        return NULL;
    if(q <= 0.){
        PyErr_SetString(PyExc_ValueError, "roche.shadow: q <= 0");
        return NULL;
    }
    if(iangle <= 0. || iangle > 90){
        PyErr_SetString(PyExc_ValueError, "roche.shadow: iangle <= 0 or > 90");
        return NULL;
    }
    if(n < 2){
        PyErr_SetString(PyExc_ValueError, "roche.shadow: n < 2");
        return NULL;
    }
    if(dist <= 0.){
        PyErr_SetString(PyExc_ValueError, "roche.shadow: dist <= 0");
        return NULL;
    }
    if(acc <= 0. || acc > 0.1){
        PyErr_SetString(PyExc_ValueError, "roche.shadow: acc <= 0 or > 0.1");
        return NULL;
    }

    // Create output array containing both x and y arrays.
    npy_intp dim[1] = {n};

    PyArrayObject* xo = NULL;
    xo = (PyArrayObject*) PyArray_SimpleNew(1, dim, PyArray_FLOAT);
    if(xo == NULL) return NULL;

    PyArrayObject* yo = NULL;
    yo = (PyArrayObject*) PyArray_SimpleNew(1, dim, PyArray_FLOAT);
    if(yo == NULL){
        Py_DECREF(xo);
        return NULL;
    }

    PyArrayObject* so = NULL;
    so = (PyArrayObject*) PyArray_SimpleNew(1, dim, PyArray_BOOL);
    if(so == NULL){
        Py_DECREF(xo);
        Py_DECREF(yo);
        return NULL;
    }

    float* x = (float*)xo->data;
    float* y = (float*)yo->data;
    bool*  s = (bool*)so->data; 

    // Compute Roche lobe
    Roche::roche_shadow(q, iangle, phi, dist, acc, x, y, s, n);

    return Py_BuildValue("NNN",xo,yo,so);

};


//----------------------------------------------------------------------------------------
// Returns tuple of x, y arrays representing the gas stream

static PyObject* 
roche_streamr(PyObject *self, PyObject *args)
{

    double q, rad;
    int n = 200;
    if(!PyArg_ParseTuple(args, "dd|i:roche.streamr", &q, &rad, &n))
        return NULL;
    if(q <= 0.){
        PyErr_SetString(PyExc_ValueError, "roche.streamr: q <= 0");
        return NULL;
    }
    if(rad < 0. || rad > 1.){
        PyErr_SetString(PyExc_ValueError, "roche.streamr: rad < 0 or > 1.");
        return NULL;
    }
    if(n < 2){
        PyErr_SetString(PyExc_ValueError, "roche.streamr: n < 2");
        return NULL;
    }

    // Create output array containing both x and y arrays.
    npy_intp dim[1] = {n};

    PyArrayObject* xo = NULL;
    xo = (PyArrayObject*) PyArray_SimpleNew(1, dim, PyArray_FLOAT);
    if(xo == NULL) return NULL;

    PyArrayObject* yo = NULL;
    yo = (PyArrayObject*) PyArray_SimpleNew(1, dim, PyArray_FLOAT);
    if(yo == NULL){
        Py_DECREF(xo);
        return NULL;
    }

    float* x = (float*)xo->data;
    float* y = (float*)yo->data;

    // Carry out the stream integration
    try{
        Roche::streamr(q, rad, x, y, n);
    }
    catch(const Roche::Roche_Error &err){
        PyErr_SetString(PyExc_ValueError, ("roche.streamr: " + err).c_str());
        return NULL;
    }

    return Py_BuildValue("NN", xo, yo);

};

static PyObject* 
roche_stream(PyObject *self, PyObject *args)
{

    double q, step;
    int n = 200;
    if(!PyArg_ParseTuple(args, "dd|i:roche.stream", &q, &step, &n))
        return NULL;
    if(q <= 0.){
        PyErr_SetString(PyExc_ValueError, "roche.stream: q <= 0");
        return NULL;
    }
    if(step <= 0. || step > 1.){
        PyErr_SetString(PyExc_ValueError, "roche.stream: step <= 0 or > 1.");
        return NULL;
    }
    if(n < 2){
        PyErr_SetString(PyExc_ValueError, "roche.stream: n < 2");
        return NULL;
    }

    // Create output array containing both x and y arrays.
    npy_intp dim[1] = {n};

    PyArrayObject* xo = NULL;
    xo = (PyArrayObject*) PyArray_SimpleNew(1, dim, PyArray_FLOAT);
    if(xo == NULL) return NULL;

    PyArrayObject* yo = NULL;
    yo = (PyArrayObject*) PyArray_SimpleNew(1, dim, PyArray_FLOAT);
    if(yo == NULL){
        Py_DECREF(xo);
        return NULL;
    }

    float* x = (float*)xo->data;
    float* y = (float*)yo->data;

    // Carry out the stream integration
    try{
        Roche::stream(q, step, x, y, n);
    }
    catch(const Roche::Roche_Error &err){
        PyErr_SetString(PyExc_ValueError, ("roche.stream: " + err).c_str());
        return NULL;
    }

    return Py_BuildValue("NN", xo, yo);

};

/* arbitrary stream integration from user-defined starting point */

static PyObject*
roche_astream(PyObject *self, PyObject *args)
{
    double q, step;
    double acc=1.e-9;
    int n = 200, type;
    Subs::Vec3 r, v;
    if(!PyArg_ParseTuple(args, "diO&O&d|id", &q, &type, rconv, (void*)&r,
                         rconv, (void*)&v, &step, &n, &acc))
        return NULL;
    if(q <= 0.){
        PyErr_SetString(PyExc_ValueError, "roche.astream: q <= 0");
        return NULL;
    }
    if(type < 0 || type > 3){
        PyErr_SetString(PyExc_ValueError, "roche.astream: type out of range 0 to 3");
        return NULL;
    }
    if(step <= 0. || step > 1.){
        PyErr_SetString(PyExc_ValueError, "roche.astream: step <= 0 or > 1.");
        return NULL;
    }
    if(n < 2){
        PyErr_SetString(PyExc_ValueError, "roche.astream: n < 2");
        return NULL;
    }
    if(acc <= 0. || acc > 0.1){
        PyErr_SetString(PyExc_ValueError, 
                        "roche.astream: acc <= 0 or acc > 0.1");
        return NULL;
    }

    // Create output array containing both x and y arrays.
    npy_intp dim[1] = {n};

    PyArrayObject* xo = NULL;
    xo = (PyArrayObject*) PyArray_SimpleNew(1, dim, PyArray_FLOAT);
    if(xo == NULL) return NULL;

    PyArrayObject* yo = NULL;
    yo = (PyArrayObject*) PyArray_SimpleNew(1, dim, PyArray_FLOAT);
    if(yo == NULL){
        Py_DECREF(xo);
        return NULL;
    }

    float* ax = (float*)xo->data;
    float* ay = (float*)yo->data;

    // Carry out the stream integration
    try{
        double xold = r.x(), yold = r.y(), apx, apy;
        if(type == 0){
            ax[0] = xold;
            ay[0] = yold;
        }else{
            Roche::vtrans(q, type, r.x(), r.y(), v.x(), v.y(), apx, apy);
            ax[0] = apx;
            ay[0] = apy;
        }
        int lp = 0;
	
        double time, dist, tdid, tnext, frac, ttry;
        double vel, smax = std::min(1.e-3,step/2.);

        vel  = sqrt(Subs::sqr(v.x()) + Subs::sqr(v.y()));
        ttry = smax/std::max(1.e20, vel);
        while(lp < n-1){
            Roche::gsint(q, r, v, ttry, tdid, tnext, time, acc);
            dist = sqrt(Subs::sqr(r.x()-xold)+Subs::sqr(r.y()-yold));
            if(dist > step){
                frac = step / dist;
                if(type == 0){
                    ax[lp+1] = ax[lp] + (r.x()-ax[lp])*frac;
                    ay[lp+1] = ay[lp] + (r.y()-ay[lp])*frac;
                }else{
                    Roche::vtrans(q, type, r.x(), r.y(), v.x(), v.y(), apx, apy);
                    ax[lp+1] = ax[lp] + (apx-ax[lp])*frac;
                    ay[lp+1] = ay[lp] + (apy-ay[lp])*frac;
                }
                xold = r.x();
                yold = r.y();
                lp++;		
            }
            vel  = sqrt(Subs::sqr(v.x()) + Subs::sqr(v.y()));
            ttry = std::min(smax/vel, tnext);
        }
    }
    catch(const Roche::Roche_Error &err){
        PyErr_SetString(PyExc_ValueError, ("roche.astream: " + err).c_str());
        return NULL;
    }

    return Py_BuildValue("NN", xo, yo);

};

//----------------------------------------------------------------------------------------
// strmnx

static PyObject* 
roche_strmnx(PyObject *self, PyObject *args)
{
    double q, acc=1.e-7;
    int n = 1;
    if(!PyArg_ParseTuple(args, "d|id", &q, &n, &acc))
        return NULL;
    if(q < 0){
        PyErr_SetString(PyExc_ValueError, "roche.strmnx: q <= 0");
        return NULL;
    }
    if(n < 1){
        PyErr_SetString(PyExc_ValueError, "roche.strmnx: n < 1");
        return NULL;
    }
    if(acc <= 0.){
        PyErr_SetString(PyExc_ValueError, "roche.strmnx: acc <= 0");
        return NULL;
    }
    Subs::Vec3 r, v;
    Roche::strinit(q, r, v);
    for(int i=0; i<n; i++)
        Roche::strmnx(q, r, v, acc);
    double tvx1, tvy1, tvx2, tvy2;
    Roche::vtrans(q, 1, r.x(), r.y(), v.x(), v.y(), tvx1, tvy1);
    Roche::vtrans(q, 2, r.x(), r.y(), v.x(), v.y(), tvx2, tvy2);
    return Py_BuildValue("dddddd", r.x(), r.y(), tvx1, tvy1, tvx2, tvy2);
};

//----------------------------------------------------------------------------------------
// Returns tuple of x, y arrays representing the primary star's Roche lobe in velocity space

static PyObject* 
roche_vlobe1(PyObject *self, PyObject *args)
{

    double q;
    int n = 200;
    if(!PyArg_ParseTuple(args, "d|i:roche.vlobe1", &q, &n))
        return NULL;
    if(q <= 0.){
        PyErr_SetString(PyExc_ValueError, "roche.vlobe1: q <= 0");
        return NULL;
    }
    if(n < 2){
        PyErr_SetString(PyExc_ValueError, "roche.vlobe1: n < 2");
        return NULL;
    }

    // Create output array containing both x and y arrays.
    npy_intp dim[1] = {n};

    PyArrayObject* vxo = NULL;
    vxo = (PyArrayObject*) PyArray_SimpleNew(1, dim, PyArray_FLOAT);
    if(vxo == NULL) return NULL;

    PyArrayObject* vyo = NULL;
    vyo = (PyArrayObject*) PyArray_SimpleNew(1, dim, PyArray_FLOAT);
    if(vyo == NULL){
        Py_DECREF(vxo);
        return NULL;
    }

    float* vx = (float*)vxo->data;
    float* vy = (float*)vyo->data;

    // Compute Roche lobe
    Roche::vlobe1(q, vx, vy, n);

    return Py_BuildValue("NN", vxo, vyo);

};


//----------------------------------------------------------------------------------------
// Returns tuple of vx, vy arrays representing the secondary star's Roche lobe in velocity space

static PyObject* 
roche_vlobe2(PyObject *self, PyObject *args)
{

    double q;
    int n = 200;
    if(!PyArg_ParseTuple(args, "d|i:roche.vlobe2", &q, &n))
        return NULL;
    if(q <= 0.){
        PyErr_SetString(PyExc_ValueError, "roche.vlobe2: q <= 0");
        return NULL;
    }
    if(n < 2){
        PyErr_SetString(PyExc_ValueError, "roche.vlobe2: n < 2");
        return NULL;
    }

    // Create output array containing both x and y arrays.
    npy_intp dim[1] = {n};

    PyArrayObject* vxo = NULL;
    vxo = (PyArrayObject*) PyArray_SimpleNew(1, dim, PyArray_FLOAT);
    if(vxo == NULL) return NULL;

    PyArrayObject* vyo = NULL;
    vyo = (PyArrayObject*) PyArray_SimpleNew(1, dim, PyArray_FLOAT);
    if(vyo == NULL){
        Py_DECREF(vxo);
        return NULL;
    }

    float* vx = (float*)vxo->data;
    float* vy = (float*)vyo->data;

    // Compute Roche lobe
    Roche::vlobe2(q, vx, vy, n);

    return Py_BuildValue("NN", vxo, vyo);

};

//----------------------------------------------------------------------------------------
// Returns tuple of vx, vy arrays representing the gas stream

static PyObject* 
roche_vstream(PyObject *self, PyObject *args)
{

    double q, step = 0.01;
    int stype=1, n = 60;
    if(!PyArg_ParseTuple(args, "d|dii:roche.vstream", &q, &step, &stype, &n))
        return NULL;
    if(q <= 0.){
        PyErr_SetString(PyExc_ValueError, "roche.vstream: q <= 0");
        return NULL;
    }
    if(step <= 0. || step >= 1.){
        PyErr_SetString(PyExc_ValueError, "roche.vstream: step <= 0 or >= 1.");
        return NULL;
    }
    if(stype < 1 && stype > 3){
        PyErr_SetString(PyExc_ValueError, "roche.vstream: stype must be 1, 2 or 3");
        return NULL;
    }
    if(n < 2){
        PyErr_SetString(PyExc_ValueError, "roche.vstream: n < 2");
        return NULL;
    }

    // Create output array containing both x and y arrays.
    npy_intp dim[1] = {n};

    PyArrayObject* vxo = NULL;
    vxo = (PyArrayObject*) PyArray_SimpleNew(1, dim, PyArray_FLOAT);
    if(vxo == NULL) return NULL;

    PyArrayObject* vyo = NULL;
    vyo = (PyArrayObject*) PyArray_SimpleNew(1, dim, PyArray_FLOAT);
    if(vyo == NULL){
        Py_DECREF(vxo);
        return NULL;
    }

    float* vx = (float*)vxo->data;
    float* vy = (float*)vyo->data;


    // Carry out the stream integration
    try{
        Roche::vstrreg(q, step, vx, vy, n, stype);
    }
    catch(const Roche::Roche_Error &err){
        PyErr_SetString(PyExc_ValueError, ("roche.vstream: " + err).c_str());
        return NULL;
    }

    return Py_BuildValue("NN", vxo, vyo);

};

//------------------------------------------------------------------------
// Returns tuple of x, y, vx, vy, time and jacobi as arrays
// along the gas stream

static PyObject*
roche_pvstream(PyObject *self, PyObject *args)
{

    double q, step = 0.01;
    int stype=1, n = 60;
    if(!PyArg_ParseTuple(args, "d|dii:roche.pvstream",
                         &q, &step, &stype, &n))
        return NULL;
    if(q <= 0.){
        PyErr_SetString(PyExc_ValueError, "roche.pvstream: q <= 0");
        return NULL;
    }
    if(step <= 0. || step >= 1.){
        PyErr_SetString(PyExc_ValueError,
                        "roche.pvstream: step <= 0 or >= 1.");
        return NULL;
    }
    if(stype < 1 && stype > 3){
        PyErr_SetString(PyExc_ValueError,
                        "roche.pvstream: stype must be 1, 2 or 3");
        return NULL;
    }
    if(n < 2){
        PyErr_SetString(PyExc_ValueError, "roche.pvstream: n < 2");
        return NULL;
    }

    // Create output arrays
    npy_intp dim[1] = {n};

    PyArrayObject* xo = NULL;
    xo = (PyArrayObject*) PyArray_SimpleNew(1, dim, PyArray_FLOAT);
    if(xo == NULL) return NULL;

    PyArrayObject* yo = NULL;
    yo = (PyArrayObject*) PyArray_SimpleNew(1, dim, PyArray_FLOAT);
    if(yo == NULL){
        Py_DECREF(xo);
        return NULL;
    }

    PyArrayObject* vxo = NULL;
    vxo = (PyArrayObject*) PyArray_SimpleNew(1, dim, PyArray_FLOAT);
    if(vxo == NULL){
        Py_DECREF(xo);
        Py_DECREF(yo);
        return NULL;
    }

    PyArrayObject* vyo = NULL;
    vyo = (PyArrayObject*) PyArray_SimpleNew(1, dim, PyArray_FLOAT);
    if(vyo == NULL){
        Py_DECREF(xo);
        Py_DECREF(yo);
        Py_DECREF(vxo);
        return NULL;
    }

    PyArrayObject* to = NULL;
    to = (PyArrayObject*) PyArray_SimpleNew(1, dim, PyArray_FLOAT);
    if(to == NULL){
        Py_DECREF(xo);
        Py_DECREF(yo);
        Py_DECREF(vxo);
        Py_DECREF(vyo);
        return NULL;
    }

    PyArrayObject* jco = NULL;
    jco = (PyArrayObject*) PyArray_SimpleNew(1, dim, PyArray_FLOAT);
    if(jco == NULL){
        Py_DECREF(xo);
        Py_DECREF(yo);
        Py_DECREF(vxo);
        Py_DECREF(vyo);
        Py_DECREF(to);
        return NULL;
    }

    float* x  = (float*)xo->data;
    float* y  = (float*)yo->data;
    float* vx = (float*)vxo->data;
    float* vy = (float*)vyo->data;
    float* t  = (float*)to->data;
    float* jc = (float*)jco->data;

    // Carry out the stream integration
    try{
        const double TLOC =  1.e-8;
        const double RLOC =  1.e-8;
        int np;
        double dt, rl1, tvx, tvy, rend, rnext;
        Subs::Vec3 r, v, rm, vm;
        int decr;

        rl1 = Roche::xl1(q);

        /* Store L1 as first point */
        Roche::vtrans(q, stype, rl1, 0., 0., 0., tvx, tvy);
        x[0]  = rl1;
        y[0]  = 0.;
        vx[0] = tvx;
        vy[0] = tvy;
        t[0]  = 0.;
        np    = 1;
        rnext = rl1*(1.-step);
        decr  = 1;

        /* Initialise stream */
        Roche::strinit(q, r, v);
        jc[0] = Roche::jacobi(q, r, v);

        while(np < n){

            /* Advance one step */
            dt = Roche::stradv(q, r, v, rnext, RLOC, 1.e-3);
            Roche::vtrans(q, stype, r.x(), r.y(), v.x(), v.y(), tvx, tvy);
            x[np]  = r.x();
            y[np]  = r.y();
            vx[np] = tvx;
            vy[np] = tvy;
            t[np]  = t[np-1] + dt;
            jc[np] = Roche::jacobi(q, r, v);
            np++;
            rnext = decr ? rnext - rl1*step : rnext + rl1*step;

            /* Locate and store next turning point */
            rm = r;
            vm = v;
            Roche::strmnx(q, rm, vm, TLOC);
            rend = rm.length();

            /* Loop over all radii wanted before next turning point */
            while(np < n && ((decr && rnext > rend) ||
                             (!decr && rnext < rend))){
                dt = Roche::stradv(q, r, v, rnext, RLOC, 1.e-3);
                Roche::vtrans(q, stype, r.x(), r.y(), v.x(), v.y(), tvx, tvy);
                x[np]  = r.x();
                y[np]  = r.y();
                vx[np] = tvx;
                vy[np] = tvy;
                t[np]  = t[np-1] + dt;
                jc[np] = Roche::jacobi(q, r, v);
                np++;
                rnext = decr ? rnext - rl1*step : rnext + rl1*step;
            }

            /* Change direction of search, and move it to start at turning
               point */
            rnext = decr ? rnext + rl1*step : rnext - rl1*step;
            r     = rm;
            v     = vm;
            decr  = !decr;
        }
    }
    catch(const Roche::Roche_Error &err){
        PyErr_SetString(PyExc_ValueError, ("roche.pvstream: " + err).c_str());
        return NULL;
    }

    return Py_BuildValue("NNNNNN", xo, yo, vxo, vyo, to, jco);

};

//----------------------------------------------------------------------------------------
// xl1

static PyObject* 
roche_xl1(PyObject *self, PyObject *args)
{
    double q;
    if(!PyArg_ParseTuple(args, "d:roche.xl1", &q))
        return NULL;
    if(q < 0){
        PyErr_SetString(PyExc_ValueError, "roche.xl1: q <= 0");
        return NULL;
    }
    double x = Roche::xl1(q);
    return Py_BuildValue("f", x);
};

//----------------------------------------------------------------------------------------
// xl2

static PyObject* 
roche_xl2(PyObject *self, PyObject *args)
{
    double q;
    if(!PyArg_ParseTuple(args, "d:roche.xl2", &q))
        return NULL;
    if(q < 0){
        PyErr_SetString(PyExc_ValueError, "roche.xl2: q <= 0");
        return NULL;
    }
    double x = Roche::xl2(q);
    return Py_BuildValue("f", x);
};

//----------------------------------------------------------------------------------------
// xl3

static PyObject* 
roche_xl3(PyObject *self, PyObject *args)
{
    double q;
    if(!PyArg_ParseTuple(args, "d:roche.xl3", &q))
        return NULL;
    if(q < 0){
        PyErr_SetString(PyExc_ValueError, "roche.cl3: q <= 0");
        return NULL;
    }
    double x = Roche::xl3(q);
    return Py_BuildValue("f", x);
};

//----------------------------------------------------------------------------------------
// xl11

static PyObject* 
roche_xl11(PyObject *self, PyObject *args)
{
    double q, spin;
    if(!PyArg_ParseTuple(args, "dd:roche.xl11", &q, &spin))
        return NULL;
    if(q < 0){
        PyErr_SetString(PyExc_ValueError, "roche.xl11: q <= 0");
        return NULL;
    }
    double x = Roche::xl11(q,spin);
    return Py_BuildValue("f", x);
};

//----------------------------------------------------------------------------------------
// xl12

static PyObject* 
roche_xl12(PyObject *self, PyObject *args)
{
    double q, spin;
    if(!PyArg_ParseTuple(args, "dd:roche.xl12", &q, &spin))
        return NULL;
    if(q < 0){
        PyErr_SetString(PyExc_ValueError, "roche.xl12: q <= 0");
        return NULL;
    }
    double x = Roche::xl12(q,spin);
    return Py_BuildValue("f", x);
};

//----------------------------------------------------------------------------------------
// The methods

static PyMethodDef RocheMethods[] = {

    {"astream", roche_astream, METH_VARARGS, 
     "(x,y) = astream(q, type, r, v, step, n=200, acc=1.e-9), returns\n"
     "arrays of the gas stream given arbitrary initial conditions, r, v\n"
     "q = M2/M1, step=distance between adjacent points, n= number of points, acc=\n."
     "accuracy of calculations. type =0, 1, 2 or 3. 0 gives x,y positions, 1,2,3 give\n"
     "different velocities -- see vstream for more detail.\n."
    },

    {"bspot", roche_bspot, METH_VARARGS, 
     "(x,y,vx,vy) = bspot(q, rad, acc=1.e-7), returns position and stream velocity on stream at radius rad\n\n"
     " q    -- mass ratio = M2/M1\n"
     " rad  -- radius to aim for.\n"
     " acc  -- computationa accuracy parameter."
    },

    {"drpot", roche_drpot, METH_VARARGS, 
     "dx,dy,dz = drpot(q, r), q = M2/M1. Returns partial derivs of Roche potential at position r"},

    {"drpot1", roche_drpot1, METH_VARARGS, 
     "dx,dy,dz = drpot1(q, spin, r), q = M2/M1. Returns partial derivs of Roche potential for star 1 at position r, spin = spin/orbital"},

    {"face", roche_face, METH_VARARGS,
     "p,d,r,g = face(q, spin, dirn, rref, pref, star=2, acc=1.e-5), returns position and direction of element of specific Roche potential.\n\n"
     " q      -- mass ratio = M2/M1\n"
     " spin   -- ratio spin/orbital frequency\n"
     " dirn   -- direction (a Vec3) to take from centre of mass of star in question.\n"
     " rref   -- reference radius greater than any radius of potential in question.\n"
     " pref   -- the potential to aim for.\n"
     " star   -- 1 or 2 for primary or secondary star."
     " acc    -- accuracy in terms of separation of location.\n"
     "Returns p = position, d = direction perpendicular to face, r = radius from centre of mass, g = gravity."
    },

    {"ref_sphere", roche_ref_sphere, METH_VARARGS, 
     "(rref,pref) = ref_sphere(q, spin, ffac, star=2), returns reference radius and potential needed for face.\n\n"
     " q      -- mass ratio = M2/M1\n"
     " spin   -- ratio spin/orbital frequency\n"
     " ffac   -- linear filling factor of star in question, defined as the radius of the star along the line of\n"
     "           centres towards its companion divided by the Roche lobe radius in that direction. For spin = 1\n"
     "           the latter is simply the distance to the L1 point, but otherwise you need to use modified L1\n"
     "           radii as returned by xl11 or xl12.\n"
     " star   -- 1 or 2 for primary or secondary star."
    },


    {"findi", roche_findi, METH_VARARGS, 
     "findi(q, deltaphi, acc=1.e-4, di=1.e-5), computes inclination for a given mass ratio and phase width"},

    {"findq", roche_findq, METH_VARARGS, 
     "findq(i, deltaphi, acc=1.e-4, dq=1.e-5, qlo=0.005, qhi=2.), computes mass ratio q for a given phase width and inclination"},

    {"findphi", roche_findphi, METH_VARARGS, 
     "findphi(q, i, delta=1.e-6), computes deltaphi for a given mass ratio and inclination"},

    {"fblink", roche_fblink, METH_VARARGS, 
     "fblink(q, i, phi, r, ffac=1., acc=1.e-4, star=2, spin=1), computes whether a point is eclipsed or not"},

    {"ineg", roche_ineg, METH_VARARGS, 
     "in,out = ineg(q, i, x, y, z=0, ffac=1., delta=1.e-7, star=2, spin=1), computes ingress and egress phase of a point"},

    {"lobe1", roche_lobe1, METH_VARARGS, 
     "x,y = lobe1(q, n=200), q = M2/M1. Returns arrays of primary star's Roche lobe."},

    {"lobe2", roche_lobe2, METH_VARARGS, 
     "x,y = lobe2(q, n=200), q = M2/M1. Returns arrays of secondary star's Roche lobe."},

    {"rpot", roche_rpot, METH_VARARGS, 
     "rp = rpot(q, r), q = M2/M1. Returns Roche potential at position r."},

    {"rpot1", roche_rpot1, METH_VARARGS, 
     "rp = rpot1(q, spin, r), q = M2/M1. Returns asynchronous Roche potential for star 1 at position r, spin = spin/orbital"},

    {"rpot2", roche_rpot2, METH_VARARGS, 
     "rp = rpot2(q, spin, r), q = M2/M1. Returns asynchronous Roche potential for star 2 at position r, spin = spin/orbital"},

    {"shadow", roche_shadow, METH_VARARGS, 
     "x,y,s = shadow(q, iangle, phi, n=200, dist=5., acc=1.e-4), q = M2/M1. Returns 2xn array of representing the eclipse shadow region."},

    {"streamr", roche_streamr, METH_VARARGS, 
     "x,y = streamr(q, rad, n=200), returns arrays of the gas stream. q = M2/M1, rad = minimum radius to aim for."},

    {"stream", roche_stream, METH_VARARGS, 
     "x,y = stream(q, step, n=200), returns arrays of the gas stream. q = M2/M1, step=distance between adjacent points."},

    {"strmnx", roche_strmnx, METH_VARARGS, 
     "x,y,vx1,vy1,vx2,vy2 = strmnx(q, n=1, acc=1.e-7), q = M2/M1. Calculates position & velocity of n-th turning point of stream. Two sets of velocities are reported, the first for the pure stream, the second for the disk at that point.\n"},

    {"vlobe1", roche_vlobe1, METH_VARARGS, 
     "(vx,vy) = vlobe1(q, n=200), q = M2/M1. Returns arrays of primary star's Roche lobe in velocity space."},

    {"vlobe2", roche_vlobe2, METH_VARARGS, 
     "vx,vy = vlobe2(q, n=200), q = M2/M1. Returns arrays of secondary star's Roche lobe in velocity space."},

    {"vstream", roche_vstream, METH_VARARGS, 
     "vx,vy = vstream(q, step=0.01, type=1, n=60), q = M2/M1. Returns arrays of positions of the gas stream in velocity space. step is measured as a fraction of the distance to the inner Lagrangian point from the primary star. type=1 is the straight velocity of the gas stream while type=2 is the velocity of the disc along the stream"},

    {"pvstream", roche_pvstream, METH_VARARGS, 
     "x, y, vx, vy, t, jac = pvstream(q, step=0.01, type=1, n=60), q = M2/M1. Returns arrays of positions, velocities time and jacobi constant along THE GAS STREAM. step is measured as a fraction of the distance to the inner Lagrangian point.\n"},

    {"xl1", roche_xl1, METH_VARARGS,
     "xl1(q), q = M2/M1. Calculate the inner Lagrangian point distance."},

    {"xl2", roche_xl2, METH_VARARGS,
     "xl2(q), q = M2/M1. Calculate the L2 point distance."},

    {"xl3", roche_xl3, METH_VARARGS,
     "xl3(q), q = M2/M1. Calculate the L3 point distance."},
    {"xl11", roche_xl11, METH_VARARGS,
     "xl11(q,spin), q = M2/M1, spin = spin/orbital of primary. Calculate the inner Lagrangian point distance with asynchronous primary."},

    {"xl12", roche_xl12, METH_VARARGS,
     "xl12(q,spin), q = M2/M1, spin = spin/orbital of secondary. Calculate the inner Lagrangian point distance with asynchronous secondary."},

    {NULL, NULL, 0, NULL} /* Sentinel */
};

PyMODINIT_FUNC
init_roche(void)
{
    (void) Py_InitModule("_roche", RocheMethods);
    import_array();
}
