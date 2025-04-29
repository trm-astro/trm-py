#ifndef DOPPLER_H
#define DOPPLER_H

#include <vector>
#include <string>
#include <complex>
#include <iostream>
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>
#include <fftw3.h>

namespace py = pybind11;

// Typedef and define to interface between C++'s native complex type and 
// fftw_complex
typedef std::complex<double> complex_fft;
#define RECAST reinterpret_cast<fftw_complex*>

// Constants
const double CKMS = 299792.458; // Speed of light in km/s
const double EFAC = 2.354820045; // Ratio FWHM/sigma for a Gaussian

// Structure for image dimensions
struct Nxyz {
    Nxyz(size_t nx, size_t ny, size_t nz = 1) : nx(nx), ny(ny), nz(nz) {}
    size_t ntot() const { return nx * ny * nz; }
    size_t nx, ny, nz;
};

// enum to clarify default selection
enum DefOpt {UNIFORM = 1, GAUSS2D = 2, GAUSS3D = 3};

// Enum for image types
enum Itype {
    PUNIT = 1,
    NUNIT = 2,
    PSINE = 3,
    NSINE = 4,
    PCOSINE = 5,
    NCOSINE = 6,
    PSINE2 = 7,
    NSINE2 = 8,
    PCOSINE2 = 9,
    NCOSINE2 = 10
};


// Namespace for Doppler-related globals
namespace Dopp {
    extern std::vector<Nxyz> nxyz;
    extern std::vector<DefOpt> def;
    extern std::vector<Itype> itype;
    extern std::vector<double> vxy, vz, bias, fwhmxy, fwhmz, squeeze, sqfwhm;
    extern std::vector<std::vector<double>> wavel;
    extern std::vector<std::vector<float>> gamma, scale;
    extern double tzero, period, quad, vfine, sfac;

    extern std::vector<size_t> nwave, nspec;
    extern std::vector<std::vector<double>> time;
    extern std::vector<std::vector<float>> expose;
    extern std::vector<std::vector<int>> nsub;
    extern std::vector<double> fwhm;
    extern double* wave;
}

// Namespace for MEMSYS-related globals and routines
namespace Mem {
    size_t memsize(size_t nipix, size_t ndpix);
    void memcore(size_t mxbuff, size_t nipix, size_t ndpix);
    void opus(const int j, const int k);
    void tropus(const int k, const int j);
    void memprm(int mode, int npar, float caim, float rmax, float acc, float& c, float& test, float& cnew, float& s, float& rnew, float& snew, float& sumf);
}

// Function declarations
bool npix_map(const py::object& Map, size_t& npix);
bool read_map(py::object& map,
              float* images,
              std::vector<Nxyz>& nxyz,
              std::vector<double>& vxy,
              std::vector<double>& vz,
              std::vector<std::vector<double>>& wave,
              std::vector<std::vector<float>>& gamma,
              std::vector<std::vector<float>>& scale,
              std::vector<Itype>& itype,
              std::vector<DefOpt>& def,
              std::vector<double>& bias,
              std::vector<double>& fwhmxy,
              std::vector<double>& fwhmz,
              std::vector<double>& squeeze,
              std::vector<double>& sqfwhm,
              double& tzero,
              double& period,
              double& quad,
              double& vfine,
              double& sfac);
bool update_map(const float* images, py::object& Map);

bool npix_data(const py::object& Data, size_t& npix);
bool read_data(const py::object& Data, float* flux, float* ferr, double* wave,
               std::vector<size_t>& nwave, std::vector<size_t>& nspec,
               std::vector<std::vector<double>>& time,
               std::vector<std::vector<float>>& expose,
               std::vector<std::vector<int>>& nsub,
               std::vector<double>& fwhm);
bool update_data(const float* flux, const py::object& Data);

void op(const float* image, const std::vector<Nxyz>& nxyz,
        const std::vector<double>& vxy, const std::vector<double>& vz,
        const std::vector<std::vector<double>>& wavel,
        const std::vector<std::vector<float>>& gamma,
        const std::vector<std::vector<float>>& scale,
        const std::vector<Itype>& itype,
        double tzero, double period, double quad, double vfine, double sfac,
        float* data, const double* wave,
        const std::vector<size_t>& nwave, const std::vector<size_t>& nspec,
        const std::vector<std::vector<double>>& time,
        const std::vector<std::vector<float>>& expose,
        const std::vector<std::vector<int>>& nsub,
        const std::vector<double>& fwhm);

void tr(float* image, const std::vector<Nxyz>& nxyz,
        const std::vector<double>& vxy, const std::vector<double>& vz,
        const std::vector<std::vector<double>>& wavel,
        const std::vector<std::vector<float>>& gamma,
        const std::vector<std::vector<float>>& scale,
        const std::vector<Itype>& itype,
        double tzero, double period, double quad, double vfine, double sfac,
        const float* data, const double* wave,
        const std::vector<size_t>& nwave, const std::vector<size_t>& nspec,
        const std::vector<std::vector<double>>& time,
        const std::vector<std::vector<float>>& expose,
        const std::vector<std::vector<int>>& nsub,
        const std::vector<double>& fwhm);

void gaussdef(const float* input, const Nxyz& nxyz, double fwhmx,
              double fwhmy, double fwhmz, float* output);

// Python-exposed functions
py::none doppler_comdat(py::object& map, py::object& data);
py::none doppler_comdef(py::object& map);
py::none doppler_datcom(py::object& data, py::object& map);
py::none doppler_memit(py::object& map, py::object& data, int niter, float caim, float tlim, float rmax);

#endif // DOPPLER_H