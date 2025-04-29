#ifndef TRM_CONSTANTS
#define TRM_CONSTANTS

#include <cmath>

//! Namespace of useful constants

namespace Constants {

  //! Offset between JD and MJD
  const double MJD2JD = 2400000.5;

  //! Speed of light, MKS, (exact value)
  const double C     = 2.99792458e8;  

  //! Gravitational constant, MKS
  const double G     = 6.673e-11;

  //! Planck's constant, MKS
  const double H     = 6.6262e-34;

  //! Boltzmann's constant, MKS
  const double K     = 1.3806e-23;  

  //! Charge on the electron (magnitude, C)
  const double E     = 1.602176565e-19;

  //! Mass of the electron, kg
  const double ME    = 9.10956e-31; 

  //! Mass of the proton, kg
  const double MP    = 1.67e-27;  

  //! Stefan-Boltzmann constant, MKS
  const double SIGMA = 5.66956e-8;

  //! Thomson cross-section, MKS
  const double SIGMAT = 6.65e-29;
  
  //! Astronomical unit, metres
  const double AU    = 1.49597870691e11;

  //! Solar luminosity, Watts
  const float  LSUN  = 3.826e26;

  //! Solar mass, kg    
  const float  MSUN  = 1.989e30;

  //! Gravitational parameter of the Sun, SI (m^3 s^-2)
  const double GMSUN  = 1.32712442099e20;

  //! Gravitational parameter of the Sun, AU^3 YR^-2
  const double GMSUNA  = 39.476927033270655;

  //! Gauss' gravitational constant sqrt(G*MSUN), AU^(3/2) day^-1
  const double KGAUSS = 0.01720209895;     

  //! G*MSUN, AU^3 day^-2 (Gauss**2)
  const double GMGAUSS = KGAUSS*KGAUSS;     

  //! Absolute visual magnitude of the Sun
  const float  MVSUN = 4.75;       

  //! Parsec, metres
  const double PC    = 3.085678e16;

  //! Solar radius, metres
  const float  RSUN  = 6.9599e8;

  //! Effective temperature of the Sun, Kelvin
  const float  TSUN  = 5700.;  
  
  //! Number of seconds in a day
  const float  DAY   = 86400.;

  //! Length of Julian year in seconds 
  const double YEAR  = 365.25*DAY;

  //! Integer number of seconds in a day
  const int    IDAY  = 86400;

  //! Number of seconds in an hour
  const float  HOUR  = 3600.;     

  //! Number of seconds in a minute
  const float  MINUTE  = 60.;     
  
  //! Pi
  const double PI    = 3.14159265358979323846264;  

  //! 2*Pi
  const double TWOPI = 2.*PI; 

  //! Ratio FWHM/sigma for a gaussian 
  const double EFAC  = sqrt(8.*log(2.)); 

  //! Wavelength of Halpha, Angstroms
  const double HALPHA  = 6562.76; 

  //! Wavelength of Hbeta, Angstroms
  const double HBETA   = 4861.327; 

  //! Wavelength of Hgamma, Angstroms
  const double HGAMMA   = 4340.465; 

  //! Wavelength of Hgamma, Angstroms
  const double HDELTA   = 4340.465;
 
};

#endif







