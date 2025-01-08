#include <cmath>
#include "trm/binary.h"

/**
 * The logarithmic derivative of radius wrt mass is an important parameter of a star
 * in determining the statbility and rate of mass transfer. This function calculates the
 * value for a white dwarf according to Nauenberg's 1972 approximation to the M-R relation
 * of degenerate stars.
 * \param m the mass of the white dwarf in solar masses
 * \return Returns d log(R) / d log(M)
 */
double Binary::zeta_wd_nauenberg(double m){
  double fac;
  fac = pow(m/1.433,4./3.);
  return -(1.+fac)/(1.-fac)/3.;
}

