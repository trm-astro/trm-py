#include <cmath>
#include <string>
#include "trm/subs.h"
#include "trm/binary.h"

/**
 * Implements Nauenberg's 1972 analytic approx to degenerate M-R 
 * relation. Given m in solar units comes back with r in solar units.
 * Not good for low masses, \sa Binary::mr_wd_eggleton(double) for
 * an alternative.
 * \param  m mass of white dwarf (solar masses)
 * \return Radius of white dwarf (solar radii) 
 */
double Binary::mr_wd_nauenberg(double m){
  if(m <= 0. || m >= 1.433)
    throw Binary_Error("Binary::mr_wd_nauenberg(double): m = " + Subs::str(m) + " out of range.");
  double fac;
  fac = pow(m/1.433,2./3.);
  return 0.0112*sqrt(1./fac-fac);
}

