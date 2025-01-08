#include <cmath>
#include <string>
#include "trm/subs.h"
#include "trm/binary.h"

/**
 * Implements Eggleton's 1986 analytic approx to degenerate M-R 
 * relation. Given m in solar units comes back with r in solar units.
 * The formula is quoted in Verbunt and Rappaport (1988) as a
 * "private communication". It is close to Nauenberg at high masses
 * but is better at low masses.
 * \param  m mass of white dwarf (solar masses)
 * \return Radius of white dwarf (solar radii) 
 */
double Binary::mr_wd_eggleton(double m){
  if(m <= 0. || m >= 1.44)
    throw Binary_Error("Binary::mr_wd_eggleton(double): m = " + Subs::str(m) + " out of range.");

  double fac1 = pow(m/1.44,2./3.);
  double fac2 = 0.00057/m;
  return 0.0114*sqrt(1./fac1-fac1)*
    pow(1.+3.5*pow(fac2,2./3.)+fac2,-2./3.);
}

