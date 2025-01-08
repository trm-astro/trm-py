#include <cstdlib>
#include <cmath>
#include "trm/binary.h"

/**
 * rcirc computes the relative size of the circularisation radius according to Verbunt and
 * Rappaport (1988) as a function of mass ratio q = m2/m1,
 * where mass transfer is from star 2 to star 1.
 * \param q mass ratio = M2/M1
 * return The circularisation radius divided by the binary separation.
 */

double Binary::rcirc(double q){
  double f = log10(q);
  return (0.0883+f*(-0.04858+f*(0.11489+f*0.02020475)));
}

/**
 * drcircdq computes the derivative wrt q of the relative size of relative size of the circularisation 
 * radius according to the formula of Verbunt and Rappaport (1988) as a function of mass ratio q = m2/m1,
 * where mass transfer is from star 2 to star 1.
 * \param q mass ratio = M2/M1
 * return The derivative wrt q of the circularisation radius divided by the binary separation.
 */

double Binary::drcircdq(double q){
  double f = log10(q);
  return ((-0.04858+f*(2.*0.11489+f*3.*0.02020475))/q/log(10.));
}


