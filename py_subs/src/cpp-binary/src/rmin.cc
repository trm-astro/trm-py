#include <cstdlib>
#include <cmath>
#include "trm/subs.h"
#include "trm/binary.h"

/**
 * This omputes the minimum radius reached by the gas stream
 * on its first orbit of the accretor, useful for deciding whether or not
 * it will hit the accretor (star 1). This is basd upon the fit of
 * Nelemans et al (2001).
 * \param q mass ratio = M2/M1
 * \return Minimum radius in units of the orbital separation
 */
double Binary::rmin(double q){
  double f = log10(q);
  return (0.04948+f*(-0.03815+f*(0.04752-f*0.006973)));
}


