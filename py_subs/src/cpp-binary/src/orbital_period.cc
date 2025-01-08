#include <cstdlib>
#include <cmath>
#include "trm/subs.h"
#include "trm/constants.h"
#include "trm/binary.h"

/*
 * orbital_period computes the ornital period of a binary star.
 * \param m1 mass of star 1 (solar masses)
 * \param m2 mass of star 2 (solar masses)
 * \param a orbital separation (solar radii)
 * \return Orbital period in seconds.
 */

double Binary::orbital_period(double m1, double m2, double a){
  return(Constants::TWOPI*pow(a*Constants::RSUN,1.5)/sqrt(Constants::G*Constants::MSUN*(m1+m2)));
}

