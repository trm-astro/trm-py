#include <stdlib.h>
#include <math.h>
#include "trm/subs.h"
#include "trm/constants.h"
#include "trm/binary.h"

/**
 * orbital_omega computes the orbital angular
 * velocity of a binary star in a circular orbit.
 * \param m1 mass of star 1 (solar masses)
 * \param m2 mass of star 2 (solar masses)}
 * \param a orbital separation (solar radii)
 * \return Orbital angular velocity in radians per second
 */
double Binary::orbital_omega(double m1, double m2, double a){
  return(sqrt(Constants::G*Constants::MSUN*(m1+m2)/(Constants::RSUN*a))/(Constants::RSUN*a));
}

