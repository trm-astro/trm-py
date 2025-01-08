#include <cstdlib>
#include "trm/constants.h"
#include "trm/subs.h"
#include "trm/binary.h"

/**
 * orbital_ang_mom computes the orbital angular momentum
 * of two stars in a circular orbit..
 * \param m1 mass of star 1 (solar masses)
 * \param m2 mass of star 2 (solar masses)
 * \param a Orbital separation (solar radii)
 * \return Orbital angular momentum in units of Msun Rsun**2 s**-1.
 */

double Binary::orbital_ang_mom(double m1, double m2, double a){
  return (m1*m2*sqrt(Constants::G*Constants::MSUN*a/(m1+m2)/pow(Constants::RSUN,3)));
}


