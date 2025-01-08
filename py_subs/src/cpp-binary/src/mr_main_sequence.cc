#include <cmath>
#include "trm/binary.h"

/**
 * Returns main-sequence radius given a mass. Completely crude
 * at the moment just uses R = M in solar units. Needs to be updated
 * for anything critical.
 * \param m mass of star (solar masses)
 * \return Radius of star (solar radii)
 */
double Binary::mr_main_sequence(double m){
  return m;
}

