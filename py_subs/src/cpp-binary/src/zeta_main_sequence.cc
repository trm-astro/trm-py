#include <cmath>
#include "trm/binary.h"

/**
 * This function returns the adiabatic logarithmic derivative of radius wrt mass for a
 * main-sequence star corresponding to the mr_main_sequence function which gives the
 * radius. Juts returns -1./3. at the moment (crude approx for v.low mass)
 * \param m the mass of the star in solar masses
 * \return Returns d log(R) / d log(M)
 */
double Binary::zeta_main_sequence(double m){
  return -1./3.;
}

