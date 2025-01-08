#include <cstdlib>
#include <cmath>
#include "trm/subs.h"
#include "trm/constants.h"
#include "trm/binary.h"
//! \file

/**
 * orbital_separation computes the orbital separation of a binary star.
 * \param m1 mass of star 1, solar
 * \param m2 mass of star 2, solar
 * \param period orbital period,  (sec)
 * \return The orbital separation in solar radii
 */
double Binary::orbital_separation(double m1, double m2, double period){
  return(
	 pow(Constants::G*Constants::MSUN*(m1+m2)*
	     Subs::sqr(period/Constants::TWOPI),1./3.)/Constants::RSUN
	 );
}

