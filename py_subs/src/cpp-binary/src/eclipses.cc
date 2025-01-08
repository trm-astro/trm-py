#include <cstdlib>
#include <math.h>
#include "trm/subs.h"
#include "trm/constants.h"
#include "trm/binary.h"

/**
 * This program computes whether a binary of specific masses,
 * radii, period and inclination will eclipse. It assumes that the
 * two stars are spherical. If eclipses do occur, the function works
 * out their length too.
 * \param m1 mass of star 1 (solar masses)
 * \param m2 mass of star 2 (solar masses)
 * \param r1 radius of star 1 (solar radii)
 * \param r2 radius of star 2 (solar radii)
 * \param period orbital priod (seconds)
 * \param iangle inclination angle (deg)
 * \param width eclipse width, from first to last contact in seconds (only if the star eclipses).
 * \return true or false according to whether the binary eclipses or not.
 */

bool Binary::eclipses(double m1, double m2, double r1, double r2, double period, double iangle, double& width){

  // compute orbital separation.

  double a = orbital_separation(m1,m2,period);

  double tot = (r1+r2)/a;
  if(cos(Constants::TWOPI*iangle/360.) < tot){
    double cphi = sqrt(1.-tot*tot)/sin(Constants::TWOPI*iangle/360.);
    width = 2.*period*acos(cphi)/Constants::TWOPI;
    return true;
  }else{
    return false;
  }
}
	
