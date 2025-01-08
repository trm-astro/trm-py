#include <cstdlib>
#include <iostream>
#include "trm/constants.h"
#include "trm/subs.h"
#include "trm/binary.h"

/** jdotgr computes the rate of change of 
 * orbital angular momentum divided by the orbital
 * angular momentum (Jdot/J) for two stars in a 
 * circular orbit. 
 * \param m1 mass of star 1 (solar masses)
 * \param m2 mass of star 2 (solar masses)
 * \param a  orbital separation (solar radii)
 * \return The value is returned has units of s^-1, and is negative.
 */

double Binary::jdotgr(double m1, double m2, double a){
  return (-32./5.*pow(Constants::G*Constants::MSUN/Constants::C/Constants::RSUN,3)/Subs::sqr(Constants::C)/
	  Constants::RSUN*m1*m2*(m1+m2)/pow(a,4));
}


