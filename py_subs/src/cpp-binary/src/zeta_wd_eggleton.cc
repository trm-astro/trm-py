#include <cmath>
#include <string>
#include "trm/subs.h"
#include "trm/binary.h"

/**
 * The logarithmic derivative of radius wrt mass is an important parameter of a star
 * in determining the statbility and rate of mass transfer. This function calculates the
 * value for a white dwarf according to Eggleton's 1986 analytic approximation to the M-R relation
 * of degenerate stars, as quoted in Verbunt and Rappaport (1988).
 * \param m the mass of the white dwarf in solar masses
 * \return Returns d log(R) / d log(M)
 */
double Binary::zeta_wd_eggleton(double m){
  if(m <= 0. || m >= 1.44)
    throw Binary_Error("Binary::mr_wd_eggleton(double): m = " + Subs::str(m) + " out of range.");

  double fac1 = pow(m/1.44,4./3.);
  double fac2 = 0.00057/m;
  double fac3 = pow(fac2,2./3.);

  return -(1.+fac1)/(1.-fac1)/3.+ 2.*(7.*fac3/3.+fac2)/(1.+3.5*fac3+fac2)/3.;
}

/**
 * This function computes the derivative wrt mass of the logarithmic derivative of radius wrt mass
 * for a white dwarf according to Eggleton's 1986 analytic approximation to the M-R relation
 * of degenerate stars, as quoted in Verbunt and Rappaport (1988). Successfully tested against
 * finite difference version.
 * \param m the mass of the white dwarf in solar masses
 * \return Returns d/d M of d log(R) / d log(M)
 */
double Binary::dzetadm_wd_eggleton(double m){
  double fac1 = m/1.44;
  double fac2 = pow(fac1,1./3.);
  double fac3 = fac1*fac2;
  double fac4 = m/0.00057;
  double fac5 = pow(fac4,2./3.);

  return -8./9.*fac2/pow(1.-fac3,2)/1.44 - 2./3.*((14./9.+7./18/fac4)/fac5+1./fac4)/fac4/
    pow(1.+7./2./fac5+1./fac4,2)/0.00057;
}

