#include <cmath>
#include <string>
#include "trm/subs.h"
#include "trm/binary.h"

/**
 * If the moment of inertia is given by k M R^2, this program computes 'k' for 
 * a Chandrasekhar white dwarf. Based upon computations I did and a slight modification
 * of a fitting formula from Gijs Nelemans, which is good to < 1.8% 
 * \param m the mass of the white dwarf in solar masses. Must be from 0 to 1.44885 (not checked)
 * \return Returns k
  */
double Binary::gyrate(double m){
  if(m <= 0. || m >= 1.44885)
    throw Binary_Error("Binary::gyrate(double): m = " + Subs::str(m) + " out of range.");
  return 0.1939*pow(1.44885 - m, 0.1917);
}

/**
 * If the moment of inertia is given by k M R^2, this program computes ' d log(k) / d log(m)' for 
 * a Chandrasekhar white dwarf. Based upon computations I did and a slight modification
 * of a fitting formula from Gijs Nelemans, which is good to < 1.8% 
 * \param m the mass of the white dwarf in solar masses. Must be from 0 to 1.44885 (not checked)
 * \return Returns d log(k) / d log(m)
  */
double Binary::zeta_gyrate(double m){
  if(m <= 0. || m >= 1.44885)
    throw Binary_Error("Binary::gyrate(double): m = " + Subs::str(m) + " out of range.");
  return -0.1917*m/(1.44885 - m);
}


