#include <cstdlib>
#include <cmath>
#include <string>
#include "trm/subs.h"

/**
 * gammp returns incomplete gamma function P(a,x)
 * \param a first parameter of incomplete gamma function
 * \param x second parameter of incomplete gamma function
 * \return P(a,x)
 * \exception Throws Subs::Subs_Error exceptions.
 * \sa   gcf(double&, double, double, double&), gser(double&, double, double, double&),
 * gammq(double, double), gammln(double). 
 */

double Subs::gammp(double a, double x){
  double gamser, gammcf, gln;

  if(x < 0.0 || a <= 0.0)
    throw Subs_Error("Invalid arguments in routine gammp, x < 0 or a <= 0");

  if(x < (a+1.0)){
    gser(gamser,a,x,gln);
    return gamser;
  }else{
    gcf(gammcf,a,x,gln);
    return 1.0-gammcf;
  }
}


