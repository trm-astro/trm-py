#include <cstdlib>
#include <cmath>
#include <string>
#include "trm/subs.h"

/**
 * gammq returns incomplete gamma function Q(a,x)
 * \param a first parameter of incomplete gamma function
 * \param x second parameter of incomplete gamma function
 * \return Q(a,x)
 * \exception Throws Subs::Subs_Error exceptions.
 * \sa   gcf(float&, float, float, float&), gser(float&, float, float, float&),
 * gammp(float, float), gammln(float). 
 */


// standard incomplete gamma function on two doubles
double Subs::gammq(double a, double x){
  double gamser, gammcf, gln;

  if(x < 0.0 || a <= 0.0)
    throw Subs_Error("Invalid arguments in routine gammq, x < 0 or a <= 0");

  if(x < (a+1.0)){
    gser(gamser,a,x,gln);
    return 1.0-gamser;
  }else{
    gcf(gammcf,a,x,gln);
    return gammcf;
  }
}

// vectorized version of gammq on a as a double and x as a vector(double)
void Subs::gammq(double a, double* x, double* out, int n){
  for(int i = 0; i < n; i++){
    out[i] = gammq(a,x[i]);
  }
}