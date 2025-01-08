#include <cstdlib>
#include <string>
#include <cmath>
#include "trm/subs.h"

/**
 * gcf returns incomplete gamma function Q(a,x) evaluated by its
 * continued fraction representation as gammcf and ln Gamma(a) as gln 
 * \param gammcf returned continued fraction value
 * \param a first parameter of incomplete gamma function
 * \param x second parameter of incomplete gamma function
 * \param gln ln(Gamma(x)), returned
 * \exception Throws Subs::Subs_Error exceptions.
 * \sa   gser(double&, double, double, double&), gammp(double, double),
 * gammq(double, double), gammln(double). 
 */

void Subs::gcf(double &gammcf, double a, double x, double &gln){
  const int ITMAX=100;
  const double EPS=3.e-7;
  const double FPMIN=1.e-30;
  int i;
  double an,b,c,d,del,h;
  
  gln = gammln(a);
  b = x+1.0-a;
  c=1.0/FPMIN;
  h=d=1.0/b;
  for(i=1;i<=ITMAX;i++){
    an = -i*(i-a);
    b += 2.0;
    d = an*d+b;
    if(fabs(d) < FPMIN) d=FPMIN;
    c=b+an/c;
    if(fabs(c) < FPMIN) c=FPMIN;
    d = 1.0/d;
    del = d*c;
    h *= del;
    if(fabs(del-1.0) < EPS) break;
  }
  if(i > ITMAX)
    throw Subs_Error("a too large, ITMAX too small in gcf");

  gammcf = exp(-x+a*log(x)-gln)*h;
}
