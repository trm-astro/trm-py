#include <cstdlib>
#include <string>
#include <cmath>
#include "trm/subs.h"

/**
 * gser returns incomplete gamma function gamma(a,x) evaluated by its
 * series representation as gamser and ln Gamma(a) as gln 
 * \param gamser returned series sum
 * \param a first parameter of incomplete gamma function
 * \param x second parameter of incomplete gamma function
 * \param gln ln(Gamma(x)), returned
 * \exception Throws Subs::Subs_Error exceptions.
 * \sa   gcf(double&, double, double, double&), gammp(double, double),
 * gammq(double, double), gammln(double). 
 */


void Subs::gser(double &gamser, double a, double x, double &gln){
  const int ITMAX=100;
  const double EPS=3.e-7;
  int n;
  double sum, del, ap;
  
  gln = gammln(a);
  if(x <= 0.0){
    if(x < 0.0)
      throw Subs_Error("x less than 0 in gser");
    gamser = 0.0;
    return;
  }else{
    ap = a;
    del = sum = 1.0/a;
    for(n=0;n<ITMAX;n++){
      ++ap;
      del *= x/ap;
      sum += del;
      if(fabs(del) < fabs(sum)*EPS){
	gamser = sum*exp(-x+a*log(x)-gln);
	return;
      }
    }
    throw Subs_Error("a too large, ITMAX too small in routine gser");
  }
}
