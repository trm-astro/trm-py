/*
!!begin
!!title xl2
!!author T.R. Marsh
!!date 7 Dec 2000
!!root  xl2
!!index xl2.cc
!!head1 xl2 computes L2 distance from primary

xl2 finds the x-cpt of the L2 Lagrangian point in terms of 
the separation of the two stars which have mass ratio q. Works by 
solving for root of a quintic polynomial, by Newton-Raphson iteration.
I define L2 as the point on the side of the primay opposite the secondary,
and therefore xl2 > 1.

See also: !!ref{xl1.html}{xl2}, !!ref{xl3.html}{xl3}

!!head2 Function call

double xl2(double q)

where q = M2/M1 is the mass ratio and it returns with the 
(distance primary --> L2)/separation

!!end
*/

#include <cstdlib>
#include <cmath>
#include <iostream>
#include <string>
#include "trm/subs.h"
#include "trm/roche.h"

double Roche::xl2(double q){

  const int NMAX = 1000;
  const double EPS = 1.e-12;

  double x, xold, f, df, mu;
  double a1, a2, a3, a4, a5, a6, d1, d2, d3, d4, d5;
  int n;

  if(q <= 0.)
      throw Roche_Error("Roche::xl2(double): q = " + Subs::str(q) + " <= 0.");

  // Set poly coefficients
  mu = q/(1.+q);
  a1 = -1.+mu;
  d1 = 1.*(a2 =  2.-2.*mu);
  d2 = 2.*(a3 = -1.-mu);
  d3 = 3.*(a4 =  1.+2.*mu);
  d4 = 4.*(a5 = -2.-mu);
  d5 = 5.*(a6 =  1.);

  n    = 0;
  xold = 0.;
  x    = 1.5;
  while(n < NMAX && fabs(x-xold) > EPS*fabs(x)){
      xold = x;
      f    = x*(x*(x*(x*(x*a6+a5)+a4)+a3)+a2)+a1;
      df   = x*(x*(x*(x*d5+d4)+d3)+d2)+d1;
      x   -= f/df;
      n++;
  }
  if(n == NMAX)
      throw Roche_Error("Roche::xl2(double): exceeded maximum iterations");
  return x;
}
