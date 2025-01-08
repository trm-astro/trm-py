/*

!!begin
!!title   xl3
!!author  T.R. Marsh
!!created 7 Dec 2000
!!root    xl3
!!index   xl3.cc
!!head1   xl3 computes L3 distance from primary

!!emph{xl3} finds the x-cpt of the L3 Lagrangian point in terms of 
the separation of the two stars which have mass ratio q. Works by 
solving for root of a quintic polynomial, by Newton-Raphson iteration.
I define L3 as the point on the side of the primary opposite the secondary
and therefore xl3 < 0.

See also: !!ref{xl1.html}{xl1}, !!ref{xl2.html}{xl2}

!!head2 Function call

double xl3(double q)

where q = M2/M1 is the mass ratio and it returns with the 
(distance primary --> inner lagrangian)/separation

!!end

*/

#include <cstdlib>
#include <cmath>
#include <iostream>
#include <string>
#include "trm/subs.h"
#include "trm/roche.h"

using Subs::operator+;

double Roche::xl3(double q){

  const int NMAX = 1000;
  const double EPS = 1.e-12;

  double x, mu, xold, f, df;
  double a1, a2, a3, a4, a5, a6, d1, d2, d3, d4, d5;
  int n;

  if(q <= 0.)
      throw Roche_Error("Roche::xl3(double): q = " + Subs::str(q) + " <= 0.");

  // Set poly coefficients
  mu = q/(1.+q);
  a1 = 1.-mu;
  d1 = 1.*(a2 = -2.+2.*mu);
  d2 = 2.*(a3 =  1.-mu);
  d3 = 3.*(a4 =  1.+2.*mu);
  d4 = 4.*(a5 = -2.-mu);
  d5 = 5.*(a6 =  1.);

  n    = 0;
  xold = 0.;
  x    = -1.;
  while(n < NMAX && fabs(x-xold) > EPS*fabs(x)){
      xold = x;
      f    = x*(x*(x*(x*(x*a6+a5)+a4)+a3)+a2)+a1;
      df   = x*(x*(x*(x*d5+d4)+d3)+d2)+d1;
      x   -= f/df;
      n++;
  }
  if(n == NMAX)
      throw Roche_Error("Roche::xl3(double): exceeded maximum iterations");
  return x;
}
