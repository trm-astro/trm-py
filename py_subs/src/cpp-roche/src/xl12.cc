/*

!!begin
!!title   xl12
!!author  T.R. Marsh
!!created 30 May 2011
!!root    xl12
!!index   xl12.cc
!!head1   xl12 computes L1 distance from primary accounting for asynchronous rotation of secondary

!!emph{xl12} finds the x-cpt of the L1 Lagrangian point in terms of 
the separation of the two stars which have mass ratio q. Works by 
solving for root of a quintic polynomial, by Newton-Raphson iteration.
L1 is the point in between the two stars, and so will lie between
0 and 1. Part of namespace Roche.

See also: !!ref{xl1.html}{xl1}, !!ref{xl2.html}{xl2}, !!ref{xl3.html}{xl3}, !!ref{xl11.html}{xl11}

!!head2 Function call

double xl12(double q, double spin)

where q = M2/M1 is the mass ratio, spin = ratio of angular/orbital frequency of secondary and 
it returns with the (distance primary --> inner lagrangian)/separation

!!end
*/

#include <cstdlib>
#include <cmath>
#include <iostream>
#include <string>
#include "trm/subs.h"
#include "trm/roche.h"

/**
 * xl12 is a function which returns the distance of the L1 point from
 * star 1 scaled by the orbital separation, allowing for asynchronous rotation of
 * the secondary
 * \param q mass ratio = m2/m1.
 * \param spin spin/orbital frequency ratio of secondary star
 */

double Roche::xl12(double q, double spin){

    const int NMAX   = 1000;
    const double EPS = 1.e-12;

    double ssq = spin*spin;
    double x, xold, f, df, mu;
    double a1, a2, a3, a4, a5, a6, d1, d2, d3, d4, d5;
    int n;

    if(q <= 0.)
        throw Roche_Error("Roche::xl12(double, double): q = " + Subs::str(q) + " <= 0.");

    // Set poly coefficients
    mu = q/(1.+q);
    a1 = -1.+mu;
    d1 = 1.*(a2 =  2.-2.*mu);
    d2 = 2.*(a3 = -ssq+mu);
    d3 = 3.*(a4 =  3.*ssq+2.*mu-2.);
    d4 = 4.*(a5 = 1.-mu-3.*ssq);
    d5 = 5.*(a6 =  ssq);

    n = 0;
    xold = 0.;
    x    = 1./(1.+ q);
    while(n < NMAX && fabs(x-xold) > EPS*fabs(x)){
        xold = x;
        f    = x*(x*(x*(x*(x*a6+a5)+a4)+a3)+a2)+a1;
        df   = x*(x*(x*(x*d5+d4)+d3)+d2)+d1;
        x   -= f/df;
        n++;
    }
    if(n == NMAX)
        throw Roche_Error("Roche::xl12(double,double): exceeded maximum iterations");
    return x;
}
