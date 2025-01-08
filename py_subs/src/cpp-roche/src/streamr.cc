#include <cstdlib>
#include <cmath>
#include <iostream>
#include "trm/roche.h"
#include "trm/vec3.h"
#include "trm/subs.h"

/**
 * streamr works by integrating the equations of motion for the Roche
 * potential using Burlisch-Stoer integration. It stops when the stream
 * reaches a target radius or a minimum radius, whichever is the larger.
 *
 * \param q mass ratio = M2/M1. Stream flows from star 2 to 1.
 * \param rad Radius to aim for. If this is less than the minimum, the stream
 * will stop at the minimum
 * \param x array of x values returned.
 * \param y array of y values returned.
 * \param n number of points to compute.
 */
void Roche::streamr(double q, double rad, float x[], float y[], int n){

  if(n < 2)
    throw Roche_Error("void Roche::streamr(double, double, float*, float*, int): need at least 2 points");

  const double EPS = 1.e-8;
  double rl1  = xl1(q);
  Subs::Vec3 r, v, rs, vs;
  strinit(q,r,v);
  rs = r;
  vs = v;
  strmnx(q,r,v,EPS);
  double rmin = r.length();
  rmin  = rmin > rad ? rmin : rad;

  r = rs;
  v = vs;
  x[0] = r.x();
  y[0] = r.y();
  double rnext;
  for(int i=1; i<n; i++){
    rnext = rl1 + (rmin-rl1)*i/(n-1);
    stradv(q, r, v, rnext, 1.e-6, 1.e-4);
    x[i] = r.x();
    y[i] = r.y();
  }
}









