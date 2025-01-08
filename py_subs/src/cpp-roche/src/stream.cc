#include <cstdlib>
#include <cmath>
#include <iostream>
#include "trm/subs.h"
#include "trm/vec3.h"
#include "trm/roche.h"

/**
 * stream works by integrating the equations of motion for the Roche
 * potential using Burlisch-Stoer integration. Every time the distance
 * from the last point exceeds step, it interpolates and stores a new 
 * point. This allows one not to spend loads of points on regions where
 * nothing is happening.
 *
 * \param q    mass ratio = M2/M1. Stream flows from star 2 to 1.
 * \param step step between points (units of separation).
 * \param x    array of x values returned.
 * \param y    array of y values returned.
 * \param n    number of points to compute.
 */

void Roche::stream(double q, double step, float x[], float y[], int n){

  const double EPS = 1.e-8;
  int lp;
  double rl1, smax, vel, ttry, tdid, time, tnext, dist, frac;
  Subs::Vec3 r, v;

  if(n < 2)
    throw Roche::Roche_Error("Need at least 2 points in stream");

  // Initialise stream

  rl1 = xl1(q);
  strinit(q,r,v);
  time = 0.0;  

  // Store L1 as first point

  x[0] = rl1;
  y[0] = 0.0;

  lp = 0;

  // Store interpolation between L1 and initial point if
  // step has been set small enough 

  dist = sqrt(Subs::sqr(r.x()-rl1)+Subs::sqr(r.y()));
  if(dist > step){
    x[1] = rl1 + (r.x()-rl1)*(frac = step/dist);
    y[1] = r.y()*frac;
    lp++;
  }

  ttry = 1.e-3;
  smax = std::min(1.e-3,step/2.);

  while(lp < n-1){
    gsint(q, r, v, ttry, tdid, tnext, time, EPS);
    dist = sqrt(Subs::sqr(r.x()-x[lp])+Subs::sqr(r.y()-y[lp]));
    if(dist > step){
      x[lp+1] = x[lp] + (r.x()-x[lp])*(frac = step/dist);
      y[lp+1] = y[lp] + (r.y()-y[lp])*frac;
      lp++;
    }
    vel  = sqrt(Subs::sqr(v.x()) + Subs::sqr(v.y()));
    ttry = std::min(smax/vel, tnext);
  }
}









