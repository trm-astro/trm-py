#include <cstdlib>
#include <cmath>
#include <iostream>
#include "trm/subs.h"
#include "trm/vec3.h"
#include "trm/roche.h"

/**
 * vstream computes the path of the gas stream in a binary in velocity space.
 * There are a few different options for the type of stream velocities produced.
 * vstream works by integrating the equations of motion for the Roche
 * potential using Burlisch-Stoer integration. Every time the speed
 * changes by step, it interpolates and stores a new point.
 *
 * The velocities are inertial frame velocities for comparison with Doppler maps.
 * 
 * \param q mass ratio = M2/M1. Stream flows from star 2 to 1.
 * \param step step between points (units of K1+K2).
 * \param vx array of x velocities returned.
 * \param vy array of y velocities returned.
 * \param vy array of y velocities returned.
 * \param rad array of radii equivalent to velocities (units of a)
 * \param n number of points to compute.
 * \param type type of velocity, see !!ref{vtrans.html}{vtrans} for supported types.
 */

void Roche::vstream(double q, double step, float vx[], float vy[], float rad[], int n, int type){

  const double EPS =  1.e-8;
  int lp;
  double rl1, smax, acc, ttry, tdid, time, tnext, dvel, frac;
  double tvx, tvy;
  Subs::Vec3 r, v, ar;

  if(n < 2) 
    throw Roche::Roche_Error("Need at least 2 points in vstream\n");

  // Initialise stream

  rl1 = xl1(q);
  strinit(q,r,v);
  time = 0.0;  

  /* Store L1 as first point */
  vtrans(q, type, rl1, 0., 0., 0., tvx, tvy);
  vx[0] = tvx;
  vy[0] = tvy;
  lp = 0;

  /* Store interpolation between L1 and initial point if
     step has been set small enough */

  vtrans(q, type, r.x(), r.y(), v.x(), v.y(), tvx, tvy);
  dvel = sqrt(Subs::sqr(tvx-vx[0])+Subs::sqr(tvy-vy[0]));
  if(dvel > step){
    vx[1]  = vx[0] + (tvx-vx[0])*(frac = step/dvel);
    vy[1]  = vy[0] + (tvy-vy[0])*frac;
    rad[1] = sqrt(Subs::sqr(r.x())+Subs::sqr(r.y()));
    lp++;
  }

  ttry = 1.e-3;
  smax = std::min(1.e-3,step/2.);

  while(lp < n-1){
    gsint(q, r, v, ttry, tdid, tnext, time, EPS);
    vtrans(q, type, r.x(), r.y(), v.x(), v.y(), tvx, tvy);
    dvel = sqrt(Subs::sqr(tvx-vx[lp])+Subs::sqr(tvy-vy[lp]));
    if(dvel > step){
      vx[lp+1]  = vx[lp] + (tvx-vx[lp])*(frac = step/dvel);
      vy[lp+1]  = vy[lp] + (tvy-vy[lp])*frac;
      rad[lp+1] = sqrt(Subs::sqr(r.x())+Subs::sqr(r.y()));
      lp++;
    }
    ar   = rocacc(q, r, v);
    acc  = sqrt(Subs::sqr(ar.x()) + Subs::sqr(ar.y()));
    ttry = std::min(smax/acc, tnext);
  }
}

