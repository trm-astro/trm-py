#include <cstdlib>
#include <cmath>
#include <iostream>
#include "trm/subs.h"
#include "trm/roche.h"

/** 
 * lobe1 returns arrays x and y for plotting an equatorial section
 * of the Roche lobe of the primary star in a binary of mass ratio q = M2/M1.
 * The arrays start and end at the inner Lagrangian point and march around 
 * uniformly in azimuth looking from the centre of mass of the primary star.
 * n is the number of points and must be at least 3. 
 *
 * \param q  mass ratio = M2/M1
 * \param x  array of x values
 * \param y  array of y values
 * \param n  number of x and y values
 */
void Roche::lobe1(double q, float x[], float y[], int n){

  // Accuracy of location of surface in terms of binary separation
  const double FRAC = 1.e-6;

  Subs::Vec3 p;
  double rl1, rad, theta, lower, upper;

  if(q <= 0.)
    throw Roche_Error("Roche::lobe1(double, float*, float*, int): q = " + Subs::str(q) + " <= 0.");

  if(n < 3) 
    throw Roche_Error("Roche::lobe1(double, float*, float*, int): n = " + Subs::str(n) + " < 3.");

  // Compute L1 point and critical potential there.

  rl1   = xl1(q);
  p.y()   = p.z() = 0.;
  p.x()   = rl1;
  double cpot  = rpot(q,p);
  upper = rl1;
  lower = upper/4.;

  // Now compute Roche lobe in regular steps of angle looking
  // from centre of Roche lobe

  for(int i=0;i<n;i++){

    // L1 point is a special case because derivative becomes zero there.
    // lambda is set so that after i=0, there is a decent starting 
    // multiplier

    if(i == 0 || i == n-1){
      x[i]   = (float)rl1;
      y[i]   = 0.;
    }else{
      theta = 2*M_PI*(double)i/(n-1);
      double dx =  cos(theta);
      double dy =  sin(theta);


      // Locate critical surface using rtsafe.
      // Based on assuming that rl1 is maximum distance
      // from centre of mass and that at no point is the
      // surface closer than 1/4 of this.

      rad   = Subs::rtsafe(Lfunc1(dx,dy,q,cpot), lower, upper, FRAC);

      x[i] = rad*dx;
      y[i] = rad*dy;
    }
  }
}



