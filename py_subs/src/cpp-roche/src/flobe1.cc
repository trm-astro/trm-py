#include <cstdlib>
#include <cmath>
#include <iostream>
#include "trm/subs.h"
#include "trm/vec3.h"
#include "trm/roche.h"

/**
 * flobe1 returns arrays x and y for plotting a Roche equipotential
 * centred on the primary star. This works for cases other than equal to the 
 * critical Roche potential. Use !!ref{lobe1.html}{lobe1} in this case.
 *
 * \param q   mass ratio = M2/M1
 * \param x   array of x values
 * \param y   array of y values
 * \param n   number of x and y values
 * \param pot the Roche potential to locate
 */

void Roche::flobe1(double q, float x[], float y[], int n, double pot){

  // Accuracy of location of surface in terms of binary separation

  const double FRAC = 1.e-6;

  Subs::Vec3 p;
  double rl1, rad, theta, lower, upper;

  if(q <= 0.)
    throw Roche_Error("Roche::flobe1(double, float*, float*, int): q = " + Subs::str(q) + " <= 0.");

  if(n < 3) 
    throw Roche_Error("Roche::flobe1(double, float*, float*, int): n = " + Subs::str(n) + " < 3.");


  // Compute L1 point and critical potential there.

  rl1   = xl1(q);
  p.y() = p.z() = 0.;
  p.x() = rl1;
  upper = rl1;
  lower = upper/1000.;

  // Now compute Roche lobe in regular steps of angle looking
  // from centre of Roche lobe

  for(int i=0;i<n;i++){

    theta = 2*M_PI*(double)i/(n-1);
    double dx =  cos(theta);
    double dy =  sin(theta);
    
    
    // Locate critical surface using rtsafe.
    // Based on assuming that rl1 is maximum distance
    // from centre of mass and that at no point is the
    // surface closer than 1/4 of this.

    try{
      rad   = Subs::rtsafe(Lfunc1(dx,dy,q,pot), lower, upper, FRAC);
    }
    catch(...){
      throw Roche_Error("Encountered problem with rtsafe in Roche::flobe1");
    }
    
    x[i] = rad*dx;
    y[i] = rad*dy;
  }
}

