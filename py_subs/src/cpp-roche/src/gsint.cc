#include <stdlib.h>
#include <math.h>
#include "trm/roche.h"
#include "trm/subs.h"

void rocder(double t, double y[], double dydt[]);

// Mass ratio for communication with rocder 
static double qp; 

/**
 * gsint carries out one step of orbit integration in a Roche potential.
 * \param q mass ratio  = M2/M1
 * \param r initial and final position (separation units)
 * \param v initial and final velocity.
 * \param ttry time step to try
 * \param tdid actual time step carried out.
 * \param tnext next time step to try.
 * \param time updated time.
 * \param eps accuracy with which to carry out step.
 */
void Roche::gsint(double q, Subs::Vec3 &r, Subs::Vec3 &v, double ttry, 
		  double &tdid, double &tnext, double &time, 
		  double eps){

  double y[6], dydt[6];
  static double yscal[6] = {1.0,1.0,1.0,1.0,1.0,1.0};

  // For communication with rocder 

  qp = q; 

  // Set initial values

  y[0] = r.x();
  y[1] = r.y();
  y[2] = r.z();
  y[3] = v.x();
  y[4] = v.y();
  y[5] = v.z();

  // Compute derivatives

  rocder(time,y,dydt);

  // Carry out a step with general B-S routine 

  Subs::bsstep(y,dydt,6,time,ttry,eps,yscal,tdid,tnext,rocder);

  // Update values for output 

  r.x() = y[0];
  r.y() = y[1];
  r.z() = y[2];
  v.x() = y[3];
  v.y() = y[4];
  v.z() = y[5];
  
  return;
}

/* 
   Function returning derivatives in correct form for bsstep etc. 
   Denoted by (*derivs) in NR. Note that the time t does not enter
   into this calculation at all but is kept for consistency.
   
   */

void rocder(double t, double y[], double dydt[]){
  
  Subs::Vec3 r, v, a;

  r.x() = y[0];
  r.y() = y[1];
  r.z() = y[2];
  dydt[0] = v.x() = y[3];
  dydt[1] = v.y() = y[4];
  dydt[2] = v.z() = y[5];

  a = Roche::rocacc(qp, r, v);

  dydt[3] = a.x();
  dydt[4] = a.y();
  dydt[5] = a.z();
  return;
}



