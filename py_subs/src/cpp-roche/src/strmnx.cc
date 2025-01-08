#include <cmath>
#include "trm/subs.h"
#include "trm/vec3.h"
#include "trm/roche.h"

/**
 * strmnx finds the next point at which stream is closest or furthest
 * from primary.
 * \param 1 mass ratio = M2/M1
 * \param r initial and final position
 * \param v initial and final velocity
 * \param acc accuracy in time to locate minimum/maximum.
 */
void Roche::strmnx(double q, Subs::Vec3 &r, Subs::Vec3 &v, double acc){

    const double EPS = 1.e-8;
    double dir1, dir, ttry, tdid, tnext, time, lo, hi;
    Subs::Vec3 ro, vo;

    // Store initial direction

    dir = dir1  = Subs::dot(r, v);

    time  = 0.0;  
    ttry  = 1.e-2;

    // Step until direction reverses

    while((dir > 0.0 && dir1 > 0.0) || (dir < 0.0 && dir1 < 0.0)){
	ro = r;
	vo = v;
	gsint(q, r, v, ttry, tdid, tnext, time, EPS);
	dir  = Subs::dot(r, v);
	ttry = tnext;
    }

    //   Now refine by reinitialising and binary chopping until
    //   close enough to requested radius.

    lo  = 0.0;
    hi  = tdid;
    while(fabs(hi-lo) > acc){
	ttry   = (lo+hi)/2.;
	r = ro;
	v = vo;
	gsint(q, r, v, ttry, tdid, tnext, time, EPS);
	dir  = Subs::dot(r, v);
	if((dir1 > 0.0 && dir < 0.0) || (dir1 < 0.0 && dir > 0.0)){
	    hi  = ttry;
	}else{
	    lo  = ttry;
	}
    }
}





