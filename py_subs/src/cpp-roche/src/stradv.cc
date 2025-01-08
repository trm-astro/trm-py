#include <cstdlib>
#include <cmath>
#include <iostream>
#include "trm/subs.h"
#include "trm/vec3.h"
#include "trm/roche.h"


/**
 * stradv advances a particle of given position and velocity until
 * it reaches a specified radius. It then returns with updated position and
 * velocity. It is up to the user not to request a value that cannot be reached.
 *
 * \param q    mass ratio = M2/M1
 * \param r    Initial and final position
 * \param v    Initial and final velocity
 * \param rad  Radius to aim for
 * \param acc  Accuracy with which to place output point at rad.
 * \param smax Largest time step allowed. It is possible that the
 * routine could take such a large step that it misses
 * the point when the stream is inside the requested
 * radius. This allows one to control this. Typical
 * value = 1.e-3.
 *
 * Returns time step taken
 */

double Roche::stradv(double q, Subs::Vec3 &r, Subs::Vec3 &v, double rad,
                   double acc, double smax){

    const  double EPS    = 1.e-8;
    const  double TMAX   = 10.;
    static double tnext  = 1.e-2;
    double rnow, ttry, tdid, time, lo, hi, rlo, rhi, rinit;
    double to;
    Subs::Vec3 ro, vo;

    // Store initial radius
    rinit = rnow = r.length();

    // Step until radius crossed
    time  = 0.0;
    while((rinit > rad && rnow > rad) || (rinit < rad && rnow < rad)){
        ro = r;
        vo = v;
        ttry = std::min(tnext, smax);
        gsint(q, r, v, ttry, tdid, tnext, time, EPS);
        rnow = r.length();
        if(time > TMAX) throw Roche_Error("Roche::stradv: taken too long without crossing radius = " + Subs::str(rad) + " in Roche::stradv.\n");
    }

    /*
      Now refine by reinitialising and binary chopping until
      close enough to requested radius.
    */

    lo  = 0.0;
    hi  = tdid;
    rlo = ro.length();
    rhi = rnow;
    to  = time;
    while(fabs(rhi-rlo) > acc){
        ttry   = (lo+hi)/2.;
        r = ro;
        v = vo;
        time = to;
        gsint(q, r, v, ttry, tdid, tnext, time, EPS);
        rnow = r.length();
        if((rhi > rad && rnow > rad) || (rhi < rad && rnow < rad)){
            rhi = rnow;
            hi  = ttry;
        }else{
            rlo = rnow;
            lo  = ttry;
        }
    }
    return time;
}



