#include <cstdlib>
#include <cmath>
#include <iostream>
#include "trm/subs.h"
#include "trm/vec3.h"
#include "trm/roche.h"

/**
 * roche_shadow returns arrays x and y for plotting the shadow of
 * a semi-detached star's Roche lobe in the equatorial plane.
 * \param q mass ratio = M2/M1
 * \param iangle inclination angle
 * \param phi orbital phase
 * \param dist maximum distance to search for shadow measured
 * from center of mass of red star, units of separation.
 * \param acc accuracy on location of shadow distance from centre of red star, units of separation. The shadow
 * edge is located by binary chop.
 * from center of mass of red star.
 * \param x     array of x values
 * \param y     array of y values
 * \param shade true/false if genuine shae or not. The array goes all the way round and when not in shade it
 * will be glued to the red star. This array allows you to see if this is the case or not.
 * \param n number of x and y values. The points will be equally spaced in angle around centre of
 * mass of the red star.
 */
void Roche::roche_shadow(double q, double iangle, double phi, double dist, double acc, float x[], float y[], bool shade[], int n){

    if(q <= 0.)
	throw Roche_Error("Roche::roche_shadow(double, double, double, double, double, float*, float*, int): q = " + Subs::str(q) + " <= 0.");

    if(n < 1) 
	throw Roche_Error("Roche::roche_shadow(double, double, double, double, double, float*, float*, int): n = " + Subs::str(n) + " < 1.");

    // Compute L1 point and critical potential there.
    double rl1     = xl1(q);
    Subs::Vec3 earth = set_earth(iangle, phi);
    Subs::Vec3 cofm(1.,0.,0.), p, dirn, p1, p2;
    p.y()   = p.z() = 0.;
    p.x()   = rl1;
    double cpot  = rpot(q,p);

    // Limits for Roche lobe computation.
    double upper = 1.-rl1;
    double lower = upper/4.;

    // Now compute Roche lobe in regular steps of angle looking
    // from centre of Roche lobe, then search fro the shadow in between
    // the Roche lobe and the maximum distance 

    double r1, r2;
    for(int i=0; i<n; i++){

	// L1 point is a special case because derivative becomes zero there.
	// lambda is set so that after i=0, there is a decent starting 
	// multiplier

	double theta = 2*M_PI*double(i)/(n-1);
	double dx    = -cos(theta);
	double dy    =  sin(theta);
	if(i == 0 || i == n-1){
	    r1 = 1-rl1 + acc;
	}else{

	    // Locate critical surface using rtsafe.
	    // Based on assuming that rl1 is maximum distance
	    // from centre of mass and that at no point is the
	    // surface closer than 1/4 of this.
	    r1   = Subs::rtsafe(Lfunc2(dx,dy,q,cpot), lower, upper, acc) + acc;
	}
	r2 = dist;

	// First check status of end points
	dirn.set(dx,dy,0.);
	p1 = cofm + r1*dirn;
	p2 = cofm + r2*dirn;

	if(!fblink(q, SECONDARY, 1., 1., acc, earth, p1)){
	    x[i]     = p1.x();
	    y[i]     = p1.y();
	    shade[i] = false;
	}else if(fblink(q, SECONDARY, 1., 1., acc, earth, p2)){
	    x[i]     = p2.x();
	    y[i]     = p2.y();
	    shade[i] = true;
	}else{
	    while(r2-r1 > acc){
		p = (p1+p2)/2;
		if(fblink(q, SECONDARY, 1., 1., acc, earth, p)){
		    p1 = p;
		    r1 = (r1+r2)/2;
		}else{
		    p2 = p;
		    r2 = (r1+r2)/2;
		}
	    }
	    p = (p1+p2)/2;
	    x[i]     = p.x();
	    y[i]     = p.y();
	    shade[i] = true;
	}
    }
}


