#include <cmath>
#include <cstdlib>
#include <iostream>
#include "trm/subs.h"
#include "trm/constants.h"
#include "trm/vec3.h"
#include "trm/roche.h"

/**
 * fblink works out whether or not a given point is eclipsed by a Roche-distorted star by searching along the line
 * of sight to the point to see if the Roche potential ever drops below the value at the stellar surface.
 *
 * \param q      mass ratio = M2/M1
 * \param star   star concerned
 * \param spin   ratio of spin to orbital frequency
 * \param star   which star is doing the eclipsing, primary or secondary
 * \param ffac   the filling factor of the star
 * \param acc    accuracy of location of minimum potential, units of separation. The accuracy in height relative to the
 * Roche potential is acc*acc/(2*R) where R is the radius of curvature of the Roche potential surface, so don't be too
 * picky. 1.e-4 would be more than good enough in most cases.
 * \param earth  vector pointing towards earth
 * \param p      point of interest
 * \return true if minimum potential is below the potential at stellar surface
 */

bool Roche::fblink(double q, STAR star, double spin, double ffac, double acc, const Subs::Vec3& earth, const Subs::Vec3& p){

    // Compute radius of reference sphere and corresponding Roche potential.
    double rref, pref;
    ref_sphere(q, star, spin, ffac, rref, pref);

    Subs::Vec3 cofm;
    if(star == PRIMARY)
	cofm.set(0.,0.,0.);
    else
	cofm.set(1.,0.,0.);

    // First compute the multipliers cutting the reference sphere (if any)
    double lam1, lam2;
    if(!sphere_eclipse(earth, p, cofm, rref, lam1, lam2)) return false;
    if(lam1 == 0.) return true;

    // Create function objects for 1D minimisation in lambda direction
    Rlpot func(q, star, spin, earth, p);

    // Now try to bracket a minimum. We just crudely compute function at regularly spaced intervals filling in the
    // gaps until the step size between the points drops below the threshold. Take every opportunity to jump out early
    // either if the potential is below the threshold or if we have bracketed a minimum.
    int nstep   = 1;
    double step = (lam2-lam1);
    double f1 = 0., f2 = 0., flam = 1., lam = lam1;

    while(step > acc){

	lam = lam1 + step/2.;
	for(int n=0; n<nstep; n++){

	    flam = func(lam);
	    if(flam <= pref) return true;

	    // Calculate these as late as possible because they may often not be needed
	    if(nstep == 1){
		f1 = func(lam1);
		f2 = func(lam2);
	    }
	    if(flam < f1 && flam < f2) break;
	    lam += step;
	}
	if(flam < f1 && flam < f2) break;
	step  /= 2.;
	nstep *= 2;
    }

    if(flam < f1 && flam < f2) {

	// OK, minimum bracketted, so finally pin it down accurately
	// Possible that multiple minima could cause problems but I have
	// never seen this in practice.
	Dlrpot dfunc(q, star, spin, earth, p);
	double xmin;
	try {
	    flam = dbrent(lam1, lam, lam2, func, dfunc, acc, true, pref, xmin);
	}
	catch(const Subs::Subs_Error& err){
	    std::cerr << "Line minimisation error inside fblink" << std::endl;
	    std::cerr << err << std::endl;
	    std::cerr << "q = " << q << ", earth = " << earth << ", point = " << p << ", star = " << int(star) 
		      << ", ffac = " << ffac << ", acc = " << acc << std::endl;
	    throw Subs::Subs_Error("fblink died");
	}
	return (flam < pref);

    }else{

	// Not bracketted even after a detailed search, and we have not jumped 
	// out either, so assume no eclipse
	return false;

    }
}

