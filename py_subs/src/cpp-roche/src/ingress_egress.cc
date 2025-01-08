#include <cstdlib>
#include <cmath>
#include <iostream>
#include <string>
#include "trm/subs.h"
#include "trm/constants.h"
#include "trm/vec3.h"
#include "trm/roche.h"

/**
 * ingress_egress tests for whether a given point is eclipsed by a Roche-distorted star. If 
 * it is, it computes the ingress and egress phases using a binary chop. The accuracy on the 
 * phase should be set to be below the expected uncertainties of the phases of your data.
 *
 * \param q       the mass ratio = M2/M1.
 * \param star    which star, primary or secondary, is doing the eclipsing
 * \param spin    ratio of spin to orbit of eclipsing star
 * \param ffac    linear filling factor
 * \param iangle  inclination angle
 * \param delta   the accuracy in phase wanted.
 * \param r       position vector of point of interest.
 * \param ingress ingress phase (if eclipsed)
 * \param egress  egress phase 
 * \return false = not eclipsed; true = eclipsed.
 */

bool Roche::ingress_egress(double q, STAR star, double spin, double ffac, double iangle, double delta, const Subs::Vec3& r, double& ingress, double& egress){

    // Compute radius of reference sphere and corresponding Roche potential.
    double rref, pref;
    ref_sphere(q, star, spin, ffac, rref, pref);

    double ri = Subs::deg2rad(iangle);
    double cosi = cos(ri), sini = sin(ri);

    Subs::Vec3 cofm;
    if(star == PRIMARY)
	cofm.set(0.,0.,0.);
    else
	cofm.set(1.,0.,0.);

    // Test for eclipse by reference sphere at any point at all
    double phi1, phi2, lam1, lam2, phi, lam;

    if(sphere_eclipse(cosi, sini, r, cofm, rref, phi1, phi2, lam1, lam2)){
	  
	// A certain accuracy of position corresponds to an accuracy in phase
	const double ACC = 2.*sqrt(2.*Constants::TWOPI*(lam2-lam1)*delta);
	
	// Now using the limits from the reference sphere, check more carefully
	if(pot_min(q, cosi, sini, star, spin, r, phi1, phi2, lam1, lam2, rref, pref, ACC, phi, lam)){
	    
	    // We now know that there is an eclipse, and one phase (=phi) when it occurs. 
	    // Refine boundaries with binary chop
	    
	    // Start with the ingress
	    double pin = phi, pout = phi1, pmid;
	    while(fabs(pin-pout) > delta){
		pmid = (pin+pout)/2.;
		if(fblink(q, star, spin, ffac, ACC, set_earth(cosi, sini, pmid), r)){
		    pin  = pmid;
		}else{
		    pout = pmid;
		}
	    }
	    ingress = (pin+pout)/2.;
	    ingress = ingress - floor(ingress);
	    
	    // Now the egress
	    pin  = phi;
	    pout = phi2;
	    while(fabs(pin-pout) > delta){
		pmid = (pin+pout)/2.;
		if(fblink(q, star, spin, ffac, ACC, set_earth(cosi, sini, pmid), r)){
		    pin  = pmid;
		}else{
		    pout = pmid;
		}
	    }
	    egress  = (pin+pout)/2.;
	    egress  = egress  - floor(egress);
	    if(egress < ingress) egress++;
	    
	    return true;
	}else{
	    // Eclipse by reference sphere, but not by star
	    return false;
	}
    }else{
	// No eclipse by reference sphere
	return false;
    }
}
