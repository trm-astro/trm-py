#include <cstdlib>
#include <cmath>
#include <iostream>
#include <string>
#include "trm/subs.h"
#include "trm/constants.h"
#include "trm/roche.h"

/**
 * sphere_eclipse tells you whether or not a given sphere will eclipse a 
 * given point or not. If the answer is yes, it will also return with four 
 * parameters to define the phase range and the multiplier range delimiting the
 * region within which the spheres surface is crossed.
 * These can then be used as the starting point for later computation.
 * (The line of sight is described as the point in question plus a scalar multiplier
 * times a unit vector pointing towards Earth -- this is the "multiplier" referred to above
 * and below). The multiplier must be positive: in other words the routine does not
 * project backwards. If the point in inside the sphere, phi1 will be set = 0, phi2 = 1,
 * lam1 = 0, and lam2 = the largest value of the multiplier lambda
 *
 * \param cosi cosine of orbital inclination
 * \param sini sine of orbital inclination
 * \param r       the position vector of the point in question (units of binary separation)
 * \param c       the centre of the sphere enclosing the star (units of binary separation)
 * \param rsphere the radius defining the sphere enclosing the star (units of binary separation)
 * \param phi1    if eclipse by the sphere, this is the start of the phase range (0 -1)
 * \param phi2    if eclipse by the sphere, this is the end of the phase range (least amount > phi1)
 * \param lam1    if eclipse by the sphere, this is the start of the multiplier range (>=0)
 * \param lam2    if eclipse by the sphere, this is the end of the multiplier range (>lam1)
 * \return false = not eclipsed; true = eclipsed.
 */

bool Roche::sphere_eclipse(double cosi, double sini, const Subs::Vec3& r, const Subs::Vec3& c, double rsphere, 
			   double& phi1, double& phi2, double& lam1, double& lam2){

    // Compute vector from centre of sphere to point
    const Subs::Vec3 d = r - c;

    // Test whether any eclipse will take place
    const double PDIST = sqrt(Subs::sqr(d.x())+Subs::sqr(d.y()));

    // This is half the minimum value of the linear coefficient in the quadratic that determines
    // whether an intersection occurs. We use the minimum to increase chance of positive roots
    // for the multiplier.
    const double BQUAD = d.z()*cosi - PDIST*sini;
    if(BQUAD >= 0.) return false;

    // If this is < 0, the point is inside the sphere
    const double CQUAD = d.sqr() - Subs::sqr(rsphere);

    double fac = Subs::sqr(BQUAD) - CQUAD;
    if(fac <= 0.) return false;

    // OK, it does cross
    fac = sqrt(fac);

    // Finally work out multiplier limits. Specialisation of NR
    // method to avoid subtraction of similar sized quantities.
    // lam2 is the usual '(-b+sqrt(b^2-4ac))/2a' root 
    lam2 = -BQUAD + fac;
    lam1 = std::max(0., CQUAD/lam2);
	
    // Now compute phase range within which eclipse occurs
    if(CQUAD < 0.){
	phi1 = 0.;
	phi2 = 1.;
    }else{
	double delta = acos((cosi*d.z()+sqrt(CQUAD))/(sini*PDIST));
	double phi   = atan2(d.y(), -d.x());
	phi1 = (phi - delta)/Constants::TWOPI;
	phi1 -= floor(phi1);
	phi2 = phi1 + 2.*delta/Constants::TWOPI;
    }
    return true;
}

/**
 * This version of sphere_eclipse tells you whether or not a given sphere will eclipse 
 * a given point at a particular phase or not. If the answer is yes,
 * it will also return with the multiplier values giving the cut points. These can then
 * be used as starting points for Roche lobe computations. These can then be used as the 
 * starting point for later computation. Points inside the sphere are regarded as being
 * eclipsed with the lower mulitplier set = 0
 *
 * \param earth   vector towards Earth
 * \param r       the position vector of the point in question (units of binary separation)
 * \param c       the centre of the sphere (units of binary separation)
 * \param rsphere the radius defining the sphere enclosing the star (units of binary separation)
 * \param lam1    if eclipse by the sphere, this is the start of the multiplier range
 * \param lam2    if eclipse by the sphere, this is the end of the multiplier range
 * \return false = not eclipsed; true = eclipsed.
 */

bool Roche::sphere_eclipse(const Subs::Vec3& earth, const Subs::Vec3& r, const Subs::Vec3& c, double rsphere, double& lam1, double& lam2){

    // Compute vector from centre of sphere to point
    const Subs::Vec3 d = r - c;

    // Work out whether the line of sight cuts sphere
    const double BQUAD = Subs::dot(earth, d);
    if(BQUAD >= 0.) return false;

    const double CQUAD = d.sqr() - Subs::sqr(rsphere);
    double fac = Subs::sqr(BQUAD) - CQUAD;
    if(fac <= 0.) return false;

    // OK it does cross
    fac = sqrt(fac);

    // Finally work out multiplier limits. Specialisation of NR
    // method to avoid subtraction of similar sized quantities.
    // lam1 must be >= 0
    lam2 = -BQUAD + fac;
    lam1 = std::max(0., CQUAD/lam2);
    return true;
}





