#include <cstdlib>
#include <cmath>
#include <iostream>
#include "trm/subs.h"
#include "trm/constants.h"
#include "trm/vec3.h"
#include "trm/roche.h"

/**
 * rpot_val computes the value of the Roche potential for a specific value of phi & lambda.
 * phi refers to the orbital phase, lambda to a multiplier that specified the position
 * of a point from an origin plus the multiplier time lambda.
 * \param q mass ratio  = M2/M1
 * \param star which star is relevant (for asynchronism)
 * \param spin ratio spin to orbital frequency
 * \param earth earth vector
 * \param p position of origin (units of separation)
 * \param lam multiplier
 * \return Roche potential at point.
 */

double Roche::rpot_val(double q, STAR star, double spin, const Subs::Vec3& earth,
                       const Subs::Vec3& p, double lam){

    if(q <= 0.)
        throw Roche_Error(Subs::str("Roche::rpot_val(double, const Subs::Vec3&") +
                          ", const Subs::Vec3&, double): q = " + Subs::str(q) + " <= 0.");

    Subs::Vec3 r = p + lam*earth;
    return star == PRIMARY ? rpot1(q, spin, r) : rpot2(q, spin, r);
}
