#include <cstdlib>
#include <cmath>
#include <iostream>
#include "trm/subs.h"
#include "trm/constants.h"
#include "trm/vec3.h"
#include "trm/roche.h"

/**
 * rpot_val_grad computes the value & gradient in phi, lambda space of the Roche potential.
 * phi orbital phase, lambda a multiplier that specified the position
 * of a point from an origin plus the multiplier times lambda.
 * \param q mass ratio  = M2/M1
 * \param star which star
 * \param spin ratio spi/orbital frequency
 * \param earth vector towards earth (defines by phase and inclination)
 * \param p position of origin (units of separation)
 * \param lam multiplier
 * \param rpot the Roche potential
 * \param dphi first derivative of Roche potential wrt phi
 * \param dlam first derivative of Roche potential wrt lambda
 */

void Roche::rpot_val_grad(double q, STAR star, double spin, const Subs::Vec3& earth,
                          const Subs::Vec3& p, double lam, double& rpot, double& dphi, double& dlam){

    if(q <= 0.)
        throw Roche_Error(Subs::str("Roche::rpot_grad(double, STAR, double, const Subs::Vec3&,") +
                          " const Subs::Vec3&, double, double& double&, double&): q = " +
                          Subs::str(q) + " <= 0.");

    Subs::Vec3 r = p + lam*earth;
    Subs::Vec3 d = (star == PRIMARY) ? Roche::drpot1(q, spin, r) : Roche::drpot2(q, spin, r);
    rpot = (star == PRIMARY) ? Roche::rpot1(q, spin, r) : Roche::rpot2(q, spin, r);

    // Derivative of earth wrt phi
    Subs::Vec3 ed(earth.y(), -earth.x(), 0.);

    // derivative wrt phi
    dphi  = Constants::TWOPI*lam*Subs::dot(d, ed);

    // derivative wrt lambda
    dlam  = Subs::dot(d, earth);

}
