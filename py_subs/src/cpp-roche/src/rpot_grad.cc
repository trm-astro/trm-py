#include <cstdlib>
#include <cmath>
#include <iostream>
#include "trm/subs.h"
#include "trm/constants.h"
#include "trm/vec3.h"
#include "trm/roche.h"

/**
 * rpot_grad computes the gradient in phi, lambda space of the Roche potential.
 * phi refers to the orbital phase, lambda to a multiplier that specified the position
 * of a point from an origin plus the multiplier time lambda.
 * \param q mass ratio  = M2/M1
 * \param star which star is relevant (to allow for asynchronism)
 * \param spin ratio of spin to orbital frequency
 * \param cosi cosine of orbital inclination
 * \param sini sine of orbital inclination
 * \param iangle  the orbital inclination, degrees. 90 = edge on.
 * \param p position of origin (units of separation)
 * \param phi phase
 * \param lam multiplier
 * \param dphi first derivative of Roche potential wrt phi
 * \param dlam first derivative of Roche potential wrt lamda
 */

void Roche::rpot_grad(double q, STAR star, double spin, const Subs::Vec3& earth, const Subs::Vec3& p, double lam, double& dphi, double& dlam){

  if(q <= 0.) 
    throw Roche_Error("Roche::rpot_grad(double, STAR, double, const Subs::Vec3&, const Subs::Vec3&, double, double& double&, double&): q = " + Subs::str(q) + " <= 0.");

  Subs::Vec3 r = p + lam*earth;
  Subs::Vec3 d = star == PRIMARY ? Roche::drpot1(q, spin, r) : Roche::drpot2(q, spin, r);

  // derivative wrt phi
  Subs::Vec3 ed(earth.y(), -earth.x(), 0.);
  dphi  = Constants::TWOPI*lam*Subs::dot(d, ed);

  // derivative wrt lambda
  dlam  = Subs::dot(d, earth);

}
