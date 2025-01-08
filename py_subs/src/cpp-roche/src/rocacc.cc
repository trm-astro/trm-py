#include <cmath>
#include "trm/subs.h"
#include "trm/vec3.h"
#include "trm/roche.h"


/**
 * rocacc calculates and returns the acceleration (in the rotating frame)
 * in a Roche potential of a particle of given position and velocity.
 *
 * \param q mass ratio = M2/M1
 * \param r position, scaled in units of separation.
 * \param v velocity, scaled in units of separation
 */

Subs::Vec3 Roche::rocacc(double q, const Subs::Vec3& r, const Subs::Vec3& v){

  double f1, f2, yzsq, r1sq, r2sq, fm1, fm2, fm3;

  f1 = 1.0/(1.0+q);
  f2 = f1*q;

  yzsq = Subs::sqr(r.y()) + Subs::sqr(r.z());
  r1sq = Subs::sqr(r.x()) + yzsq;
  r2sq = Subs::sqr(r.x()-1.0) + yzsq;
  fm1  = f1/(r1sq*sqrt(r1sq));
  fm2  = f2/(r2sq*sqrt(r2sq));
  fm3  = fm1+fm2;

  Subs::Vec3 tmp;
  tmp.x() = -fm3*r.x() + fm2 + 2.0*v.y() + r.x() - f2;
  tmp.y() = -fm3*r.y()       - 2.0*v.x() + r.y();
  tmp.z() = -fm3*r.z();
  return tmp;
}
