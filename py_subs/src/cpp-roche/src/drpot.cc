#include <cstdlib>
#include <cmath>
#include <iostream>
#include "trm/subs.h"
#include "trm/vec3.h"
#include "trm/roche.h"

/**
 * drpot computes partial derivatives of Roche potential with respect
 * to position at p for mass ratio q.
 * \param q mass ratio  = M2/M1
 * \param p position (units of separation)
 */

Subs::Vec3 Roche::drpot(double q, const Subs::Vec3& p){

  double mu, comp, r1, r2, r1sq, r2sq, mu1;

  if(q <= 0.) throw Roche_Error("Roche::drpot(double, const Subs::Vec3&): q = " + Subs::str(q) + " <= 0.");

  r1    = sqrt(r1sq = p.sqr());
  r2    = sqrt(r2sq = r1sq + 1. - 2.*p.x());
  mu    = q/(1+q);
  mu1   = mu/r2/r2sq;
  comp  = (1.-mu)/r1/r1sq;
  Subs::Vec3 tmp;
  tmp.x() = comp*p.x() + mu1*(p.x()-1.) -  p.x() + mu;
  tmp.y() = comp*p.y() + mu1*p.y()      -  p.y();
  tmp.z() = comp*p.z() + mu1*p.z();
  return tmp;
}

/**
 * drpot1 computes partial derivatives of asynchronous Roche potential with respect
 * to position at p for mass ratio q, star 1.
 * \param q mass ratio  = M2/M1
 * \param spin ratio of spin to orbital frequency.
 * \param p position (units of separation)
 */

Subs::Vec3 Roche::drpot1(double q, double spin, const Subs::Vec3& p){

    double mu, comp, r1, r2, r1sq, r2sq, mu1, ssq;

  if(q <= 0.) throw Roche_Error("Roche::drpot1(double, double, const Subs::Vec3&): q = " +
                                Subs::str(q) + " <= 0.");

  r1    = sqrt(r1sq = p.sqr());
  r2    = sqrt(r2sq = r1sq + 1. - 2.*p.x());
  mu    = q/(1+q);
  mu1   = mu/r2/r2sq;
  comp  = (1.-mu)/r1/r1sq;
  ssq   = Subs::sqr(spin);
  Subs::Vec3 tmp;
  tmp.x() = comp*p.x() + mu1*(p.x()-1.) -  ssq*p.x() + mu;
  tmp.y() = comp*p.y() + mu1*p.y()      -  ssq*p.y();
  tmp.z() = comp*p.z() + mu1*p.z();
  return tmp;
}

/**
 * drpot2 computes partial derivatives of asynchronous Roche potential with respect
 * to position at p for mass ratio q, star 2.
 * \param q mass ratio  = M2/M1
 * \param spin ratio of spin to orbital frequency.
 * \param p position (units of separation)
 */

Subs::Vec3 Roche::drpot2(double q, double spin, const Subs::Vec3& p){

    double mu, comp, r1, r2, r1sq, r2sq, mu1, ssq;

  if(q <= 0.) throw Roche_Error("Roche::drpot2(double, double, const Subs::Vec3&): q = " + Subs::str(q) + " <= 0.");

  r1    = sqrt(r1sq = p.sqr());
  r2    = sqrt(r2sq = r1sq + 1. - 2.*p.x());
  mu    = q/(1+q);
  mu1   = mu/r2/r2sq;
  comp  = (1.-mu)/r1/r1sq;
  ssq   = Subs::sqr(spin);
  Subs::Vec3 tmp;
  tmp.x() = comp*p.x() + mu1*(p.x()-1.) -  ssq*(p.x()-1.) + mu - 1.;
  tmp.y() = comp*p.y() + mu1*p.y()      -  ssq*p.y();
  tmp.z() = comp*p.z() + mu1*p.z();
  return tmp;
}
