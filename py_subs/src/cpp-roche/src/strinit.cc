#include <cmath>
#include "trm/subs.h"
#include "trm/vec3.h"
#include "trm/roche.h"


/** strinit sets a particle just inside the L1 point with the 
 * correct velocity as given in Lubow and Shu.
 *
 * \param q mass ratio = M2/M1
 * \param r start position returned
 * \param v start velocity returned
 */

void Roche::strinit(double q, Subs::Vec3& r, Subs::Vec3& v){

  const double SMALL = 1.e-5;
  double rl1, mu, a, lambda1, m1;

  rl1 = Roche::xl1(q);
  mu  = q/(1.0+q);
  a   = (1.0-mu)/pow(rl1,3)+mu/pow((1.-rl1),3);
  lambda1 = sqrt(((a-2.0) + sqrt(a*(9.0*a-8.0)))/2.0);
  m1  = (lambda1*lambda1-2.0*a-1.0)/2.0/lambda1;

  r.x()  = rl1-SMALL;
  r.y()  = -m1*SMALL;
  r.z()  = 0.0;
  v.x()  = -lambda1*SMALL;
  v.y()  = -lambda1*m1*SMALL;
  v.z()  = 0.0;
}


