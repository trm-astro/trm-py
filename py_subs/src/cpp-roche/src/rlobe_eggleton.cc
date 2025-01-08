#include <cmath>
#include <string>
#include "trm/subs.h"
#include "trm/roche.h"

/**
 * rlobe_eggleton is a function which returns Peter Eggleton's formula
 * for the volume-averaged Roche lobe radius
 * \param q mass ratio = m2/m1.
 */

double Roche::rlobe_eggleton(double q){
  if(q <= 0.)
    throw Roche_Error("Roche::rlobe_eggleton(double): q = " + Subs::str(q) + " <= 0.");
  double q3 = pow(q,1./3.);
  return (0.49*q3*q3/(0.6*q3*q3+log(1.+q3)));
}

