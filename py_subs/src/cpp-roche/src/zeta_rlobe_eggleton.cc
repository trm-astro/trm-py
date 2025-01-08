#include <math.h>
#include "trm/roche.h"

/**
 * zeta_rlobe_eggleton returns d log(rl) / d log (m2)
 * where rl is Peter Eggleton's formula for the volume-averaged 
 * Roche lobe radius divided by the orbital separation. This assumes
 * that m1+m2 = constant.
 *
 * \param q mass ratio = M2/M1
 * \return Returns d log(rl) / d log (m2) where rl = Roche lobe radius of
 * secondary star divided by separation according to Eggleton's formula.
 */

double Roche::zeta_rlobe_eggleton(double q){
  double q1 = pow(q,1./3.);
  double loneq = log(1+q1);
  return ((1.+q)/3.*(2.*loneq-q1/(1.+q1))/(0.6*q1*q1+loneq));
}


/**
 * dzetadq_rlobe_eggleton returns d zeta / d q where zeta is the result 
 * of zeta_rlobe_eggleton(double q). This has been tested successfully
 * against finite difference value.
 *
 * \param q mass ratio = M2/M1
 * \return Returns d zeta d q
 */

double Roche::dzetadq_rlobe_eggleton(double q){
  double q1    = pow(q,1./3.);
  double q2    = q1*q1;
  double opq1  = 1.+q1;
  double loneq = log(opq1);
  double denom = 0.6*q2+loneq;
  double numer = 2.*loneq-q1/opq1;
  return (numer/denom/3. + (1.+q)/3.*((1.+2.*q1)/3./pow(q1*opq1,2) - 
				      numer*(0.4/q1+1./(3.*q2*(1.+q1)))/denom)/denom);
}





