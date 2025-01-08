#include <cstdlib>
#include <cmath>
#include "trm/constants.h"
#include "trm/subs.h"
#include "trm/roche.h"

/**
 * eprob  computes probability that accretor is eclipsed by donor in a Roche lobe
 * filling binary assuming randomly oriented orbits. Accretor assumed to have zero size.
 * \param q mass ratio = M2/M1
 * \return Returns the probability that the centre of accretor is eclipsed.
 */

double Roche::eprob(double q){
  double i1=0., i2=90., a, i, ir;
  coord r, e;

  e.y  = r.x = r.y = r.z = 0.;
  while(i2-i1 > 0.001){
    ir = Constants::PI*(i = (i1+i2)/2.)/180.;
    e.x = sin(ir);
    e.z = cos(ir);
    if(blink(q,r,e,0.01)){
      i2 = i;
    }else{
      i1 = i;
    }
  }
  return e.z;
}



  
