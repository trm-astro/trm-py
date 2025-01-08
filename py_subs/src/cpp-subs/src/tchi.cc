#include <cmath>
#include "trm/subs.h"

/**
 * tchi returns a number above which there is a 10**lfprob of 
 * getting chi**2 of ndof degrees of freedom by random chance alone.
 * \param lfprob log10(false alarm probability) i.e. -2 = 1%
 * \param ndof number of degrees of freedom
 */

double Subs::tchi(double lfprob, int ndof){

  // First, get upper limit to threshold
  
  double t1=0.,t2=ndof, t, sig = sqrt(2.*ndof);
  do{
    t2 += sig;
  }while(log10(gammq(ndof/2.,t2/2.)) > lfprob);
  
  // Second, binary chop

  do{
    t = (t1+t2)/2.;
    if(log10(gammq(ndof/2.,t/2.)) > lfprob){
      t1 = t;
    }else{
      t2 = t;
    }
  }while(fabs(t2-t1) > 0.00001*sig);
  return (t1+t2)/2.;
}






