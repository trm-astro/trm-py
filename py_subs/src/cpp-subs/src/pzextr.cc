#include <cstdlib>
#include "trm/subs.h"

extern double **d, *x; // Defined in bsstep

/**
 * pzextr is a polynomial extrapolation routine used by bsstep. It
 * extrapolates a sequence of estimates with progressively smaller
 * values of x = xest.
 * \param iest the call we are now on
 * \param xest the present value of x
 * \param yest the present values of y
 * \param yz extrapolated values of y
 * \param dy estimates errors on the yz
 * \param nv number of y values
 */

void Subs::pzextr(int iest, double xest, double yest[], double yz[], double dy[], int nv){
  int k1, j;
  double q, f2, f1, delta, c[nv];
  
  x[iest] = xest;
  for(j=0;j<nv;j++) dy[j] = yz[j] = yest[j];
  if(iest == 0){
    for(j=0;j<nv;j++) d[j][0] = yest[j];
  }else{
    for(j=0;j<nv;j++) c[j]=yest[j];
    for(k1=0;k1<iest;k1++){
      delta = 1.0/(x[iest-k1-1]-xest);
      f1 = xest*delta;
      f2 = x[iest-k1-1]*delta;
      for(j=0;j<nv;j++){
	q        = d[j][k1];
	d[j][k1] = dy[j];
	delta    = c[j]-q;
	dy[j]    = f1*delta;
	c[j]     = f2*delta;
	yz[j]   += dy[j];
      }
    }
    for(j=0;j<nv;j++) d[j][iest]=dy[j];
  }
  return;
}




