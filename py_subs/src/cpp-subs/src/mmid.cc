/*

!!begin 
!!title  Modified mid-point routine
!!author T.R.Marsh
!!head1 Modified mid-point routine

Modified mid-point method routine from NR. This advances a vector of
variables y from x to x+htot in a series of nstep sub-steps.

!!head2 Function call

void mmid(double y[], double dydx[], int nvar, double xs, double htot, 
int nstep, double yout[], void (*derivs)(double, double[], double[]));

!!head2 Arguments

!!table
!!arg y[nvar]    !! y values
!!arg dydx[nvar] !! derivative values
!!arg nvar       !! number of y values and derivatives
!!arg xs         !! start x value
!!arg htot       !! total step to take
!!arg nstep      !! number of sub-steps
!!arg yout[nvar] !! output y values
!!arg derivs     !! derivative computation function

!!end

*/

#include <cstdlib>
#include "trm/subs.h"

void Subs::mmid(double y[], double dydx[], int nvar, double xs, double htot, 
		int nstep, double yout[], 
		void (*derivs)(double, double[], double[])){
  int n, i;
  double x,swap,h2,h, *ym = new double[nvar], *yn = new double[nvar];

  h = htot/nstep;
  for(i=0;i<nvar;i++){
      ym[i] = y[i];
      yn[i] = y[i]+h*dydx[i];
  }
  x=xs+h;
  (*derivs)(x, yn, yout);
  h2 = 2.0*h;
  for(n=1;n<nstep;n++){
      for(i=0;i<nvar;i++){
          swap  = ym[i]+h2*yout[i];
          ym[i] = yn[i];
          yn[i] = swap;
      }
      x += h;
      (*derivs)(x, yn, yout);
  }
  for(i=0;i<nvar;i++)
      yout[i]=0.5*(ym[i]+yn[i]+h*yout[i]);
  delete[] ym;
  delete[] yn;
  return;
}

void Subs::mmid(double y[], double dydx[], int nvar, double xs, double htot, 
		int nstep, double yout[], const Bsfunc& derivs){

  int n, i;
  double x,swap,h2,h, *ym = new double[nvar], *yn = new double[nvar];

  h = htot/nstep;
  for(i=0;i<nvar;i++){
    ym[i] = y[i];
    yn[i] = y[i]+h*dydx[i];
  }
  x=xs+h;
  derivs(x, yn, yout);
  h2 = 2.0*h;
  for(n=1;n<nstep;n++){
    for(i=0;i<nvar;i++){
      swap  = ym[i]+h2*yout[i];
      ym[i] = yn[i];
      yn[i] = swap;
    }
    x += h;
    derivs(x, yn, yout);
  }
  for(i=0;i<nvar;i++)
    yout[i]=0.5*(ym[i]+yn[i]+h*yout[i]);
  delete[] ym;
  delete[] yn;
  return;
}



