/*
 
Routine from NR for speeding up integration of ODEs when dealing with second-order
conservative equations. Replaces mmid used by more general method. 

*/
#include <cstdlib>
#include "trm/subs.h"

void Subs::stoerm(double y[], double d2y[], int nvar, double xs, double htot,
                  int nstep, double yout[], const Bsfunc& derivs){

    int i, n, neqns, nn;
    double h, h2, halfh, x, *ytemp = new double[nvar];

    h = htot / nstep;
    halfh = 0.5*h;
    neqns = nvar/2;
    for(i=0;i<neqns;i++){
        n        = neqns+i;
        ytemp[n] = h*(y[n]+halfh*d2y[i]);
        ytemp[i] = y[i] + ytemp[n];
    }
    x = xs + h;
    derivs(x,ytemp,yout);
    h2 = h*h;
    for(nn=1;nn<nstep;nn++){
        for(i=0;i<neqns;i++){
            n         = neqns+i;
            ytemp[n] += h2*yout[i];
            ytemp[i] += ytemp[n];
        }
        x += h;
        derivs(x,ytemp,yout);
    }
    for(i=0;i<neqns;i++){
        n       = neqns + i;
        yout[n] = ytemp[n]/h + halfh*yout[i];
        yout[i] = ytemp[i];
    }
    delete[] ytemp;
}
