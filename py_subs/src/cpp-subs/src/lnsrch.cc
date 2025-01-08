//
// Given an n-dimensional point xold[n], the value of the function there,
// fold, and its gradient, g[n], a direction p[n], lnsrch finds a new point
// x[n] along the direction p[n] where the function func has decreased
// sufficiently.
//

#include <math.h>
#include <iostream.h>
#include "useful.h"

void lnsrch(int n, double xold[], double fold, double g[], double p[],
	    double x[], double &f, double stpmax, int &check, 
	    double (*func)(double [])){
  int i;
  double a, alam, alam2, alamin, b, disc, f2, fold2, rhs1, rhs2, slope;
  double sum, temp, test, tmplam;
  const double TOLX = 1.e-10;
  const double ALF  = 1.e-4;

  check = 0;
  for(sum=0.,i=0;i<n;i++) sum += p[i]*p[i];
  sum = sqrt(sum);
  if(sum > stpmax)
    for(i=0;i<n;i++) p[i] *= stpmax/sum;
  for(slope=0.,i=0;i<n;i++)
    slope += g[i]*p[i];
  
  test = 0.;
  for(i=0;i<n;i++){
    temp = fabs(p[i])/max(fabs(xold[i]),1.);
    test = max(test, temp);
  }
  alamin = TOLX/test;
  alam   = 1.;
  for(;;){
    for(i=0;i<n;i++) x[i] = xold[i]+alam*p[i];
    f = (*func)(x);
    if(alam < alamin){
      for(i=0;i<n;i++) x[i] = xold[i];
      check = 1;
      return;
    }else if(f <= fold+ALF*alam*slope) return;
    else{
      if(alam == 1.)
	tmplam = -slope/(2.*(f-fold-slope));
      else{
	rhs1 = f-fold-alam*slope;
	rhs2 = f2-fold2-alam2*slope;
	a=(rhs1/(alam*alam)-rhs2/(alam2*alam2))/(alam-alam2);
	b=(-alam2*rhs1/(alam*alam)+alam*rhs2/(alam2*alam2))/(alam-alam2);
	if(a == 0.)
	  tmplam = -slope/(2.*b);
	else{
	  disc = b*b-3.*a*slope;
	  if(disc < 0.){
	    std::cout << "Roundoff problem in lnsrch.\n";
	    exit(1);
	  }
	  else
	    tmplam = (-b+sqrt(disc))/(3.*a);
	}
	tmplam = min(tmplam, 0.5*alam);
      }
    }  
    alam2 = alam;
    f2    = f;
    fold2 = fold;
    alam  = max(tmplam, 0.1*alam);
  }
}

