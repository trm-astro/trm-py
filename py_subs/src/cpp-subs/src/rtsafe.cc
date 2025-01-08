#include "trm/subs.h"

/**
 * rtsafe is a Numerical Recipes-based routine to find roots
 * of a function using bisection or Newton-Raphson as appropriate.
 * \param func function object inherited from the abstract class Subs::RTfunc which declares the
 * function that must be defined.
 * \param x1 value to the left of the root
 * \param x2 value to the right of the root
 * \param xacc minimum accuracy in returned root
 * \return Returns the x value of the root.
 */

double Subs::rtsafe(const RTfunc& func, double x1, double x2, double xacc){
    
  const int MAXIT = 100;
  int j;
  double xl, xh, fl, fh, df, rts, temp, dx, dxold, f;
  
  func(x1, fl, df);
  func(x2, fh, df);
  
  if((fl > 0. && fh > 0.) || (fl < 0. && fh < 0.))
    throw Subs_Error("double Subs::rtsafe(const RTfunc&, double, double, double): root not bracketed. x1,x2,fl,fh = " +
		     Subs::str(x1) + ", " + Subs::str(x2) + ", " + Subs::str(fl) + ", " + Subs::str(fh) );
  
  if(fl == 0.0) return x1;
  if(fh == 0.0) return x2;
  
  if(fl < 0.0){
    xl = x1;
    xh = x2;
  }else{
    xh = x1;
    xl = x2;
  }
  
  rts = 0.5*(x1+x2);
  dx  = dxold = fabs(x2-x1);
  func(rts, f, df);
  for(j=0;j<MAXIT;j++){
    if((((rts-xh)*df-f)*((rts-xl)*df-f) >= 0.0)
       || (fabs(2.0*f) > fabs(dxold*df))){
      dxold = dx;
      dx = 0.5*(xh-xl);
      rts = xl + dx;
      if(xl == rts) return rts;
    }else{
      dxold = dx;
      dx = f/df;
      temp = rts;
      rts -= dx;
      if(temp == rts) return rts;
    }
    if(fabs(dx) < xacc) return rts;
    func(rts,f,df);
    
    if(f < 0.0)
      xl = rts;
    else
      xh = rts;
  }
  throw Subs_Error("double Subs::rtsafe(const RTfunc&, double, double, double): too many iterations.");
}

