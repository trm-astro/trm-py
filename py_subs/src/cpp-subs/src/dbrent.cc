#include <cmath>
#include <cstdlib>
#include <iostream>
#include "trm/subs.h"

/* Given a bracketted minimum, this routine refines the minimum, making use of derivatives. 
 * It comes from NR, with the added option of a quick bail out. ax to cx must bracket a 
 * minimum otherwise you will be in trouble.
 * \param ax       one extreme of x
 * \param bx       mid x value
 * \param cx       other extreme x value
 * \param func     function object returning function value at x
 * \param dfunc    function object retunring derivative value at x
 * \param acc      (absolute) accuracy in x
 * \param stopfast retrun as soon as function value goes below a reference value or not
 * \param fref     reference function value if stopfast = true
 * \param xmin     x at minimum
 * \return function value at xmin
 */

double Subs::dbrent(double ax, double bx, double cx, Sfunc& func, Sfunc& dfunc, double acc,
		     bool stopfast, double fref, double& xmin){

  double a = (ax < cx ? ax : cx);  
  double b = (ax > cx ? ax : cx);
  double x, w, v;
  x = w = v = bx;
  double fw, fv, fx;
  fw = fv = fx = func(x);
  
  if(stopfast && fx < fref){
    xmin = x;
    return fx;
  }

  double dw, dv, dx;
  dw = dv = dx = dfunc(x);
  
  const int ITMAX=100;

  double xm, tol1, tol2, e=0., d1, d2, u1, u2, olde, d = 0, u, fu, du;
  bool ok1, ok2;
  for(int iter=0; iter<ITMAX; iter++){
    xm = 0.5*(a+b);
    tol1 = acc;
    tol2 = 2.*tol1;
    if(fabs(x-xm) <= (tol2 -0.5*(b-a))){
      xmin = x;
      return fx;
    }
    if(fabs(e) > tol1){
      d2 = d1 = 2.*(b-a);
      if(dw != dx) d1 = (w-x)*dx/(dx-dw);
      if(dv != dx) d2 = (v-x)*dx/(dx-dv);

      u1  = x + d1;
      u2  = x + d2;
      ok1 = (a-u1)*(u1-b) > 0. && dx*d1 <= 0.;
      ok2 = (a-u2)*(u2-b) > 0. && dx*d2 <= 0.;
      olde = e;
      e = d;
      if(ok1 || ok2){
	if(ok1 && ok2)
	  d = (fabs(d1) < fabs(d2) ? d1 : d2);
	else if(ok1)
	  d = d1;
	else
	  d = d2;
	if(fabs(d) <= fabs(0.5*olde)){
	  u = x + d;
	  if(u - a < tol2 || b-u < tol2)
	    d = Subs::sign(tol1, xm-x);
	}else{
	  d = 0.5*(e = (dx >= 0. ? a-x : b-x));
	}
      }else{
	d = 0.5*(e = (dx >= 0. ? a-x : b-x));
      }
    }else{
      d = 0.5*(e = (dx >= 0. ? a-x : b-x));
    }
      
    if(fabs(d) >= tol1){
      u  = x + d;
      fu = func(u);
      if(stopfast && fu < fref){
	xmin = u;
	return fu;
      }
    }else{
      u = x + Subs::sign(tol1,d);
      fu = func(u);
      if(stopfast && fu < fref){
	xmin = u;
	return fu;
      }
      if(fu > fx){
	xmin = x;
	return fx;
      }
    }
    du = dfunc(u);
    if(fu <= fx){
      if(u >= x) 
	a=x;
      else
	b=x;
      v=w; fv=fw; dv=dw;	
      w=x; fw=fx; dw=dx;
      x=u; fx=fu; dx=du;
    }else{
      if(u < x) 
	a=u;
      else
	b=u;
      if(fu <= fw || w == x){
	v=w; fv=fw; dv=dw;	
	w=u; fw=fu; dw=du;
      }else if(fu < fv || v == x || v == w){
	v=u; fv=fu; dv=du;
      }
    }
  }
  throw Subs_Error("Subs::dbrent: too many iterations");
  return 0.;
}
