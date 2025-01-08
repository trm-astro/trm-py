#include "trm/subs.h"

/** brent implements NR 1D function minimisation method. Be careful! The order of
 * the arguments has been swaaped from the NR version and it is important to
 * get them right.
 *
 * \param xstart initial guess at minimum
 * \param x1     first point bracketing minimum
 * \param x2     last point bracketing minimum
 * \param f      a function object inherited from Subs::Sfunc
 * which returns the value of the function to be minimised at any position.
 * \param tol    tolerance in x to locate minimum
 * \param xmin   the returned minimum
 * \return brent returns the minimum value
 */
double Subs::brent(double xstart, double x1, double x2, Sfunc& f, double tol, double& xmin){
  
  const int ITMAX = 100;
  const double CGOLD = 0.3819660;
  
  int iter;
  double a,b,d=0,etemp,fu,fv,fw,fx,p,q,r,tol2,u,v,w,x,xm;
  double e=0.0;
  
  a = x1 < x2 ? x1 : x2;
  b = x1 > x2 ? x1 : x2;
  x  = w  = v  = xstart;
  fw = fv = fx = f(x);
  
  tol2 = 2.*tol;
  for(iter=0;iter<ITMAX;iter++){
    xm   = 0.5*(a+b);
    if(fabs(x-xm) < (tol2-0.5*(b-a))){
      xmin = x;
      return fx;
    }
    if(fabs(e) > tol){
      r = (x-w)*(fx-fv);
      q = (x-v)*(fx-fw);
      p = (x-v)*q-(x-w)*r;
      q = 2.0*(q-r);
      if(q > 0.0) p = -p;
      q     = fabs(q);
      etemp = e;
      e     = d;
      if(fabs(p) >= fabs(0.5*q*etemp) || p <= q*(a-x) || p >= q*(b-x))
	d = CGOLD*(e=(x >= xm ? a-x : b-x));
      else{
	d = p/q;
	u = x+d;
	if(u-a < tol2 || b-u < tol2) d = sign(tol,xm-x);
      }
    }else{
      d=CGOLD*(e=(x >= xm ? a-x : b-x));
    }
    u  = fabs(d) >= tol ? x+d : x+sign(tol,d);
    
    fu = f(u); // the one function evaluation/iteration
    
    if(fu <= fx){
      if(u >= x) a=x; else b=x;
      v  = w;
      w  = x;
      x  = u;
      fv = fw;
      fw = fx;
      fx = fu;
    }else{
      if(u < x) a=u; else b=u;
      if(fu <= fw || w==x){
	v  = w;
	w  = u;
	fv = fw;
	fw = fu;
      }else if(fu <= fv || v == x || v == w){
	v  = u;
	fv = fu;
      }
    }
  }
  std::cerr << "Too many iterations in brent" << std::endl;
  xmin = x;
  return fx;
}
