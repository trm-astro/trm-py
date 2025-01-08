#include <cmath>
#include "trm/subs.h"

/** Routine for bracketting the minimum of a function. Starts
 * from two points then carries on in a downhill direction until it
 * brackets a minimum.
 * \param ax   Input and returned. Initial position to define where to start from
 * \param bx   Input and returned. Another initial position to define where to start from. Head off in a downhill direction from these.
 * \param cx   Return value; on completion ax, bx, cx will bracket the minimum
 * \param fa   function value at ax on return
 * \param fb   function value at bx on return
 * \param fc   function value at cx on return
 * \param func the function object whose minimum is to be bracketted.
 */

namespace Subs {
  void shift(double& a, double& b, double& c, double d);
}

void Subs::mnbrak(double& ax, double& bx, double& cx, double& fa, double& fb, double& fc, Sfunc& func){

  double ulim, ux, r, q, fu;

  const double GOLD   = 1.618034;
  const double GLIMIT = 100.;
  const double TINY   = 1.e-20;
  
  fa = func(ax);
  fb = func(bx);
  if(fb > fa){
    std::swap(ax,bx);
    std::swap(fa,fb);
  }

  // First guess at c
  cx = bx + GOLD*(bx-ax);
  fc = func(cx);

  while(fb > fc){

    r    = (bx-ax)*(fb-fc);
    q    = (bx-cx)*(fb-fa);
    ux   = bx - ((bx-cx)*q-(bx-ax)*r)/(2.*sign(std::max(abs(q-r),TINY),q-r));
    ulim = bx + GLIMIT*(cx-bx);

    if((bx-ux)*(ux-cx) > 0.){
      fu = func(ux);
      if(fu < fc){
	ax = bx;
	bx = ux;
	fa = fb;
	fb = fu;
	return;
      }else if(fu > fb){
	cx = ux;
	fc = fu;
	return;
      }
      ux = cx + GOLD*(cx-bx);
      fu = func(ux);

    }else if((cx-ux)*(ux-ulim) > 0.){
      fu = func(ux);
      if(fu < fc){
	shift(bx, cx, ux, cx+GOLD*(ux-bx));
	shift(fb, fc, fu, func(ux));
      }

    }else if( (ux-ulim)*(ulim-cx) >= 0.){
      ux = ulim;
      fu = func(ux);

    }else{
      ux = cx + GOLD*(cx-bx);
      fu = func(ux);
    }

    shift(ax, bx, cx, ux);
    shift(fa, fb, fc, fu);
  }
}

namespace Subs {
  void shift(double& a, double& b, double& c, double d){
    a = b;
    b = c;
    c = d;
  }
}
