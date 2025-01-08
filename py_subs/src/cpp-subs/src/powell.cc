#include "trm/subs.h"
#include "trm/format.h"
#include "trm/array1d.h"

namespace Subs {
    void linmin(Subs::Array1D<double>& p, Subs::Array1D<double>& d, const Subs::Array1D<double>& scale, 
		Subs::Afunc& func, double& fret);
}

/** Implementation of Powell's method of function minimisation.
 * 
 * \param p    initial position, set to the best point found on output
 * \param xi   initial directions to minimise along, returned as current best set of directions. These directions should N
 * Array1Ds, each of N elements where N equals number of variable parameters. Initialise as a set of unit vectors.
 * \param scale factors showing the typical scale length of motion in each parameter, needed to help the routine if they are
 * very different. This is a modification of the NR routine because I found that mixing very different scales screws the line 
 * minimisation up.
 * \param ftol fractional tolerance on function value to indicate done minimising
 * \param itmax maximum number of iterations to perform
 * \param iter number of iterations carried out
 * \param fret returned function value 
 * \param func the function to minimise supports a call to an Array1D of parameters.
 */

void Subs::powell(Array1D<double>& p, Buffer1D<Array1D<double> >& xi, const Array1D<double>& scale, double ftol, int itmax, int& iter, double& fret, Afunc& func){

  const double TINY = 1.e-25;

  // calculate function value at starting point
  fret = func(p);

  // Save the starting point
  Array1D<double> pt = p, ptt(p.size()), d(p.size());

  int ibig;
  double fp, fptt, del;
  for(iter=0; iter<itmax; iter++){
    fp   = fret;
    ibig = 0;
    del  = 0.;

    // minimize along each direction in turn
    for(int i=0; i<p.size(); i++){
      fptt = fret;
      d    = xi[i];
      linmin(p, d, scale, func, fret);
      if(fptt - fret > del){
	del  = fptt - fret;
	ibig = i;
      }
    }
     
    if(2.*(fp-fret) <= ftol*(fabs(fp)+fabs(fret))+TINY)
      return;

    // Construct extrapolated point and the average direction moved
    // Save old starting point
    ptt = 2.*p - pt;
    d   = (p - pt)/scale;
    pt  = p;

    fptt = func(ptt);

    if(fptt < fp){
      double t = 2.*(fp-2.*fret+fptt)*sqr(fp-fret-del)-del*sqr(fp-fptt);
      if(t < 0.){
	linmin(p, d, scale, func, fret);
	xi[ibig]       = xi[p.size()-1];
	xi[p.size()-1] = d;
      }
    }
  }    
}

namespace Subs {
    // this minimises along a line

    void linmin(Subs::Array1D<double>& p, Subs::Array1D<double>& d, const Subs::Array1D<double>& scale, Subs::Afunc& func, double& fret){
  
	Subs::Safunc oned(func, p, d, scale);
	double ax = 0., bx = 1., cx, fa, fb, fc;
	Subs::mnbrak(ax, bx, cx, fa, fb, fc, oned);
	double xmin;
	fret = brent(bx, ax, cx, oned, 1.e-6, xmin);
	d *= xmin;
	p += scale*d;
	
    }
   
}
