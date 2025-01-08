#include "trm/subs.h"
#include "trm/array1d.h"

/** Default constructor
 * \param func the multi-D function
 * \param p the starting point
 * \param d the direction to take
 */
Subs::Safunc::Safunc(Afunc& func, const Array1D<double>& p, const Array1D<double>& d, const Array1D<double>& scale) : 
  func(func), p(p), d(d), scale(scale) {}

/** This routine should return the function value after it has been passed 
 * a vector of parameters. It is not defined as 'const' since you may well want to 
 * alter the object during the computation.
 */
double Subs::Safunc::operator()(double x) {
  Array1D<double> r = p + x*scale*d;
  return func(r);
}
