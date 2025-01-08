#include "trm/subs.h"

/* Numerical differentiation (Stu Littlefair). Takes centred differences for all but the end points. 
 * i.e. derivative of point n = (y[n+1]-y[n-1])/(x[n+1]-x[n-1]). The end points are 
 * set to the forward and backward versions of this, e.g. (y[1]-y[0])/(x[1]-x[0])
 * \param x - the x data
 * \param y - the y data
 * \param deriv - the returned y values
 */
void Subs::numdiff(const Buffer1D<double>& x, const Buffer1D<double>& y, Buffer1D<double>& deriv){

  if (x.size() != y.size())
    throw Subs_Error("Subs::numdiff: x and y data sizes must match");

  if (x.size() < 2)
    throw Subs_Error("Subs::numdiff: must have at least 2 points");

  deriv.resize(x.size());
  int np     = x.size();
  
  // Calculate differential for all but first and last points
  for(int j=1; j<np-1; j++)
    deriv[j] = (y[j+1]-y[j-1])/(x[j+1]-x[j-1]);

  // now set first and last points
  deriv[0]    = (y[1]-y[0])/(x[1]-x[0]);
  deriv[np-1] = (y[np-1]-y[np-2])/(x[np-1]-x[np-2]);

}


