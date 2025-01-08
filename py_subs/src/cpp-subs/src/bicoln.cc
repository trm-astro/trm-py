#include <cstdlib>
#include <cmath>
#include "trm/subs.h"

/**
 * factln returns the natural logarithm of factorial n (ln(n!)).
 * This uses a small amount of static storage to speed up the computation
 * for small n.
 * \param n the value of to calculate the factorial of.
 */

double Subs::bicolnctln(int n){
    static double a[101];
    if(n < 0)
	throw Subs_Error("Subs::factln: n = " + Subs::str(n) + " is out of range");
    if(n <= 1) return 0.;
    if(n <= 100) 
	return a[n] ? a[n] : (a[n] = gammln(n+1.0));
    else
	return gammln(n+1.0);
}



