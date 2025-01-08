#include "trm/subs.h"

class Zaghloul {
    // function object to compute integrand of zaghloul equation 4
public:
    Zaghloul(double a, double v) : a(a), v(v) {}

    double operator()(double u) const {
	double diff = v-u;
	return std::exp(-diff*(v+u))*std::sin(2.0*a*diff);
    }
    
private:
    double a, v;
};

/*
 * Computes Voigt function following method of Zaghloul & Ali 2012
 */

double Subs::voigt(double a, double v, double eps) {
    
    const int IMAX = 10000;

    double ulo, uhi, t1, ff, gg, hh, t2;
    double sum, half_period;

    // for the expansion when a is greater than 26.6
    const double C2  = 1.0/2.0;
    const double C4  = 3.0/4.0;
    const double C6  = 15.0/8.0;
    const double C8  = 105.0/16.0;  
    const double C10 = 945.0/32.0; 
    const double C12 = 10395.0/64.0;
    
    // RTPI = sqrt(PI), FWHM = 2*sqrt(ln(2)) = 
    // full width half max of gaussian with variance 1/RTPI
    const double PI     = 3.1415926535897932384;
    const double RTPI   = 1.7724538509055159;
    const double FACTOR = 2.0/RTPI;
    const double FWHM   = 1.6651;
    const double ZERO   = 0.;
    const double EPS2   = 1.0e-4*eps;

    // The first term of zaghloul's equation 4, 
    if (a < 26.6) {
	t1 = std::exp(a*a) * erfc(a); 
    }else{
	double ainv2 = 1.0/(a*a);
	t1 = 1.0/(RTPI*a)*(1.0-ainv2*(C2-ainv2*(C4-ainv2*(C6-ainv2*(C8-ainv2*(C10-ainv2*C12))))));
    }

    ff  = std::exp(-(v*v));
    gg  = std::cos(2.0*a*v);
    hh  = std::sin(2.0*a*v);

    double ans = t1 * ff * gg;

    // the second term of zaghloul's eq 4
    // the sinusoids in the integrands have a period of pi/a_ratio.
    // we'll be integrating over each half-period so that all
    // contributions have the same sign if the half-period.
    // we'll also be integrating over the full-width-half-max
    // of the gaussian. we'll take the smallest length scale
    // as the one to integrate over.

    half_period = std::min(0.5*PI/a, FWHM);

    // first the voigt function itself
    // initialize the integral sum and the upper integration limit
    sum  = 0.0;
    uhi  = v;

    // form a lower integration bound, 
    // allow for positive of negative values of v

    // Declare function object for integration
    Zaghloul zaghloul(a,v);

    for(int i=0; i<IMAX; i++){

	if(uhi < 0.){
	    ulo = uhi + half_period;
	}else{
	    ulo = uhi - half_period;
	}

	if (uhi == 0.0){
	    ulo = 0.0;
	}else if(uhi < 0.0){
	    ulo = std::min(ZERO,ulo);
	}else{
	    ulo = std::max(ZERO,ulo);
	}

	// add in these contributions to the total
	t2   = Subs::qromb(zaghloul, ulo, uhi, eps, 5, 20, false); 
	sum += t2;
       
	// convergence check
	if (ulo == 0.0 || std::abs(t2/sum) <= EPS2) break;

	// swap the limits for the next half period, and end of loop
	uhi = ulo;

	if(i == IMAX - 1)
	    std::cerr << "WARNING: voigt did not converge" << std::endl;
    }

    return ans + FACTOR*sum;
}


// overload with a vectorized version
void Subs::voigt(double a, const double* v, double* out, int n, double eps) {
    for(int i=0; i<n; i++){
        out[i] = voigt(a, v[i], eps);
    }
}
