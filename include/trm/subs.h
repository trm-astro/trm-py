#ifndef TRM_SUBS_H
#define TRM_SUBS_H

#include <cstdlib>
#include <cmath>
#include <stdint.h>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <sstream>
#include <list>
#include <vector>
#include "plstream.h"

//! Namespace for workhorse functions

/**
 * This namespace contains generally useful routines and classes, such as sorting routines and the like.
 */

namespace Subs {

    // Start with some typedefs to define particular lengths of data as needed for
    // reading/writing binary data from disk in a reasonably machine-independent manner
    // and for some particular routines (ran4). These may need altering on different machines

    //! 1-byte signed integer
    typedef char CHAR;

    //! 1-byte unsigned integer
    typedef unsigned char UCHAR;

    //! 2-byte signed integer
    typedef int16_t INT2;

    //! 2-byte unsigned integer
    typedef uint16_t UINT2;

    //! 4-byte signed integer
    typedef int32_t INT4;

    //! 4-byte unsigned integer
    typedef uint32_t UINT4;

    //! 4-byte floating point number
    typedef float REAL4;

    //! 8-byte floating point number
    typedef double REAL8;

    // Forward declarations
    template <class T> class Array1D;
    template <class T> class Array2D;
    template <class T> class Buffer1D;
    template <class T> class Buffer2D;
    class Poly;

    //! Error class for any of the routines in subs to inherit.

    /**
     * Subs_Error objects are inherited from the string class. They can
     * therefore be printed as usual.
     */

    class Subs_Error : public std::string {
    public:

	//! Default constructor
	Subs_Error() : std::string() {}

	//! Constructor which stores a message
	Subs_Error(const std::string& str) : std::string(str) {}
    };

    template <class X, class Y, class Z>
    struct xyz;

    template <class X, class Y, class Z>
    std::istream& operator>>(std::istream& ist, xyz<X,Y,Z>& obj);

    template <class X, class Y, class Z>
    std::ostream& operator<<(std::ostream& ost, const xyz<X,Y,Z>& obj);

    //! General 3 variable structure

    /** This is completely general structure for storing triples of numbers,
     * for instance time, velocity, error or whatever.
     */
    template <class X, class Y, class Z>
    struct xyz{

	//! Default constructor
	xyz() : x(0), y(0), z(0) {}

	//! Constructor from values
	xyz(const X& xx, const Y& yy, const Z& zz) : x(xx), y(yy), z(zz) {}

	//! X value
	X x;

	//! Y value

	Y y;

	//! Z value
	Z z;

	//! ASCII input
	friend std::istream& operator>><>(std::istream& ist, xyz<X,Y,Z>& obj);

	//! ASCII output
	friend std::ostream& operator<<<>(std::ostream& ost, const xyz<X,Y,Z>& obj);
    };
  
    /** ASCII input operator. Reads three numbers separated by spaces.
     * \param ist input stream
     * \param obj the 3 parameter object to read the data into
     * \return the input stream
     */
    template <class X, class Y, class Z>
    std::istream& operator>>(std::istream& ist, xyz<X,Y,Z>& obj){
	ist >> obj.x >> obj.y >> obj.z;
	return ist;
    }
  
    /** ASCII output operator. Writes three numbers separated by spaces.
     * \param ost output stream
     * \param obj the 3 parameter object to write out
     * \return the output stream
     */
    template <class X, class Y, class Z>
    std::ostream& operator<<(std::ostream& ost, const xyz<X,Y,Z>& obj){
	ost << obj.x << " " << obj.y << " " << obj.z;
	return ost;
    }

    //! Defines a particular instance of an xyz<X,Y,Z> object suited to radial velocity work
    typedef xyz<double,float,float> rv;

    //! Defines a particular instance of an xyz<X,Y,Z> object suited for more accurate work
    typedef xyz<double,double,float> ddat;

    template <class X, class Y>
    struct xy;

    template <class X, class Y>
    std::istream& operator>>(std::istream& ist, xy<X,Y>& obj);

    template <class X, class Y>
    std::ostream& operator<<(std::ostream& ost, const xy<X,Y>& obj);

    //! General 2 variable structure

    /** This is completely general structure for storing pairs of numbers,
     * for instance time, velocity, error or whatever.
     */
    template <class X, class Y>
    struct xy{

	//! Default constructor
	xy() : x(0), y(0) {}

	//! Constructor from values
	xy(const X& xx, const Y& yy) : x(xx), y(yy) {}

	//! X value
	X x;

	//! Y value
	Y y;

	//! ASCII input
	friend std::istream& operator>><>(std::istream& ist, xy<X,Y>& obj);

	//! ASCII output
	friend std::ostream& operator<<<>(std::ostream& ost, const xy<X,Y>& obj);
    };

    /** ASCII input operator. Reads two numbers separated by spaces.
     * \param ist input stream
     * \param obj the 2 parameter object to read the data into
     * \return the input stream
     */
    template <class X, class Y>
    std::istream& operator>>(std::istream& ist, xy<X,Y>& obj){
	ist >> obj.x >> obj.y;
	return ist;
    }

    /** ASCII output operator. Writes two numbers separated by spaces.
     * \param ost output stream
     * \param obj the 2 parameter object to write out
     * \return the output stream
     */
    template <class X, class Y>
    std::ostream& operator<<(std::ostream& ost, const xy<X,Y>& obj){
	ost << obj.x << " " << obj.y;
	return ost;
    }

    //! Square a value
    template <class X> 
    inline X sqr(const X& a){
	return (a*a);
    }

    //! Take the nearest integer
    template <class X> 
    inline X nint(const X& a){
	return X(floor(a+0.5));
    }

    //! Return a value clamped to lie between lower and upper values
    template <class X> 
    X clamp(const X& low, const X& value, const X& high){
	if(low > high) 
	    throw Subs_Error("limit(X&, X&, X&): low value is greater than high value");
	if(value < low)  return low;
	if(value > high) return high;
	return value;
    }

    //! Return a with the same sign as b
    template <class X, class Y> 
    inline X sign(const X& a, const Y& b){
	return(b >= 0. ? (X)fabs(a) : -(X)fabs(a));
    } 

    //! Return the absolute value
    template <class X> 
    inline X abs(const X& a){
	if(a < 0)
	    return (-a);
	else
	    return a;
    }

    //! Compute sqrt(a*a+b*b) avoiding under/over flows 
    template <class X> X pythag(const X& a, const X& b){
      X absa, absb, temp;
      absa = fabs(a);
      absb = fabs(b);
      if(absa > absb){
	temp = absb / absa;
	return absa*sqrt(1.+temp*temp);
      }else if(absb == 0.){
	return 0.;
      }else{
	temp = absa / absb;
	return absb*sqrt(1.+temp*temp);
      }
    }

    //! Structure for return of useful observing info.
    struct Altaz{
	//! hour angle (hours after meridian)
	double ha;       
	//! altitude above horizon, degrees (no refraction)
	double alt_true; 
	//! altitude above horizon, degrees (including refraction)
	double alt_obs;  
	//! azimuth from north through east
	double az;       
	//! column of air relative to zenith at sea-level.
	double airmass;
	//! parallactic angle
	double pa;       
    };

    //! Computes the centroid of a peak in a 1D array
    void centroid(const float data[], const float var[], int p1, int p2, float fwhm, float start, bool emission, double& pos, float& epos);

    //! Generates uniform random deviates
    double ran1(INT4 &seed);

    //! Generates uniform random deviates
    double ran2(INT4 &seed);

    //! Generates uniform random deviates
    double ran3(INT4 &seed);

    //! Generates uniform random deviates
    double ran4(INT4 &seed);

    //! Generates gaussian random deviates
    double gauss1(INT4 &seed);

    //! Generates gaussian random deviates
    double gauss2(INT4 &seed);

    //! Generates gaussian random deviates
    double gauss3(INT4 &seed);

    //! Generates gaussian random deviates
    double gauss4(INT4 &seed);

    //! Poisson deviates
    float poisson1(float mu, INT4& seed);

    //! Poisson deviates
    float poisson2(float mu, INT4& seed);

    //! Poisson deviates
    float poisson3(float mu, INT4& seed);

    //! Poisson deviates
    float poisson4(float mu, INT4& seed);

    //! Adds extension onto name if not already present
    std::string filnam(std::string name, const std::string& extens);

    //! Trapezoidal intgration routine
    float trapzd(float (*func)(float x), float a, float b, int n);

    //! Put a program to sleep
    void sleep(double seconds);

    //! Convert degrees to radians
    inline double deg2rad(double deg){
      return M_PI*deg/180.;
    }

    //! Convert radians to degrees
    inline double rad2deg(double rad){
      return 180.*rad/M_PI;
    }

    //! Simpson's rule integration routine
    float qsimp(float (*func)(float x), float a, float b);

    //! Singular value decomposition
    /**
     * svdcmp performs singular value decomposition.
     * Given an M x N matrix A this routine computes its
     * singular value decomposition A = U.W.V^t. The matrix U
     * replaces A on output. The diagonal matrix W is returned
     * as a vector and the matrix V is returned in the last argument (not
     * the transpose).
     * \param a M x N matrix A. Declare e.g. a(M,N). M might be the number of data point and N the
     * number of polynomials for example. Returns with elements of the matrix U.
     * \param w N element vector of diagonal elements of centre matrix W
     * \param v N by N elements of matrix V
     */
    template <class X>
    void svdcmp(Buffer2D<X>& a,  Buffer1D<X>& w, Buffer2D<X>& v){
  
	int m = a.nrow(), n = a.ncol();
	w.resize(n);
	v.resize(n,n);
	int flag, i, its, j, jj, k, l=0, nm=0;
	X anorm, c, f, g, h, s, scale, x, y, z;
    
	Buffer1D<X> rv1(n); // Work space array
	g = scale = anorm = 0.;
    
	for(i=0; i<n; i++){
	    l = i+1;
	    rv1[i] = scale*g;
	    g=s=scale=0.;
	    if(i < m){
		for(k=i; k<m;k++) scale += fabs(a[k][i]);
		if(scale){
		    for(k=i;k<m;k++){
			a[k][i] /= scale;
			s += a[k][i]*a[k][i];
		    }
		    f = a[i][i];
		    g = -sign(sqrt(s),f);
		    h = f*g-s;
		    a[i][i] = f-g;
		    for(j=l;j<n;j++){
			for(s=0.,k=i;k<m;k++) s += a[k][i]*a[k][j];
			f=s/h;
			for(k=i;k<m;k++) a[k][j] += f*a[k][i];
		    }
		    for(k=i;k<m;k++) a[k][i] *= scale;
		}
	    }
	    w[i] = scale*g;
	    g=s=scale=0.;
	    if(i<m && i != n-1){
		for(k=l;k<n;k++) scale += fabs(a[i][k]);
		if(scale){
		    for(k=l;k<n;k++){
			a[i][k] /= scale;
			s += a[i][k]*a[i][k];
		    }
		    f=a[i][l];
		    g=-sign(sqrt(s),f);
		    h=f*g-s;
		    a[i][l]=f-g;
		    for(k=l;k<n;k++) rv1[k]=a[i][k]/h;
		    for(j=l;j<m;j++){
			for(s=0.,k=l;k<n;k++) s += a[j][k]*a[i][k];
			for(k=l;k<n;k++) a[j][k] += s*rv1[k];
		    }
		    for(k=l;k<n;k++) a[i][k] *= scale;
		}
	    }
	    anorm = std::max(anorm, X(fabs(w[i])+fabs(rv1[i])));
	}
	for(i=n-1;i>=0;i--){
	    if(i < n-1){
		if(g){
		    for(j=l;j<n;j++)
			v[j][i] = (a[i][j]/a[i][l])/g;
		    for(j=l;j<n;j++){
			for(s=0.,k=l;k<n;k++) s += a[i][k]*v[k][j];
			for(k=l;k<n;k++) v[k][j] += s*v[k][i];
		    }
		}
		for(j=l;j<n;j++) v[i][j]=v[j][i]=0.;
	    }
	    v[i][i]=1.;
	    g=rv1[i];
	    l=i;
	}
	for(i=std::min(m,n)-1;i>=0;i--){
	    l=i+1;
	    g=w[i];
	    for(j=l;j<n;j++) a[i][j] = 0.;
	    if(g){
		g=1./g;
		for(j=l;j<n;j++){
		    for(s=0.,k=l;k<m;k++) s += a[k][i]*a[k][j];
		    f=(s/a[i][i])*g;
		    for(k=i;k<m;k++) a[k][j] += f*a[k][i];
		}
		for(j=i;j<m;j++) a[j][i] *= g;
	    }else for(j=i;j<m;j++) a[j][i] = 0.;
	    ++a[i][i];
	}
	for(k=n-1;k>=0;k--){
	    for(its=1;its<=30;its++){
		flag=1;
		for(l=k;l>=0;l--){
		    nm=l-1;
		    if(X(fabs(rv1[l])+anorm) == anorm){
			flag=0;
			break;
		    }
		    if(X(fabs(w[nm])+anorm) == anorm) break;
		}
		if(flag){
		    c=0.;
		    s=1.;
		    for(i=l;i<=k;i++){
			f=s*rv1[i];
			rv1[i] *= c;
			if(X(fabs(f)+anorm) == anorm) break;
			g=w[i];
			h=pythag(f,g);
			w[i]=h;
			h=1./h;
			c=g*h;
			s = -f*h;
			for(j=0;j<m;j++){
			    y=a[j][nm];
			    z=a[j][i];
			    a[j][nm]=y*c+z*s;
			    a[j][i]=z*c-y*s;
			}
		    }
		}
		z=w[k];
		if(l==k){
		    if(z < 0.){
			w[k] = -z;
			for(j=0;j<n;j++) v[j][k] = -v[j][k];
		    }
		    break;
		}
		if(its == 30) 
		    std::cerr << "No convergence in svdcmp in 30 iterations\n";
		x=w[l];
		nm=k-1;
		y=w[nm];
		g=rv1[nm];
		h=rv1[k];
		f=((y-z)*(y+z)+(g-h)*(g+h))/(2.*h*y);
		g=pythag(f,X(1));
		f=((x-z)*(x+z)+h*((y/(f+sign(g,f)))-h))/x;
		c=s=1.;
		for(j=l;j<=nm;j++){
		    i=j+1;
		    g=rv1[i];
		    y=w[i];
		    h=s*g;
		    g=c*g;
		    z=pythag(f,h);
		    rv1[j]=z;
		    c=f/z;
		    s=h/z;
		    f=x*c+g*s;
		    g=g*c-x*s;
		    h=y*s;
		    y *= c;
		    for(jj=0;jj<n;jj++){
			x=v[jj][j];
			z=v[jj][i];
			v[jj][j]=x*c+z*s;
			v[jj][i]=z*c-x*s;
		    }
		    z=pythag(f,h);
		    w[j]=z;
		    if(z){
			z=1./z;
			c=f*z;
			s=h*z;
		    }
		    f=c*g+s*y;
		    x=c*y-s*g;
		    for(jj=0;jj<m;jj++){
			y=a[jj][j];
			z=a[jj][i];
			a[jj][j]=y*c+z*s;
			a[jj][i]=z*c-y*s;
		    }
		}
		rv1[l]=0.;
		rv1[k]=f;
		w[k]=x;
	    }
	}
    }
	
    //! Singular value decomposition, back substitution
    /** Solves u x = b for vectors x (N elements) and b (M elements),
     *  where M x N matrix u has been transformed into its singular value 
     *  decomposition by svdcmp
     *
     * \param u M x N matrix from svdcmp. Declare e.g. u(M,N). M might be the number of data point and N the number of polynomials for example.
     * \param w N element vector from svdcmp. Should have been edited to set small values = 0
     * \param v N x N matrix from svdcmp
     * \param b M element target vector
     * \param x  element solution vector
     */
    template <class X>
    void svbksb(const Buffer2D<X>& u, const Buffer1D<X>& w, const Buffer2D<X>& v, const Buffer1D<X>& b, Buffer1D<X>& x){
  
	size_t jj, j, i, m = u.nrow(), n = u.ncol();
	X s;
    
	x.resize(n);
	Buffer1D<X> tmp(n); // Work space array
    
	for(j=0;j<n;j++){
	    s=0.0;
	    if(w[j]){
		for(i=0;i<m;i++) s += u[i][j]*b[i];
		s /= w[j];
	    }
	    tmp[j] = s;
	}
    
	for(j=0;j<n;j++){
	    s=0.0;
	    for(jj=0;jj<n;jj++) s += v[j][jj]*tmp[jj];
	    x[j]=s;
	}
    }
  
    //! Singular value decomposition fitting
    double svdfit(const Buffer1D<rv>& data, Buffer1D<float>& a,
		  const Buffer2D<float>& vect, Buffer2D<float>& u,
		  Buffer2D<float>& v, Buffer1D<float>& w);

    //! Singular value decomposition fitting
    double svdfit(const Buffer1D<ddat>& data, Buffer1D<double>& a, const Buffer2D<double>& vect, 
		  Buffer2D<double>& u, Buffer2D<double>& v, Buffer1D<double>& w);

    //! Singular value decomposition fitting
    double svdfit(const Buffer1D<rv>& data, Buffer1D<float>& a,
		  const Buffer1D<double>& cosine, const Buffer1D<double>& sine, 
		  Buffer2D<float>& u, Buffer2D<float>& v, Buffer1D<float>& w);

    //! Evaluates incomplete gamma function
    void  gser(double &gamser, double a, double x, double &gln);

    //! Evaluates incomplete gamma function
    void  gcf(double &gammcf, double a, double x, double &gln);

    //! Evaluates incomplete gamma function
    double gammp(double a, double x);

    //! Evaluates incomplete gamma function
    double gammq(double a, double x); // standard incomplete gamma function on two doubles
	void gammq(double a, double* x, double* out, int n); // vectorized version of gammq on a as a double and x as a vector(double)
  
    //! Evaluates ln of the gamma function
    double gammln(double xx);

    //! Evaluates ln of factorial n
    double factln(int n);

    //! FFT routine of complex single precision array
    void  fft(float *data, unsigned long nump, int flag);

    //! FFT routine of real single precision array
    void  fftr(float *data, unsigned long nump, int flag);

    //! FFT of complex double precision array
    void  fft(double *data, unsigned long nump, int flag);

    //! FFT of real double precision array
    void  fftr(double *data, unsigned long nump, int flag);

    //! Simultaneous FFT of two real double precision arrays
    void twofft(double data1[], double data2[], double fft1[], double fft2[], unsigned long int n);

    //! Lomb-Scargle periodogram, Press & Rybicki's fast method.
    void fasper(double *x, float *y, float *e,  int n, double ofac, double hifac, Buffer1D<double>& freq, Buffer1D<double>& pgram);

    //! Lomb-Scargle periodogram, Press & Rybicki's fast method.
    void fasper(double *x, float *y, float *e, int n, double fmax, int nfreq, Buffer1D<double>& freq, Buffer1D<double>& pgram);

    //! Amplitude spectrum, Press & Rybicki-type method
    void famp(double *x, float *y, float *e,  int n, double fmax, int nfreq, Subs::Buffer1D<double>& freq, Subs::Buffer1D<double>& amps);

    //! Burlisch-Stoer routine
    bool bsstep(double y[], double dydx[], int nv, double& xx, 
		double htry, double eps, double yscal[], double &hdid, 
		double &hnext, 
		void (*derivs)(double, double [], double []));

    //! Modified mid-point routine
    void  mmid(double y[], double dydx[], int nvar, double xs, 
	       double htot, int nstep, double yout[], 
	       void (*derivs)(double, double[], double[]));

    //! Abstract class for bsstep function object
    class Bsfunc {

    public:

	//! The function call
	/** Evaluates derivatives dydt given current time t and coordinates y
	 * use function object to store other parameters needed
	 */
      virtual void operator()(double t, double y[], double dydt[]) const = 0;
      virtual ~Bsfunc(){}
    };

    bool bsstep(double y[], double dydx[], int nv, double &xx, 
		double htry, double eps, double yscal[], 
		double &hdid, double &hnext, const Bsfunc& derivs);

    //! Modified mid-point routine
    void  mmid(double y[], double dydx[], int nvar, double xs, 
	       double htot, int nstep, double yout[], 
	       const Bsfunc& derivs);

    // Alternative for mmid for conservative 2nd order equations
    bool bsstepst(double y[], double dydx[], int nv, double &xx, 
		  double htry, double eps, double yscal[], 
		  double &hdid, double &hnext, const Bsfunc& derivs);

    void stoerm(double y[], double d2y[], int nvar, double xs, double htot,
		int nstep, double yout[], const Bsfunc& derivs);

    //! Polynomial extrapolation routine
    void  pzextr(int iest, double xest, double yest[], double yz[], 
		 double dy[], int nv);

    //! Computes eigenvalues and vectors 
    void  jacob(Buffer2D<double>& a, Buffer1D<double>& d, Buffer2D<double>& v, int &nrot);

    //! Compute Chi**2 level corresponding to a certain chance
    double tchi(double lfprob, int ndof);

    //! LU decomposition routine
    void ludcmp(double **a, unsigned int n, unsigned *indx, double &d);

    //! LU decomposition routine
    void ludcmp(Buffer2D<double>& a, Buffer1D<size_t>& indx, double &d);

    //! LU decomposition back substitution routine
    void lubksb(double **a, unsigned int n, unsigned int *indx, double *b);

    //! LU decomposition back substitution routine
    void lubksb(const Buffer2D<double>& a, const Buffer1D<size_t>& indx, Buffer1D<double>& b);

    //! Gauss-Jordan elimination
    void gaussj(Buffer2D<double>& a, Buffer2D<double>& b);

    //! Gauss-Jordan elimination
    void gaussj(int nover, Buffer2D<double>& a, Buffer2D<double>& b);

    //! Gauss-Jordan elimination
    void gaussj(double** a, int n, double** b, int m);
  
    //! Sigma clipping
    void sigma_reject(const float* data, int n, float thresh, bool careful,
		      double& rawmean, double& rawrms, double& mean, double& rms, int& nrej);

    //! Numerical derivative
    void numdiff(const Buffer1D<double>& x, const Buffer1D<double>& y, Buffer1D<double>& deriv);

    //! Boxcar averaging
    void boxcar(const Buffer1D<double>& data, Buffer1D<double>& filt, int width);

    //! Astronomical extinction (Cardelli et al 1989)
    double extinct(double lambda, double ratio);

    //! Voigt function
    double voigt(double a, double v, double eps);

	//! Overloaded version of voigt function for vectorized operations
	void voigt(double a, const double* v, double* out, int n, double eps);

    //! Converts an integer to a char*
    void strint(const unsigned int num, const unsigned int nd, char* intstr);

    //! Makes a string from a numerical value
    template <class T>
	std::string str(const T& con, int prec=15){
	std::ostringstream ss;
	ss.precision(prec);
	if(!(ss << con))
	    throw Subs_Error("Subs::str: conversion error. Something must be very odd.");
	return ss.str();
    }

    //! Returns character equivalent to a single digit
    char digit_to_char(int digit);

    //! Converts a string to upper case
    std::string toupper(const std::string& str);

    //! Converts a string to upper case
    std::string tolower(const std::string& str);

    //! Makes an ndigit long string from an int
    std::string str(const int& con, int ndig);

    //! Makes an ndigit long std::string from a long int
    std::string str(const long int& con, int ndig);

    //! Writes a string to a binary file
    void write_string(std::ofstream& s, const std::string& str);

    //! Reads a string from a binary file
    void read_string(std::ifstream& s, std::string& str, bool swap_bytes);

    //! Skips a string while reading a binary file
    void skip_string(std::ifstream& s, bool swap_bytes);

    //! Reads a line from an input stream and splits into a vector of vectors of strings
    std::vector<std::vector<std::string> > read_multi_string_line(std::istream& s);

    //! Reads a line from an input stream and split into a vector of strings
    std::vector<std::string> read_line(std::istream& s);

    //! Strips off trailing trailing whitespace from a string
    std::string strip_trailing_whitespace(const std::string& str);

    //! Strips off trailing trailing whitespace from a string
    std::string strip_leading_whitespace(const std::string& str);

    //! Modes for rebin
    enum REBIN_MODE {
	AVERAGE,    /**< Normalises by the number of input pixels binned into each output pixel */
	INTEGRATE   /**< Sums up input pixels to form each output pixel */
    };

    //! Types of interpolation used by rebin
    enum REBIN_TYPE {
	LINEAR /**< Linear interpolation of flux in a bin */
    };

    //! Rebins data
    void rebin(const Buffer1D<double>& input, int infirst, int inlast, const Poly& inpoly, 
	       Buffer1D<double>& output, int outfirst, int outlast, const Poly& outpoly, 
	       REBIN_MODE mode, REBIN_TYPE type);

    //! Rebins data and variances
    void rebin(const Buffer1D<double>& indat, const Buffer1D<float>& invar, int infirst, int inlast, const Poly& inpoly, 
	       Buffer1D<double>& outdat, Buffer1D<float>& outvar, int outfirst, int outlast, const Poly& outpoly, 
	       REBIN_MODE mode, REBIN_TYPE type);

    //! Tests for little endian
    bool is_little_endian();

    //! Tests for big endian
    bool is_big_endian();

    //! Byte swaps an integer
    int byte_swap(int i);

    //! Byte swaps a unsigned integer
    unsigned int byte_swap(unsigned int i);

    //! Byte swaps a long unsigned integer
    long unsigned int byte_swap(long unsigned int i);

    //! Byte swaps a long integer
    long int byte_swap(long int i);

    //! Byte swaps a float
    float byte_swap(float f);

    //! Byte swaps a double
    double byte_swap(double f);

    //! Reverse the order of an array of bytes
    void reverse_bytes(char* buffer, const int& nbuff);

    //! Reverse 2 byte array
    void byte_swap2(char* buffer);

    //! Reverse 4 byte array
    void byte_swap4(char* buffer);

    //! Reverse 8 byte array
    void byte_swap8(char* buffer);

    //! Byte-swap array of data of arbitrary type
    /** Carries out a byte-swap on an array as needed when translating between 
     * big- and little-endian machines. 
     *
     * \param x the array to be byte swapped
     * \param the number of X values
     */
    template <class X> void byte_swap(X* x, int nx){
    
	switch(sizeof(X)){
      
	    case 1:
		break;
      
	    case 2:
		for(int i=0;i<nx;i++,x++)
		    byte_swap2((char*)x);
		break;
      
	    case 4:
		for(int i=0;i<nx;i++,x++)
		    byte_swap4((char*)x);
		break;
      
	    case 8:
		for(int i=0;i<nx;i++,x++)
		    byte_swap8((char*)x);
		break;
      
	    default:
		throw Subs_Error("Subs::swap_bytes(X*, int): unrecognised number of bytes/X = " + Subs::str(sizeof(X)));
	}
    }

    //! Linear interpolation
    /** Workhorse to carry out linear interpolation. Given two (x,y) points,
     * this routine computes the x value equivalent to an arbitrary y value for a line that runs through
     * the two points. the two X values must be different.
     * \param x1 X value of first point
     * \param y1 Y value of first point
     * \param x2 X value of second point
     * \param y2 Y value of second point
     * \param x  X value to try to compute equivalent Y
     * \return The y value equivalent to x
     */
    inline double linterp(double x1, double y1, double x2, double y2, double x){
	return (y1*(x2-x)+y2*(x-x1))/(x2-x1);
    }

    //! locate a value in an ordered list

    /** locate locates the position of x in an ordered array xx[n]. It returns
     * an index j such that x is between xx[j-1] and xx[j]. If j=0 or n 
     * then x is out of range. Use 'locate' if you think that x will jump
     * around a lot.
     *
     * \param xx ordered array of values (either increasing or decreasing)
     * \param n number of elements
     * \param x value to locate
     * \return Returned index such that x lies between xx[j-1] and xx[j]
     */

    template <class X> 
    unsigned long int locate(const X* xx, unsigned long int n, X x){

	if(x == xx[0]){
	    return 1;
	}else if(x == xx[n-1]){
	    return n-1;
	}

	unsigned long int ju, jm, jl;
	bool ascnd;  

	jl = 0;
	ju = n+1;
	ascnd = (xx[n-1] >= xx[0]);
	while(ju - jl > 1){
	    jm = (ju+jl) >> 1;
	    if((x >= xx[jm-1]) == ascnd)
		jl = jm;
	    else
		ju = jm;
	}

	// because of need to work with positive indices (owing to unsigned int)
	// jl and ju are 1 more than one might expect, so return jl as the answer
	return jl;
    }


    //! hunts for a value in an ordered list 

    /** 'hunt' locates the position of x in an ordered array xx[n]. It returns
     * value jhi such that x is between xx[jhi-1] and xx[jhi]. If jhi=0 or n 
     * then x is off one end or the other of the range. 'Use 'hunt' when x does not change much between
     * calls.
     * \param xx ordered array of values (either increasing or decreasing)
     * \param n number of elements
     * \param x value to locate
     * \param jhi Index such that x lies between xx[jhi-1] and xx[jhi], initialised
     * at approx expected position.
     */

    template <class X> 
    void hunt(const X* xx, unsigned long int n, X x, unsigned long int& jhi){

	unsigned long int jm, jl, ju, inc;
	bool ascnd;  

	ascnd = (xx[n-1] > xx[0]);

	// We will look for jl and ju that satisy xx[jl-1] < x < xx[ju-1]
	// 1 offset to avoid negative jl given that jl-1 = 0 at left end.
	// Because of this, set jl=jhi, not ju.
	jl = jhi; 

	if(jl < 1 || jl > n){
	    // If initial value out of range, reset to full range (i.e. like 'locate')
	    // with one off each end of valid range for jl and ju.
	    jl = 0;
	    ju = n+1;

	}else{

	    inc = 1; // initial increment
	    if((x >= xx[jl-1]) == ascnd){ // hunt up
		if(jl == n) return; // already at max
		ju = jl + 1;
		while((x >= xx[ju-1]) == ascnd){
		    jl = ju;
		    inc <<= 1; // double increment
		    ju = jl + inc;
		    if(ju > n){ // ok, run off end
			ju = n+1;
			break; 
		    }
		}

	    }else{ // Hunt down
		if(jl == 1){ // off end
		    jhi = 0;
		    return;
		} 
		ju = jl--;
		while((x < xx[jl-1]) == ascnd){
		    ju = jl;
		    inc <<= 1; // double increment
		    if(inc >= ju){ // fallen off lower end
			jl = 0;
			break;
		    }else{
			jl = ju - inc;
		    }
		}
	    }
	}

	// A couple of special cases can be dealt with quickly.
	if(x == xx[n-1]){
	    jhi = n-1;
	    return;
	}
	if(x == xx[0]){
	    jhi = 1;
	    return;
	}

	// Now have initial values of jl and ju, just like 'locate', except,
	// with any luck they are a bit tighter at the start. 
	while(ju-jl > 1){
	    jm = (ju+jl) >> 1;
	    if((x > xx[jm-1]) == ascnd)
		jl = jm;
	    else
		ju = jm;
	}

	// because of need to work with positive indices (owing to unsigned int)
	// jl and ju are 1 more than one might expect, so return jl as the answer
	jhi = jl;
    }

    //! Sorts an array

    /**
     * quicksort sorts an array into ascending order with the quicksort
     * algorithm. It works on any type for which the operators 
     * '=', '<=', '>' and '<' are defined. quicksort is fast on average 
     * but can vary substantially in speed according to the ordering of the array. 
     * \param arr array to sort. Returned sorted.
     * \param n number of elements
     */

    template <class X> 
    void quicksort(X arr[], long int n){
	const int NSTACK=50;
	const int M=7;
	long int ir=n-1, j, k, istack[NSTACK];
	long int jstack=-1,i,l=0;
	X a;

	for(;;){
	    if(ir-l < M){
		for(j=l+1;j<=ir;j++){
		    a=arr[j];
		    for(i=j-1;i>=l;i--){
			if(arr[i] <= a) break;
			arr[i+1] = arr[i];
		    }
		    arr[i+1] = a;
		}
		if(jstack == -1) break;
		ir = istack[jstack--];
		l  = istack[jstack--];
	    }else{
		k = ((unsigned long)l+ir) >> 1;
		std::swap(arr[k],arr[l+1]);
		if(arr[l] > arr[ir]){
		    std::swap(arr[l],arr[ir]);
		}
		if(arr[l+1] > arr[ir]){
		    std::swap(arr[l+1],arr[ir]);
		}
		if(arr[l] > arr[l+1]){
		    std::swap(arr[l],arr[l+1]);
		}
		i=l+1;
		j=ir;
		a=arr[l+1];
		for(;;){
		    do i++; while (arr[i] < a);
		    do j--; while (arr[j] > a);
		    if(j < i) break;
		    std::swap(arr[i],arr[j]);
		}
		arr[l+1] = arr[j];
		arr[j]   = a;
		jstack  += 2;
		if(jstack >= NSTACK){
		    std::cerr << "NSTACK too small in qsort" << std::endl;
		    exit(EXIT_FAILURE);
		}
		if(ir-i+1 >= j-l){
		    istack[jstack]   = ir;
		    istack[jstack-1] = i;
		    ir=j-1;
		}else{
		    istack[jstack]   = j-1;
		    istack[jstack-1] = l;
		    l=i;
		}
	    }
	}
    }

    //! Returns a sorted index array

    /** heaprank' returns an index 'key' of an array 'arr' such that 
     * arr[key[i]] for i=0,1,2 etc ascends. It can work on any type
     * for which the '<' (less than) operator is defined.
     * \param arr array of values
     * \param key indexing array
     * \param n   number of elements
     */

    template <class X> 
    void heaprank(X arr[], unsigned long int key[], unsigned long int n){
	X x;
	unsigned long int i, l, ir, ik, j;
	for(i = 0; i < n; i++){
	    key[i] = i;
	}
	if(n == 1) return;

	l  = (n >> 1);
	ir = n-1;
	for(;;){
	    if(l > 0){
		ik = key[--l];
		x  = arr[ik];
	    }else{
		ik = key[ir];
		x  = arr[ik];
		key[ir] = key[0];
		if(--ir == 0){
		    key[0] = ik;
		    break;
		}
	    }
	    i = l;
	    j = l + l + 1;
	    while(j <= ir){
		if(j < ir && arr[key[j]] < arr[key[j+1]]) j++;
		if(x < arr[key[j]]){
		    key[i] = key[j];
		    i = j;
		    j = j + j + 1;
		}else break;
	    }
	    key[i] = ik;
	}
	return;
    }

    //! heapsort -- returns sorted array
    /** heapsort sorts an array into ascending order
     * \param arr array of values
     * \param n   number of elements
     */
    template <class X> 
    void heapsort(X arr[], unsigned long int n){
  
	X x;
	unsigned long int i, l, ir, j;

	if(n < 2) return;

	l  = (n >> 1);
	ir = n-1;
	for(;;){
	    if(l > 0){
		x  = arr[--l];
	    }else{
		x  = arr[ir];
		arr[ir] = arr[0];
		if(--ir == 0){
		    arr[0] = x;
		    break;
		}
	    }
	    i = l;
	    j = l + l + 1;
	    while(j <= ir){
		if(j < ir && arr[j] < arr[j+1]) j++;
		if(x < arr[j]){
		    arr[i] = arr[j];
		    i = j;
		    j = j + j + 1;
		}else break;
	    }
	    arr[i] = x;
	}
	return;
    }

    template <class X> 
    X qtrap(X (*func)(X x), X a, X b, X eps, int nmax){
	X trapzd(X (*func)(X x), X a, X b, int n);
	int j;
	X s, olds;
  
	olds = -1.e30;
	for(j=0; j < nmax; j++){
	    s = trapzd(func,a,b,j);
	    if(fabs(s-olds) < eps*fabs(olds)) return s;
	    olds = s;
	}
	std::cerr << "Too many steps in qtrap\n";
	return 0.0;
    }

    template <class X> 
    X trapzd(X (*func)(X x), X a, X b, unsigned int n){
	X x, tnm, sum, del;
	static X s;
	unsigned long int it, j;

	if(n == 0){

	    // start with end values

	    return (s=0.5*(b-a)*((*func)(a)+(*func)(b)));

	}else{

	    // compute it=2**(n-1). This will be the
	    // the number of function evaluations needed.

	    for(it=1, j=0; j<n-1; j++) it <<= 1;
	    tnm = it;
	    del = (b-a)/tnm;
	    x   = a+0.5*del;
	    for(sum=0., j=0; j<it; j++, x+=del) sum += (*func)(x);
	    s = 0.5*(s+(b-a)*sum/tnm);
	    return s;
	}
    }

    template <class Func, class X> 
    X trapzd(const Func& f, X a, X b, unsigned int n){
	X x, tnm, sum, del;
	static X s;
	unsigned long int it, j;

	if(n == 0){

	    // start with end values

	    return (s=0.5*(b-a)*(f(a)+f(b)));

	}else{

	    // compute it=2**(n-1). This will be the
	    // the number of function evaluations needed.

	    for(it=1, j=0; j<n-1; j++) it <<= 1;
	    tnm = it;
	    del = (b-a)/tnm;
	    x   = a+0.5*del;
	    for(sum=0., j=0; j<it; j++, x+=del) sum += f(x);
	    s = 0.5*(s+(b-a)*sum/tnm);
	    return s;
	}
    }

    //! Simpson rule integrator

    /**
     * qsimp carries out trapezoidal integration of a function
     * until a specified level of accuracy has been reached. It does
     * this by repeatedly calling trapzd until
     * the integral changes by less than a specified fraction of its
     * previous value, until a maximum counter has been reached or
     * if the integral stays at zero twice in a row.
     *
     * qsimp also takes an argument to force it to take a minimum
     * number of steps to avoid it stopping too early.
     * \param func  a 1D function 
     * \param a     lower limit to integrate from
     * \param b     upper limit to integrate to
     * \param eps   accuracy
     * \param nmin  mininum number to reach when calling trapzd
     * \param nmax  maximum sub-division factor.
     * \param print print diagnostic info 
     * \return Returns the integral.
     */

    template <class X> 
    X qsimp(X (*func)(X x), X a, X b, X eps, int nmin, int nmax,  bool print){
	int n;
	X s, st, ost, os;
  
	ost = os = -1.e30;
	for( n=0; n<nmax; n++){
	    st =trapzd(func,a,b,n);
	    s  = (4.0*st-ost)/3.0;
	    if(n > nmin)
		if(print) std::cerr << n << " " << s << " " << os << std::endl;
	    if(fabs(s-os) < eps*fabs(os) || 
	       (s == 0.0 && os == 0.0)) return s;
	    os  = s;
	    ost = st;
	}
	throw Subs_Error("qsimp<X>(X (*func)(X), X, X, X, int, int,  bool): too many steps");
    }

    //! Polynomial interpolation routine
    /** polint carries out Lagrangian interpolation of a series of x,y
     * values and returns the interpolated value and an estimate of
     * the error on it. For limited numbers of point interpolation,
     * call this with a short section of the array.
     *
     * \param xa  array of x values
     * \param ya  array of y values
     * \param n   number of values
     * \param x   value to extrapolate/interpolate to
     * \param y   extrapolated/interpolated value
     * \param dy estimated error on y
     */

    template <class X> 
    void polint(X xa[], X ya[], int n, X x, X &y, X &dy){
	int i, m, ns=0;
	X den, dif, dift, ho, hp, w;
	X *c = new X[n];
	X *d = new X[n];

	dif = fabs(x-xa[0]);
	for(i=0;i<n;i++){
	    if((dift=fabs(x-xa[i])) < dif){
		ns  = i;
		dif = dift;
	    }
	    c[i] = ya[i];
	    d[i] = ya[i];
	}
	y = ya[ns--];
	for(m=1; m<n; m++){
	    for(i=0; i<n-m; i++){
		ho = xa[i]   - x;
		hp = xa[i+m] - x;
		w  = c[i+1]  - d[i];
		if( (den=ho-hp) == 0.0) 
		    throw Subs_Error("Error in polint; two identical x values");
		den  = w/den;
		d[i] = hp*den;
		c[i] = ho*den;
	    }
	    y += (dy = (2*(ns+1) < (n-m) ? c[ns+1] : d[ns--]));
	}
	delete[] c;
	delete[] d;
	return;
    }

    //! Abstract class for basic function object
    /** This class is the base class for usage by any routine reaquiring a basic function
     * object representing a 1D function. It declares one function that is needed which 
     * returns the function value given the position. This is such a standard usage that it is 
     * called Sfunc
     */
    class Sfunc {

    public:

	//! The function call
	/** This routine should return the function value at position x.
	 * \param x  the position to evaluate the function and derivative
	 * \return the function value
	 */
	virtual double operator()(double x) = 0;
	virtual ~Sfunc(){}
    };

    //! Romberg integration

    /**
     * qromb carries out trapezoidal integration of a function
     * until a specified level of accuracy has been reached. It does
     * this by repeatedly calling trapzd until the integral changes 
     * by less than a specified fraction of its
     * previous value, until a maximum counter has been reached or
     * if the integral stays at zero twice in a row. It does this
     * by extrapolating the value and seeing how the extrapolation
     * changes with successively smaller steps. It is good for
     * relatively smooth function.
     *
     * qromb also takes an argument to force it to take a minimum
     * number of steps to avoid it stopping too early.
     * \param func  a 1D function 
     * \param a     lower limit to integrate from
     * \param b     upper limit to integrate to
     * \param eps   accuracy
     * \param nmin  mininum number to reach when calling trapzd
     * \param nmax  maximum sub-division factor.
     * \param print print diagnostic info 
     * \return Returns the integral.
     */

    template <class X> 
    X qromb(X (*func)(X x), X a, X b,  X eps, int nmin, 
	    int nmax, bool print){
	X ss, dss;
	const int NMAX=50;
	X s[NMAX], h[NMAX+1];
	int n;
    
	if(nmax > NMAX) throw Subs_Error("nmax > NMAX inside qromb");
	if(nmin >= nmax) throw Subs_Error("nmin >= nmax inside qromb");
    
	h[0] = 1.0;
	for(n=0; n<nmax; n++){
	    s[n] = trapzd(func, a, b, n);
	    if(n >= nmin){
		polint(h+n-nmin,s+n-nmin,nmin,0.0,ss,dss);
		if(print) std::cerr << n << " " << s[n] << " " << ss 
				    << " " << dss << std::endl;
		if(fabs(dss) <= eps*fabs(ss)) return ss;
	    }else if(print){
		std::cerr << n << " " << s[n] << std::endl;
	    }
	    h[n+1] = 0.25*h[n];
	}
	throw Subs_Error("Too many steps in qromb");
    }

    //! Romberg integration
    template <class Func, class X>
    X qromb(const Func& f, X a, X b,  X eps, int nmin, int nmax, bool print){
	X ss, dss;
	const int NMAX=50;
	X s[NMAX], h[NMAX+1];
	int n;

	if(nmax > NMAX) throw "nmax > NMAX inside qromb";
	if(nmin >= nmax) throw "nmin >= nmax inside qromb";

	h[0] = 1.0;
	for(n=0; n<nmax; n++){
	    s[n] = trapzd(f, a, b, n);
	    if(n >= nmin){
		polint(h+n-nmin,s+n-nmin,nmin,0.0,ss,dss);
		if(print) std::cerr << n << " " << s[n] << " " << ss 
				    << " " << dss << std::endl;
		if(fabs(dss) <= eps*fabs(ss)) return ss;
	    }else if(print){
		std::cerr << n << " " << s[n] << std::endl;
	    }
	    h[n+1] = 0.25*h[n];
	}
	std::cerr << "Too many steps in qromb\n";
	return ss;
    }

    //! 1D minimisation routine without derivatives
    double brent(double xstart, double x1, double x2, Sfunc& f, double tol, double& xmin);

    //! 1D minimisation with derivatives
    double dbrent(double ax, double bx, double cx, Sfunc& func, Sfunc& dfunc, double acc,
		  bool stopfast, double pref, double& xmin);

    //! Abstract class for rtsafe function
    /** This class is the base class for usage by the 'rtsafe'. It declares the one 
     * function that is needed.
     */
    class RTfunc {

    public:

	//! The function call
	/** This routine should compute the function value and derivative at position x.
	 * \param x  the poisition to evaluate the function and derivative
	 * \param f  the function value (returned)
	 * \param fd the derivative (returned)
	 */
	virtual void operator()(double x, double& f, double& fd) const = 0;
	virtual ~RTfunc(){}

    };

    //! Find root of a function
    double rtsafe(const RTfunc& func, double x1, double x2, double xacc);

    //! Selects the k-th smallest value in an array arr[n].
    /**
     * 'select' finds the k-th smallest element of an array. Although
     * this can be done by sorting, it turns out that there is a speed
     * advantage in just going for the one element of interest. It scales
     * linearly with the number of elements plus overheads.
     * \param arr array to select from. \b NB It is returned in a scrambled order!
     * \param n number of elements
     * \param k the element to choose, starting with k=0 as the smallest,
     * and running up to n-1 at the largest (if entered outside this range it will be truncated at these limits).
     */

    template <class T>
    T select(T* arr, int n, int k){
	k = std::max(0, std::min(n-1, k));
	int i, ir, j, l, mid;
	T a;
    
	l  = 0;
	ir = n-1;
	for(;;){
	    if(ir <= l+1){
		if(ir == l+1 && arr[ir] < arr[l]) std::swap(arr[l],arr[ir]);
		return arr[k];
	    }else{
		mid = (l+ir) >> 1;
		std::swap(arr[mid],arr[l+1]);
		if(arr[l] > arr[ir]) std::swap(arr[l],arr[ir]);
		if(arr[l+1] > arr[ir]) std::swap(arr[l+1],arr[ir]);
		if(arr[l] > arr[l+1]) std::swap(arr[l],arr[l+1]);
		i = l+1;
		j = ir;
		a = arr[l+1];
		for(;;){
		    do i++; while(arr[i] < a);
		    do j--; while(arr[j] > a);
		    if(j < i) break;
		    std::swap(arr[i],arr[j]);
		}
		arr[l+1] = arr[j];
		arr[j]   = a;
		if(j >= k) ir = j-1;
		if(j <= k) l  = i;
	    }
	}
    }

    //! Plot colours (PGPLOT)
    /** These are written in the same order as they appear in PGPLOT
     * since they are integer codes in PGPLOT. They make it more obvious
     * what colour is being plotted
     */
    enum PLOT_COLOUR {
	NONE=-1, BLACK, WHITE, RED, GREEN, BLUE, CYAN, PURPLE, YELLOW, 
	ORANGE, LIGHT_GREEN, DARK_GREEN, LIGHT_BLUE, DARK_BLUE, 
	PINK, DARK_GREY, LIGHT_GREY
    };

    //! Converts a string to equivalent colour flag
    PLOT_COLOUR what_colour(const std::string& colour);

    //! Planck function Bnu.
    double planck(double wave, double temp);

    //! Logarithmic derivative of Planck function Bnu wrt wavelength
    double dplanck(double wave, double temp);

    //! Logarithmic derivative of Planck function Bnu wrt T
    double dlpdlt(double wave, double temp);

    //! Abstract class for Runge-Kutta integration
    /** This class is the base class for usage by the Runge-Kutta integration routines
     *  which need a function which computes derivatives
     */
    class RKfunc {
    public:
	//! The function call
	/** This routine should compute the time derivative of a vector of y values given the time and y values
	 * \param t the time
	 * \param y the y values
	 * \param dydt the derivatives (returned)
	 */
	virtual void operator()(double t, double y[], double dydt[]) = 0;
	virtual ~RKfunc(){}
    };

    //! Fifth order Runge-Kutta integrator
    void rkqs(double &x, double y[], double dydx[], double yscale[], int n, 
	      double htry, double eps, double& hdid, double &hnext, RKfunc& derivs);


    //! Writes out a vector in binary format
    template <class X>
    void write(const std::vector<X>& vec, std::ostream& ostr){
	Subs::UINT4 nvec = vec.size();
	ostr.write((char*)&nvec, sizeof(Subs::UINT4));
	X v;
	for(Subs::UINT4 i=0; i<vec.size(); i++)
	    v = vec[i];
	ostr.write((char*)&v, sizeof(X));    
    }

    //! Read a vector in binary format
    template <class X>
    void read(std::vector<X>& vec, std::istream& istr){
	Subs::UINT4 nvec;
	istr.read((char*)&nvec, sizeof(nvec));
	vec.resize(nvec);
	X v;
	for(Subs::UINT4 i=0; i<vec.size(); i++){
	    istr.read((char*)&v, sizeof(X));
	    vec[i] = v;
	}
    }

    //! Median filter

    /* Median filter routine.
     * \param data the data to be filtered
     * \param filt the output filtered data
     * \param width the width in bins of the filter (odd)
     */
    template <class T>
    void medfilt(const Buffer1D<T>& data, Buffer1D<T>& filt, int width){
	
	if(width == 1){
	    filt = data;
	    return;
	}
	
	filt.resize(data.size());
	
	if(width % 2 == 0)
	    throw Subs_Error("void Subs::medfilt(const Buffer1D<T>&, Buffer1D<T>&, int): filter width must be odd");
	
	// Get work space
	unsigned long int *iptr = new unsigned long int[width];
	T                 *tptr = new T[width];
	
	int np     = data.size();
	int jstart = 0;
	int jstop  = std::min(width/2, np-1);
	unsigned long int nact = jstop+1;
	
	// Initialise, extending point at end which should not affect median.
	for(size_t j=0; j<nact; j++)
	    tptr[jstop+j] = data[j];
	
	Subs::heaprank(tptr+jstop, iptr, nact);
	
	// Sort
	for(size_t j=0; j<nact; j++)
	    tptr[j] = tptr[jstop+iptr[j]];
	
	filt[0] = tptr[nact/2];
	for(int i=1; i<np; i++){
      int jnew1 = std::max(i-width/2, 0);    
      int jnew2 = std::min(i+width/2, np-1);
	    
      // Update key index: if start has changed, delete first point
      if(jnew1 == jstart + 1){
		nact--;
		for(size_t j=0, jadd=0; j<nact; j++){
          if(iptr[j] == 0) jadd = 1;
          size_t jn = j + jadd;
          iptr[j] = iptr[jn] - 1;
          tptr[j] = tptr[jn];
		}
      }
      
      if(jnew2 == jstop + 1){
		T test = data[jnew2];
		int jtest   = locate(tptr, nact, test);
		if(jtest == int(nact)){
          iptr[nact] = nact;
          tptr[nact] = test;
		}else{
          for(int j=nact; j>jtest; j--){
			iptr[j] = iptr[j-1];
			tptr[j] = tptr[j-1];
          }
          tptr[jtest] = test;
          iptr[jtest] = nact;
		}
		nact++;
      }
      filt[i] = tptr[nact/2];
      jstart = jnew1;
      jstop  = jnew2;
	}
	delete[] iptr;
	delete[] tptr;
    }

    //! Returns a double from a string
    double string_to_double(const std::string& entry);

    //! Returns an int from a string
    int string_to_int(const std::string& entry);

    //! Returns a char from a string
    char string_to_char(const std::string& entry);

    //! Returns a bool from a string
    bool string_to_bool(const std::string& entry);

    //! Abstract class for powell and amoeba
    /** This class is the base class for usage by the simplex minimisation routine amoeba and
     * similar routines that need to know a function value given an Array1D of parameters. 
     */
    class Afunc {
    public:
	//! The function call
	/** This routine should return the function value after it has been passed 
	 * a vector of parameters. It is not defined as 'const' since you may well want to 
	 * alter the object during the computation.
	 */
	virtual double operator()(const Array1D<double>& vec) = 0;
	virtual ~Afunc(){}
    };

    //! Wrapper class to make 1D function out of vector one
    /** This converts a multi-D function into a 1D function given a starting point and direction
     * to move along from it. A 1D value x is converted to a position using r = p + x*d where
     * r, p and d are vectors.
     */
    class Safunc : public Sfunc {
    public:

	//! Constructor
	/** \param func mult-D function
	 *  \param p initial point
	 * \param  d direction to travel from p, normalised
	 * \param  scale the normalisation factors used for each dimension of d
	 */
	Safunc(Afunc& func, const Array1D<double>& p, const Array1D<double>& d, const Array1D<double>& scale);

	//! The function call
	double operator()(double x);

    private:

	Afunc& func;
	const Array1D<double>& p;
	const Array1D<double>& d;
	const Array1D<double>& scale;

    };

    //! Namespace for genetic algorithm routines
    namespace Genetic {

	//! Translates vector of parameters to a string
	std::string model_to_string(const Subs::Array1D<double>& vals);

	//! Translates a string to a vector of parameters
	void string_to_model(const std::string& model, Subs::Array1D<double>& vals);

	//! Mutates a model represented by a string
	void mutate(std::string& model, INT4& seed, double rate);
  
	//! Crosses two models to make another
	std::string cross(const std::string& model1, const std::string& model2, INT4& seed);

    };

    //! Routine for bracketting a minimum
    void mnbrak(double& ax, double& bx, double& cx, double& fa, double& fb, double& fc, Sfunc& func);

    //! Simplex minimisation routine
    void amoeba(std::vector<std::pair<Array1D<double>, double> >& params, double ftol, int nmax, Afunc& func, int& nfunc);

    //! Powell's 'direction set' method
    void powell(Array1D<double>& p, Buffer1D<Array1D<double> >& xi, const Array1D<double>& scale, double ftol, int itmax, int& iter, double& fret, Afunc& func);

    // Linear least-squares

    //! Abstract class that for llsqr function object routines
    class Llfunc {
    public:

	//! Evaluate the nfunc function values at x, returning them in v
	/** \param x the ordinate to evaluate the functions at
	 * \param v the values
	 */
	virtual void eval(double x, double* v) const = 0;

	//! Returns the number of functions
	virtual int get_nfunc() const = 0;

	//! Destructor
	virtual ~Llfunc() {}
    };

    //! General linear least square fitter
    void llsqr(int ndata, const double* y, const float* e, int nfunc, double** func, double* coeff, double** covar);

    //! Evaluates Chi**2 after application of llsqr
    double llsqr_chisq(int ndata, const double* y, const float* e, int nfunc, double** func, const double* coeff);

    //! Evaluates reduced Chi**2 after application of llsqr
    double llsqr_reduced_chisq(int ndata, const double* y, const float* e, int nfunc, double** func, const double* coeff);

    //! Evaluates fitted values after application of llsqr
    void llsqr_eval(int ndata, int nfunc, double** func, const double* coeff, double* fit);

    //! Carries out sigma clipping after llsqr
    int llsqr_reject(int ndata, const double* y, float* e, int nfunc, double** func, const double* coeff, double thresh, bool slow);

    //! General linear least square fitter
    void llsqr(int ndata, const double* x, const double* y, const float* e, const Llfunc& func, double* coeff, double** covar);

    //! Evaluates Chi**2 after application of llsqr
    double llsqr_chisq(int ndata, const double* x, const double* y, const float* e, const Llfunc& func, const double* coeff);

    //! Evaluates reduced Chi**2 after application of llsqr
    double llsqr_reduced_chisq(int ndata, const double* x, const double* y, const float* e, const Llfunc& func, const double* coeff);

    //! Evaluates fitted values after application of llsqr
    void llsqr_eval(int ndata, const double* x, const Llfunc& func, const double* coeff, double* fit);

    //! Carries out sigma clipping after llsqr
    int llsqr_reject(int ndata, const double* x, const double* y, float* e, const Llfunc& func, const double* coeff, double thresh, bool slow);

}; // end of subs namespace

// 1D buffer template class
namespace Subs {

    //! Error class
    class Buffer1D_Error : public Subs_Error {
    public:
	Buffer1D_Error() : Subs_Error("") {};
	Buffer1D_Error(const std::string& str) : Subs_Error(str) {};
    };

    //! ASCII output
    template <class X>
    std::ostream& operator<<(std::ostream& s, const Buffer1D<X>& vec);

    //! ASCII input
    template <class X>
    std::istream& operator>>(std::istream& s, Buffer1D<X>& vec);

    //! Buffer class for handling memory.
    /** A class designed to supply a safe 1D array, 'safe' in
     * that it is deallocated whenever the object goes out of scope.
     * It creates a pointer to an array which can then be used
     * in the usual way as an array. The pointer is automatically
     * deleted when its Buffer1D goes out of scope. It also stores the number
     * of pixels and the number of memory elements allocated. The latter can
     * be larger than the number of pixels to make extension of the array size
     * more efficient. Buffer1D is designed to handle any type of data and therefore
     * does not supply operations such as addition; look at Array1D for such 
     * specialisation which inherits Buffer1D and adds such facilities. Buffer1D can
     * therefore contain complex objects as its elements. Such objects will need to
     * support the operations of assignment and ASCII I/O. Binary I/O functions
     * write, read and skip are provided but should only be used on objects with
     * no pointers because they will store the pointers but not whatever they
     * point to.
     */
    template <class X>
    class Buffer1D {
    public:

	//! Default constructor, makes NULL pointer.
	Buffer1D() : buff(NULL), npix(0), nmem(0) {}

	//! Constructor grabbing space for npix points
	Buffer1D(int npix);

	//! Constructor grabbing space for npix points but with nmem elements allocated
	Buffer1D(int npix, int nmem);

	//! Constructor from a file name
	Buffer1D(const std::string& file);

	//! Copy constructor
	Buffer1D(const Buffer1D& obj);

	//! Constructor from a vector
	Buffer1D(const std::vector<X>& obj);

	//! Destructor
	virtual ~Buffer1D(){
	    if(buff != NULL)
		delete[] buff;
	}

	//! Assignment
	Buffer1D& operator=(const Buffer1D& obj);

	//! Assignment
	Buffer1D& operator=(const X& con);

	//! Returns the number of pixels
	int get_npix() const {return npix;}

	//! Returns the number of pixels
	int size() const {return npix;}

	//! Returns the number of pixels allocated in memory
	int mem() const {return nmem;}

	//! Change number of pixels
	void resize(int npix);

	//! Zeroes number of elements in the array (but does not change memory allocation)
	void clear(){npix = 0;}

	//! Change number of pixels directly, no other effect (experts only)
	void set_npix(int npix){
	    if(npix > nmem)
		throw Buffer1D_Error("Subs::Buffer1D::set_npix(int): attempt to set npix > nmem not allowed.");
	    this->npix = npix;
	}

	//! Element access
	const X& operator[](int i) const {
	    return buff[i];
	}

	//! Element access
	X& operator[](int i) {
	    return buff[i];
	}

	//! Add another value to the end of the buffer, preserving all other data 
	void push_back(const X& value);

	//! Remove a pixel
	void remove(int index);

	//! Insert a pixel
	void insert(int index, const X& value);

	//! Conversion operator for interfacing with normal C style arrays.
	operator X*(){return buff;}

	//! Returns pointer to the buffer to interface with functions requiring normal C style arrays.
	X* ptr() {return buff;}

	//! Returns pointer to the buffer to interface with functions requiring normal C style arrays.
	const X* ptr() const {return buff;}

	//! Read from an ASCII file
	void load_ascii(const std::string& file);

	//! Write out a Buffer1D to a binary file
	void write(std::ostream& s) const;

	//! Skip a Buffer1D in a binary file
	static void skip(std::istream& s, bool swap_bytes);

	//! Read a poly from a binary file
	void read(std::istream& s, bool swap_bytes);

	friend std::ostream& operator<<<>(std::ostream& s, const Buffer1D& vec);

	friend std::istream& operator>><>(std::istream& s, Buffer1D& vec);

    protected:

	//! The pointer; used extensively in Array1D, change at your peril!
	X* buff;

	//! For derived class ASCII input
	virtual void ascii_input(std::istream& s);

	//! For derived class ASCII output
	virtual void ascii_output(std::ostream& s) const;

    private:
    
	// number of pixels and number of memory elements
	int npix;
	int nmem;

    };


    /** This constructor gets space for exactly npix points
     * \param npix the number of points to allocate space for
     */
    template <class X> 
    Buffer1D<X>::Buffer1D(int npix) : npix(npix), nmem(npix) {
	if(npix < 0)
	    throw Buffer1D_Error("Subs::Buffer1D<>(int): attempt to allocate < 0 point");
	if(npix == 0){
	    buff = NULL;
	}else{
	    if((buff = new(std::nothrow) X [nmem]) == NULL){
		this->npix = nmem = 0;
		throw Buffer1D_Error("Subs::Buffer1D::Buffer1D(int): failed to allocate memory");
	    }
	}
    }

    /** This constructor gets space for exactly npix points
     * \param npix the number of pixels
     * \param nmem the number of memory elements
     */
    template <class X> 
    Buffer1D<X>::Buffer1D(int npix, int nmem) : npix(npix), nmem(nmem) {
	if(npix < 0)
	    throw Buffer1D_Error("Subs::Buffer1D<>(int, int): attempt to set < 0 pixels");
	if(nmem < npix)
	    throw Buffer1D_Error("Subs::Buffer1D<>(int, int): must allocate at least as many memory elements as pixels");
	if(nmem == 0){
	    buff = NULL;
	}else{
	    if((buff = new(std::nothrow) X [nmem]) == NULL){
		this->npix = nmem = 0;
		throw Buffer1D_Error("Subs::Buffer1D<>(int, int): failure to allocate " + Subs::str(nmem) + " points.");
	    }
	}
    }

    /** Constructor by reading a file
     * \param file the file to read, an ASCII file.
     */
    template <class X>
    Buffer1D<X>::Buffer1D(const std::string& file) : buff(NULL), npix(0), nmem(0) {
	try{
	    load_ascii(file);
	}
	catch(const Buffer1D_Error& err){
	    throw Buffer1D_Error("Buffer1D<X>::Buffer1D(const std::string&): error constructing from a file " + err);
	}
    }

    /** Copy constructor to make an element by element copy of an object
     */
    template <class X> 
    Buffer1D<X>::Buffer1D(const Buffer1D<X>& obj) : npix(obj.npix), nmem(obj.npix) {
	if(nmem == 0){
	    buff = NULL;
	}else{
	    if((buff = new(std::nothrow) X [nmem]) == NULL){
		npix = nmem = 0;
		throw Buffer1D_Error("Subs::Buffer1D<>(const Buffer1D<>&): failure to allocate " + Subs::str(nmem) + " points.");
	    }
	    for(int i=0; i<npix; i++)
		buff[i] = obj.buff[i];
	}
    }

    /** Constructor to make an element by element copy of a vector
     */
    template <class X> 
    Buffer1D<X>::Buffer1D(const std::vector<X>& obj) : npix(obj.size()), nmem(obj.size()) {
	if(nmem == 0){
	    buff = NULL;
	}else{
	    if((buff = new(std::nothrow) X [nmem]) == NULL){
		npix = nmem = 0;
		throw Buffer1D_Error("Subs::Buffer1D<>(const std::vector<>&): failure to allocate " + Subs::str(nmem) + " points.");
	    }
	    for(int i=0; i<npix; i++)
		buff[i] = obj[i];
	}
    }

    /** Sets one Buffer1D equal to another.
     */
    template <class X> 
    Buffer1D<X>& Buffer1D<X>::operator=(const Buffer1D<X>& obj){

	if(this == &obj) return *this;

	// First check whether we can avoid reallocation of memory
	if(buff != NULL){
	    if(obj.npix <= nmem){
		npix = obj.npix;
		for(int i=0; i<npix; i++)
		    buff[i] = obj.buff[i];
		return *this;
	    }else{
		delete[] buff;
	    }
	}
	
	// Allocate memory
	npix = nmem = obj.npix;
	if(nmem == 0){
	    buff = NULL;
	}else{
	    if((buff = new(std::nothrow) X [nmem]) == NULL){
		npix = nmem = 0;
		throw Buffer1D_Error("Subs::Buffer1D<>(const Buffer1D<>&): failure to allocate " + Subs::str(nmem) + " points.");
	    }
	}

	// Finally copy
	for(int i=0; i<npix; i++)
	    buff[i] = obj.buff[i];

	return *this;
    }

    /** Sets a Buffer1D to a constant
     */
    template <class X> 
    Buffer1D<X>& Buffer1D<X>::operator=(const X& con){

	for(int i=0; i<npix; i++)
	    buff[i] = con;

	return *this;
    }

    /** This changes the number of pixels. It does not preserve the data in general.
     * \param npix the new array size
     */
    template <class X>
    void Buffer1D<X>::resize(int npix){
	if(buff != NULL){
	    if(npix <= nmem){
		this->npix = npix;
		return;
	    }else{
		delete[] buff;
	    }
	}
    
	this->npix = nmem = npix;
	if(nmem < 0)
	    throw Buffer1D_Error("Subs::Buffer1D::resize(int): attempt to allocate < 0 points");
	if(nmem == 0){
	    buff = NULL;
	    return;
	}
	if((buff = new(std::nothrow) X [nmem]) == NULL){
	    this->npix = nmem = 0;
	    throw Buffer1D_Error("Subs::Buffer1D::resize(int): failed to allocate new memory");
	}
    }

    /** This routine adds a new value to the end of a buffer, increasing the 
     * memory allocated if need be.
     * \param value new value to add to the end
     */ 
    template <class X>
    void Buffer1D<X>::push_back(const X& value){
	if(npix < nmem){
	    buff[npix] = value;
	    npix++;
	}else{
	    nmem *= 2;
	    nmem  = (nmem == 0) ? 1 : nmem;  
	    X* temp;
	    if((temp = new(std::nothrow) X [nmem]) == NULL){
		nmem /= 2;
		throw Buffer1D_Error("Subs::Buffer1D::push_back(const X&): failed to extend memory");
	    }
	    for(int i=0; i<npix; i++)
		temp[i] = buff[i];
	    temp[npix] = value;
	    npix++;
	    delete[] buff;
	    buff = temp;
	}
    }

    /** This routine removes a pixel at a given index. The memory allocated
     * is not changed.
     * \param index the pixel to be removed
     */ 
    template <class X>
    void Buffer1D<X>::remove(int index){
	npix--;
	for(int i=index; i<npix; i++)
	    buff[i] = buff[i+1];
    } 

    /** This routine inserts a pixel at a given index. The memory allocated
     * may have to increase
     * \param index the pixel to be removed
     */ 
    template <class X>
    void Buffer1D<X>::insert(int index, const X& value){

	if(npix < nmem){
	    for(int i=npix; i>index; i--)
		buff[i] = buff[i-1];
	    buff[index] = value;
	    npix++;
	}else{
	    nmem *= 2;
	    nmem  = (nmem == 0) ? 1 : nmem;  
	    X* temp;
	    if((temp = new(std::nothrow) X [nmem]) == NULL){
		nmem /= 2;
		throw Buffer1D_Error("Subs::Buffer1D::insert(int, const X&): failed to extend memory");
	    }
	    for(int i=0; i<index; i++)
		temp[i] = buff[i];
	    for(int i=npix; i>index; i--)
		temp[i] = buff[i-1];
	    temp[index] = value;
	    npix++;
	    delete[] buff;
	    buff = temp;
	}
    } 

    //! Binary output
    template <class X>
    void Buffer1D<X>::write(std::ostream& s) const {
	s.write((char*)&npix, sizeof(int));
	s.write((char*)buff,  sizeof(X[npix]));
    }

    //! Binary input
    template <class X>
    void Buffer1D<X>::read(std::istream& s, bool swap_bytes) {
	s.read((char*)&npix, sizeof(int));
	if(!s) return;
	if(swap_bytes) npix = Subs::byte_swap(npix);
	this->resize(npix);
	s.read((char*)buff,  sizeof(X[npix]));
	if(swap_bytes) Subs::byte_swap(buff, npix);
    }

    //! Binary skip
    template <class X>
    void Buffer1D<X>::skip(std::istream& s, bool swap_bytes) {
	int npixel;
	s.read((char*)&npixel, sizeof(int));
	if(!s) return;
	if(swap_bytes) npixel = Subs::byte_swap(npixel);
	s.ignore(npixel*sizeof(X));
    }

    /* Loads data into a Buffer1D from an ASCII file with one
     * element per line. The elements must support ASCII input.
     * Define a suitable structure for complex input. Lines starting with
     * # are skipped.
     * \param file the file name to load
     */
    template <class Type> 
    void Buffer1D<Type>::load_ascii(const std::string& file){

	std::ifstream fin(file.c_str());  
	if(!fin)
	    throw Buffer1D_Error("void Buffer1D<>::load~_ascii(const std::string&): could not open " + file);

	// Clear the buffer
	this->resize(0);
	Type line;
	char c;
	while(fin){
	    c = fin.peek();
	    if(!fin) break;
	    if(c == '#' || c == '\n'){
		while(fin.get(c)) if(c == '\n') break;
	    }else{
		if(fin >> line) this->push_back(line);
		while(fin.get(c)) if(c == '\n') break; // ignore the rest of the line
	    }
	}
	fin.close();

    }

    template <class X>
    std::ostream& operator<<(std::ostream& s, const Buffer1D<X>& vec){
	vec.ascii_output(s);
	return s;
    }

    template <class X>
    void Buffer1D<X>::ascii_output(std::ostream& s) const {
	if(!s) return;
	s << this->size();
	for(int i=0; i<this->size(); i++)
	    s << " " << (*this)[i];
    }

    template <class X>
    std::istream& operator>>(std::istream& s, Buffer1D<X>& vec){
	vec.ascii_input(s);
	return s;
    }

    template <class X>
    void Buffer1D<X>::ascii_input(std::istream& s){
	if(!s) return;
	int nelem;
	s >> nelem;
	this->resize(nelem);
	for(int i=0; i<this->size(); i++)
	    s >> (*this)[i];
    }

}

// more subs namespace
namespace Subs {
  
    //! Returns a sorted index array

    /** heaprank' returns an index 'key' of an array 'arr' such that 
     * arr[key[i]] for i=0,1,2 etc ascends. It can work on any type
     * for which the '<' (less than) operator is defined.
     * \param arr array of values
     * \param key indexing array
     */

    template <class X> 
    void heaprank(const Buffer1D<X>& arr, Buffer1D<int>& key){
	X x;
	int i, l, ir, ik, j, n = arr.size();
	key.resize(n);
	for(i = 0; i < n; i++){
	    key[i] = i;
	}
	if(n == 1) return;

	l  = (n >> 1);
	ir = n-1;
	for(;;){
	    if(l > 0){
		ik = key[--l];
		x  = arr[ik];
	    }else{
		ik = key[ir];
		x  = arr[ik];
		key[ir] = key[0];
		if(--ir == 0){
		    key[0] = ik;
		    break;
		}
	    }
	    i = l;
	    j = l + l + 1;
	    while(j <= ir){
		if(j < ir && arr[key[j]] < arr[key[j+1]]) j++;
		if(x < arr[key[j]]){
		    key[i] = key[j];
		    i = j;
		    j = j + j + 1;
		}else break;
	    }
	    key[i] = ik;
	}
	return;
    }

}


#endif
