#include "trm/subs.h"

namespace Subs {
    void rkck(double &x, double y[], double dydx[], int n, double h, double yout[], double yerr[], RKfunc& derivs);
}

//! Fifth order Runge-Kutta integrator
/** Runge-Kutta integrator from Numerical Recipes, adapted
 * for function objects
 * \param x      the independent variable (input)
 * \param y      array of dependent variables (input)
 * \param dydx   array of derivatives of dependent variables (input)
 * \param yscale scaling factors for measuring the error
 * \param n      the number of variables
 * \param htry   the stepsize to attempt
 * \param eps    accuracy for each equation after allowing for the scaling factors
 * \param hdid   the step actually carried out
 * \param hnext  suggested next step size
 * \param derivs Function object with arguments (x,y,dydx) which computes the derivatives
 * given x and y
 */
void Subs::rkqs(double &x, double y[], double dydx[], double yscale[], int n, 
	  double htry, double eps, double& hdid, double &hnext, RKfunc& derivs){

    // Constants
    const double SAFETY  = 0.9;
    const double PGROW   = -0.2;
    const double PSHRINK = -0.25;
    const double ERRCON  = 1.89e-4;
  
    // variables
    int i;
    double errmax, h, htemp, xnew;
    double *yerr  = new double [n];
    double *ytemp = new double [n];
  
    // code
    h = htry;
    for(;;){
	rkck(x,y,dydx,n,h,ytemp,yerr,derivs);
	errmax = 0.;
	for(i=0;i<n;i++)
	    errmax = std::max(errmax, fabs(yerr[i]/yscale[i]));
	errmax /= eps;
	if(errmax <= 1.0) break;
	htemp = SAFETY*h*pow(errmax,PSHRINK);
    
	// Truncation error too large, reduce stepsize
	h = (h >= 0. ? std::max(htemp, 0.1*h) : std::min(htemp, 0.1*h));
	xnew = x + h;
	if(xnew == x) 
	    throw Subs_Error("void Subs::rkqs(double&, double, double, double, int, double, double, double&, double&, RKfunc&): stepsize underflow.");
    }
    if(errmax > ERRCON){
	hnext = SAFETY*h*pow(errmax, PGROW);
    }else{
	hnext = 5.0*h;
    }
    x += (hdid=h);
    for(i=0; i<n; i++)
	y[i] = ytemp[i];
    delete[] yerr;
    delete[] ytemp;
}    

namespace Subs {
    //! Takes a Cash-Karp Runge-Kutta step
    /** Runge-Kutta integrator from Numerical Recipes, adapted
     * for function objects
     * \param x      the independent variable (input)
     * \param y      array of dependent variables (input)
     * \param dydx   array of derivatives of dependent variables (input)
     * \param n      the number of variables
     * \param h      the stepsize
     * \param yout   incremented variables
     * \param yerr   estimate of truncation error from fourth-order method
     * \param derivs Function object with arguments (x,y,dydx) which computes the derivatives
     * given x and y
     */
    void rkck(double &x, double y[], double dydx[], int n, double h, double yout[], double yerr[], Subs::RKfunc& derivs){
	
	// Constants, lots of them
	const double a2=0.2, a3=0.3, a4=0.6, a5=1.0, a6=0.875;
	const double b21=0.2, b31=3./40., b32=9./40., b41=0.3, b42=-0.9, b43=1.2;
	const double b51=-11./54., b52=2.5, b53=-70./27., b54=35./27., b61=1631./55296.;
	const double b62=175./512., b63=575./13824., b64=44275./110592., b65=253./4096.;
	const double c1=37./378., c3=250./621., c4=125./594., c6=512./1771., dc5=-277./14336.;
	const double dc1=c1-2825./27648.,dc3=c3-18575./48384., dc4=c4-13525/55296., dc6=c6-0.25;
	
	int i;
	
	// Arrays
	double *ak2   = new double[n];
	double *ak3   = new double[n];
	double *ak4   = new double[n];
	double *ak5   = new double[n];
	double *ak6   = new double[n];
	double *ytemp = new double[n];
	
	// Take 6 steps
	for(i=0; i<n; i++)
	    ytemp[i] = y[i] + h*b21*dydx[i];
	derivs(x+a2*h, ytemp, ak2);
	for(i=0; i<n; i++)
	    ytemp[i] = y[i] + h*(b31*dydx[i]+b32*ak2[i]);
	derivs(x+a3*h, ytemp, ak3);
	for(i=0; i<n; i++)
	    ytemp[i] = y[i] + h*(b41*dydx[i]+b42*ak2[i]+b43*ak3[i]);
	derivs(x+a4*h, ytemp, ak4);
	for(i=0; i<n; i++)
	    ytemp[i] = y[i] + h*(b51*dydx[i]+b52*ak2[i]+b53*ak3[i]+b54*ak4[i]);
	derivs(x+a5*h, ytemp, ak5);
	for(i=0; i<n; i++)
	    ytemp[i] = y[i] + h*(b61*dydx[i]+b62*ak2[i]+b63*ak3[i]+b64*ak4[i]+b65*ak5[i]);
	derivs(x+a6*h, ytemp, ak6);
	
	// Accumulate increments
	for(i=0; i<n; i++)
	    yout[i] = y[i] + h*(c1*dydx[i]+c3*ak3[i]+c4*ak4[i]+c6*ak6[i]);
	
	// Estimate error from difference between fourth and fifth
	for(i=0; i<n; i++)
	    yerr[i] = h*(dc1*dydx[i]+dc3*ak3[i]+dc4*ak4[i]+dc5*ak5[i]+dc6*ak6[i]);
	
	delete[] ytemp;
	delete[] ak6;
	delete[] ak5;
	delete[] ak4;
	delete[] ak3;
	delete[] ak2;
    }
}
