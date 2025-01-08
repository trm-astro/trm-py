#include "trm/subs.h"
#include "trm/buffer2d.h"

/** llsqr carries out a general linear least squares
 * fit. It uses the normal equations approach and thus is senstive to degeneracy
 * but fast. This version used precompted values of the functions so is memory hungry
 * but fast.
 * \param ndata the number of data points
 * \param y    the data 
 * \param e    uncertainties on the data, <= 0 to mask a point
 * \param nfunc the number of functions
 * \param func  the functions as the columns of a matrix, i.e. func[2][4] is the value of the 3rd pixel of the 5th function. 
 * Dimensions: ny=ndata, nx=nfunc
 * \param coeff the fitted coefficients, i.e. the nfunc multipliers of the functions that lead to the
 * best (in a least squares sense) fit to the data (returned)
 * \param covar the covariance matrix (nfunc by nfunc), returned
 */

void Subs::llsqr(int ndata, const double* y, const float* e, int nfunc, double** func, double* coeff, double** covar){

  if(ndata < 1)
    throw Subs_Error("Subs::llsqr(int, const double*, const float*, int, const double**, double*) : < 1 data point");
  if(nfunc < 1)
    throw Subs_Error("Subs::llsqr(int, const double*, const float*, int, const double**, double*) : < 1 function");

  Buffer2D<double> beta(nfunc, 1);

  // Intialise
  for(int j=0; j<nfunc; j++){
    beta[j][0] = 0.;
    for(int k=0; k<nfunc; k++)
      covar[j][k] = 0.;
  }

  // Accumulate matrices
  for(int i=0; i<ndata; i++){
    if(e[i] > 0.){
      double weight = 1./Subs::sqr(e[i]);
      for(int j=0; j<nfunc; j++){
	double wt = weight*func[i][j];
	for(int k=0; k<=j; k++)
	  covar[j][k] += wt*func[i][k];
	beta[j][0] += wt*y[i];
      } 
    }
  }

  // fill rest of matric by symmetry
  for(int j=1; j<nfunc; j++)
    for(int k=0; k<j; k++)
      covar[k][j] = covar[j][k];

  // Solve normal equations
  gaussj(covar, nfunc, beta, 1);

  // return for coefficients
  for(int i=0; i<nfunc; i++)
    coeff[i] = beta[i][0];

}

/** llsqr_eval calculates the fit to the data after a run of llsqr
 * \param ndata the number of data points
 * \param nfunc the number of functions
 * \param func  the functions as the columns of a matrix, i.e. func[2][4] is the value of the 3rd pixel of the 5th function. 
 * Dimensions: ny=ndata, nx=nfunc
 * \param coeff the fitted coefficients, i.e. the nfunc multipliers of the functions that lead to the
 * best (in a least squares sense) fit to the data
 * \param fit the fitted values
 * \return the Chi**2
 */

void Subs::llsqr_eval(int ndata, int nfunc, double** func, const double* coeff, double* fit){

  if(ndata < 1)
    throw Subs_Error("Subs::llsqr_eval(int, int, const double**, double*) : < 1 data point");
  if(nfunc < 1)
    throw Subs_Error("Subs::llsqr_eval(int, int, const double**, double*) : < 1 function");
  
  for(int i=0; i<ndata; i++){
    double ymodel = 0.;
    for(int j=0; j<nfunc; j++)
      ymodel += coeff[j]*func[i][j];
    fit[i] = ymodel;
  }
}

/** llsqr_chisq evaluates the chi**2 value after a run of llsqr
 * \param ndata the number of data points
 * \param y    the data (dimension ndata)
 * \param e    uncertainties on the data, <= 0 to mask a point (dimension ndata)
 * \param nfunc the number of functions 
 * \param func  the functions as the columns of a matrix, i.e. func[2][4] is the value of the 3rd pixel of the 5th function. 
 * Dimensions: ny=ndata, nx=nfunc
 * \param coeff the fitted coefficients, i.e. the nfunc multipliers of the functions that lead to the
 * best (in a least squares sense) fit to the data
 * \return the Chi**2
 */

double Subs::llsqr_chisq(int ndata, const double* y, const float* e, int nfunc, double** func, const double* coeff){

  if(ndata < 1)
    throw Subs_Error("Subs::llsqr_chisq(int, const double*, const float*, int, const double**) : < 1 data point");
  if(nfunc < 1)
    throw Subs_Error("Subs::llsqr_chisq(int, const double*, const float*, int, const double**) : < 1 function");

  double chisq = 0.;
  
  for(int i=0; i<ndata; i++){
    if(e[i] > 0.){
      double ymodel = 0.;
      for(int j=0; j<nfunc; j++)
	ymodel += coeff[j]*func[i][j];
      chisq += Subs::sqr((y[i]-ymodel)/e[i]);
    }
  }
  return chisq;
}

/** llsqr_reduced_chisq evaluates the reduced chi**2 value after a run of llsqr
 * \param ndata the number of data points
 * \param y    the data (dimension ndata)
 * \param e    uncertainties on the data, <= 0 to mask a point (dimension ndata)
 * \param nfunc the number of functions 
 * \param func  the functions as the columns of a matrix, i.e. func[2][4] is the value of the 3rd pixel of the 5th function. 
 * Dimensions: ny=ndata, nx=nfunc
 * \param coeff the fitted coefficients, i.e. the nfunc multipliers of the functions that lead to the
 * best (in a least squares sense) fit to the data
 * \return the Chi**2
 */

double Subs::llsqr_reduced_chisq(int ndata, const double* y, const float* e, int nfunc, double** func, const double* coeff){

  if(ndata < 1)
    throw Subs_Error("Subs::llsqr_reduced_chisq(int, const double*, const float*, int, const double**) : < 1 data point");
  if(nfunc < 1)
    throw Subs_Error("Subs::llsqr_reduced_chisq(int, const double*, const float*, int, const double**) : < 1 function");

  double chisq = 0.;
  int ndof = - nfunc;
  for(int i=0; i<ndata; i++){
    if(e[i] > 0.){
      double ymodel = 0.;
      for(int j=0; j<nfunc; j++)
	ymodel += coeff[j]*func[i][j];
      chisq += Subs::sqr((y[i]-ymodel)/e[i]);
      ndof++;
    }
  }
  if(ndof < 1)
    throw Subs_Error("Subs::llsqr_reduced_chisq(int, const double*, const float*, int, const double**) : < 1 degree of freedom");
  return chisq/ndof;
}

/** llsqr_reject carries out a rejection cycle after a run of llsqr. Data points
 * are masked by making error bars negative.
 *
 * \param ndata the number of data points
 * \param y    the data 
 * \param e    uncertainties on the data, will be modified 
 * \param nfunc the number of functions
 * \param func  the functions as the columns of a matrix, i.e. func[2][4] is the value of the 3rd pixel of the 5th function. 
 * Dimensions: ny=ndata, nx=nfunc
 * \param coeff the fitted coefficients, i.e. the nfunc multipliers of the functions that lead to the
 * best (in a least squares sense) fit to the data
 * \param thresh the threshold in terms of sigma for rejection.
 * \param slow true to just reject the worst point, else all points above threshold will go
 * \return number of points rejected.
 */

int Subs::llsqr_reject(int ndata, const double* y, float* e, int nfunc, double** func, const double* coeff, double thresh,
		       bool slow){

  if(ndata < 1)
    throw Subs_Error("Subs::llsqr_reject(int, const double*, float*, int, const double**, double, bool) : < 1 data point");
  if(nfunc < 1)
    throw Subs_Error("Subs::llsqr_reject(int, const double*, float*, int, const double**, double, bool) : < 1 function");
  if(thresh <= 0.)
    throw Subs_Error("Subs::llsqr_reject(int, const double*, float*, int, const double**, double, bool) : thresh <= 0.");

  double reduced_chisq = llsqr_reduced_chisq(ndata, y, e, nfunc, func, coeff);

  int nrej = 0, iworst = -1;
  double worst = -1., limit = sqrt(reduced_chisq)*thresh;
  for(int i=0; i<ndata; i++){
    if(e[i] > 0.){
      double ymodel = 0.;
      for(int j=0; j<nfunc; j++)
	ymodel += coeff[j]*func[i][j];
      double dev = fabs(y[i]-ymodel)/e[i];
      if(dev > limit){
	if(slow){
	  if(dev > worst){
	    worst  = dev;
	    iworst = i;
	  }
	}else{
	  e[i] = - e[i];
	  nrej++;
	}
      }
    }
  }
  if(slow && iworst >= 0){
    e[iworst] = -e[iworst];
    nrej = 1;
  }

  return nrej;
}


//! General linear least square fitter

/** llsqr carries out a general linear least squares fit. It uses the normal equations 
 * approach and thus is sensitive to degeneracy but fast. This version uses a function object 
 * to evaluate the functions to be fitted to save memory at the expense of speed.
 * \param ndata the number of data points
 * \param x    the X values
 * \param y    the Y values
 * \param e    uncertainties on the Y values (<= 0 to mask)
 * \param func the functions in the form of a function object that has member functions of the form
 * 'void eval(double x, double* v)' to evaluate the nfunc function values at x, returning them in v,  and 
 * 'int get_nfunc()' const which returns the number of functions.
 * \param coeff the fitted coefficients, i.e. the nfunc multipliers of the functions that lead to the
 * best (in a least squares sense) fit to the data (returned)
 * \param covar the covariance matrix (nfunc by nfunc), returned
 */

void Subs::llsqr(int ndata, const double* x, const double* y, const float* e, const Llfunc& func, double* coeff, double** covar){
  
  if(ndata < 1)
    throw Subs_Error("Subs::llsqr(int, const double*, const double*, const float*, const Llfunc&, double*, double**) : < 1 data point");
  if(func.get_nfunc() < 1)
    throw Subs_Error("Subs::llsqr(int, const double*, const double*, const float*, const Llfunc&, double*, double**) : < 1 function");
  
  Buffer2D<double> beta(func.get_nfunc(), 1);
  
  // Intialise
  for(int j=0; j<func.get_nfunc(); j++){
    beta[j][0] = 0.;
    for(int k=0; k<func.get_nfunc(); k++)
      covar[j][k] = 0.;
  }
  
  // Accumulate matrices
  double* v = new double[func.get_nfunc()];
  for(int i=0; i<ndata; i++){
    if(e[i] > 0.){
      double weight = 1./Subs::sqr(e[i]);
      func.eval(x[i],v);
      for(int j=0; j<func.get_nfunc(); j++){
	double wt = weight*v[j];
	for(int k=0; k<=j; k++)
	  covar[j][k] += wt*v[k];
	beta[j][0] += wt*y[i];
      } 
    }
  }
  delete[] v;
  
  // fill rest of matrix by symmetry
  for(int j=1; j<func.get_nfunc(); j++)
    for(int k=0; k<j; k++)
      covar[k][j] = covar[j][k];
  
  // Solve normal equations
  gaussj(covar, func.get_nfunc(), beta, 1);
  
  // return for coefficients
  for(int i=0; i<func.get_nfunc(); i++)
    coeff[i] = beta[i][0];
  
}

/** llsqr_eval calculates the fit to the data after a run of the function object version of llsqr
 * \param ndata the number of data points
 * \param x    the X values
 * \param func the functions in the form of a function object that has member functions of the form
 * 'void eval(double x, double* v)' to evaluate the nfunc function values at x, returning them in v,  and 
 * 'int get_nfunc()' const which returns the number of functions.
 * \param coeff the fitted coefficients, i.e. the nfunc multipliers of the functions that lead to the
 * best (in a least squares sense) fit to the data
 * \param fit the ndata fitted values
 * \return the Chi**2
 */

void Subs::llsqr_eval(int ndata, const double* x, const Llfunc& func, const double* coeff, double* fit){
  
  if(ndata < 1)
    throw Subs_Error("Subs::llsqr_eval(int, const double*, int, const Llfunc&, double*, double*) : < 1 data point");
  if(func.get_nfunc() < 1)
    throw Subs_Error("Subs::llsqr_eval(int, const double*, int, const Llfunc&, double*, double*) : < 1 function");
  
  double* v = new double[func.get_nfunc()];
  for(int i=0; i<ndata; i++){
    func.eval(x[i],v);
    double ymodel = 0.;
    for(int j=0; j<func.get_nfunc(); j++)
      ymodel += coeff[j]*v[j];
    fit[i] = ymodel;
  }
  delete[] v;
}

/** llsqr_chisq evaluates the chi**2 value after a run of the function object version of llsqr
 * \param ndata the number of data points
 * \param x    the X values
 * \param y    the Y values
 * \param e    uncertainties on the Y values (<= 0 to mask)
 * \param func the functions in the form of a function object that has member functions of the form
 * 'void eval(double x, double* v)' to evaluate the nfunc function values at x, returning them in v,  and 
 * 'int get_nfunc()' const which returns the number of functions.
 * \param coeff the fitted coefficients, i.e. the nfunc multipliers of the functions that lead to the
 * best (in a least squares sense) fit to the data
 * \return the Chi**2
 */
double Subs::llsqr_chisq(int ndata, const double* x, const double* y, const float* e, const Llfunc& func, const double* coeff){
  
  if(ndata < 1)
    throw Subs_Error("Subs::llsqr_chisq(int, const double*, const double*, const float*, const Llfunc&, const double*) : < 1 data point");
  if(func.get_nfunc() < 1)
    throw Subs_Error("Subs::llsqr_chisq(int, const double*, const double*, const float*, const Llfunc&, const double*) : < 1 function");
  
  double chisq = 0.;
  double* v = new double[func.get_nfunc()];
  for(int i=0; i<ndata; i++){
    if(e[i] > 0.){
      func.eval(x[i],v);
      double ymodel = 0.;
      for(int j=0; j<func.get_nfunc(); j++)
	ymodel += coeff[j]*v[j];
      chisq += Subs::sqr((y[i]-ymodel)/e[i]);
    }
  }
  delete[] v;
  return chisq;
}

/** llsqr_reduced_chisq evaluates the chi**2 value after a run of the function object version of llsqr
 * \param ndata the number of data points
 * \param x    the X values
 * \param y    the Y values
 * \param e    uncertainties on the Y values (<= 0 to mask)
 * \param func the functions in the form of a function object that has member functions of the form
 * 'void eval(double x, double* v)' to evaluate the nfunc function values at x, returning them in v,  and 
 * 'int get_nfunc()' const which returns the number of functions.
 * \param coeff the fitted coefficients, i.e. the nfunc multipliers of the functions that lead to the
 * best (in a least squares sense) fit to the data
 * \return the Chi**2
 */
double Subs::llsqr_reduced_chisq(int ndata, const double* x, const double* y, const float* e, const Llfunc& func, const double* coeff){
  
  if(ndata < 1)
    throw Subs_Error("Subs::llsqr_reduced_chisq(int, const double*, const double*, const float*, const Llfunc&, const double*) : < 1 data point");
  if(func.get_nfunc() < 1)
    throw Subs_Error("Subs::llsqr_reduced_chisq(int, const double*, const double*, const float*, const Llfunc&, const double*) : < 1 function");
  
  int ndof = - func.get_nfunc();
  double chisq = 0.;
  double* v = new double[func.get_nfunc()];
  for(int i=0; i<ndata; i++){
    if(e[i] > 0.){
      func.eval(x[i],v);
      double ymodel = 0.;
      for(int j=0; j<func.get_nfunc(); j++)
	ymodel += coeff[j]*v[j];
      chisq += Subs::sqr((y[i]-ymodel)/e[i]);
      ndof++;
    }
  }
  delete[] v;
  if(ndof < 1)
    throw Subs_Error("Subs::llsqr_reduced_chisq(int, const double*, const double*, const float*, const Llfunc&, const double*) : < 1 degree of freedom");
  return chisq/ndof;
}

/** llsqr_reject carries out a rejection cycle after a run of the function object version of llsqr. Data points
 * are masked by making error bars negative.
 *
 * \param ndata the number of data points
 * \param x    the X values
 * \param y    the Y values
 * \param e    uncertainties on the Y values, will be modified 
 * \param func the functions in the form of a function object that has member functions of the form
 * 'void eval(double x, double* v)' to evaluate the nfunc function values at x, returning them in v,  and 
 * 'int get_nfunc()' const which returns the number of functions.
 * \param coeff the fitted coefficients, i.e. the nfunc multipliers of the functions that lead to the
 * best (in a least squares sense) fit to the data
 * \param thresh the threshold in terms of sigma for rejection.
 * \param slow true to just reject the worst point, else all points above threshold will go
 * \return number of points rejected.
 */

int Subs::llsqr_reject(int ndata, const double* x, const double* y, float* e, const Llfunc& func, const double* coeff, double thresh,
			 bool slow){

    if(ndata < 1)
      throw Subs_Error("Subs::llsqr_reject(int, const double*, const double*, float*, const Llfunc&, const double*, double, bool) : < 1 data point");
    if(func.get_nfunc() < 1)
      throw Subs_Error("Subs::llsqr_reject(int, const double*, const double*, float*, const Llfunc&, const double*, double, bool) : < 1 function");
    if(thresh <= 0.)
      throw Subs_Error("Subs::llsqr_reject(int, const double*, const double*, float*, const Llfunc&, const double*, double, bool) : thresh <= 0.");

    double reduced_chisq = llsqr_reduced_chisq(ndata, x, y, e, func, coeff);

    int nrej = 0, iworst = -1;
    double worst = -1., limit = sqrt(reduced_chisq)*thresh;
    double* v = new double[func.get_nfunc()];
    for(int i=0; i<ndata; i++){
      if(e[i] > 0.){
	double ymodel = 0.;
	func.eval(x[i],v);
	for(int j=0; j<func.get_nfunc(); j++)
	  ymodel += coeff[j]*v[j];
	double dev = fabs(y[i]-ymodel)/e[i];
	if(dev > limit){
	  if(slow && dev > worst){
	    worst  = dev;
	    iworst = i;
	  }else{
	    e[i] = - e[i];
	    nrej++;
	  }
	}
      }
    }
    delete[] v;
    if(slow && iworst >= 0){
      e[iworst] = -e[iworst];
      nrej = 1;
    }
    return nrej;
}

  
