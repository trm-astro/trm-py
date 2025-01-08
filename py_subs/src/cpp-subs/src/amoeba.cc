#include "trm/subs.h"
#include "trm/array1d.h"

namespace Subs {
  double amoeba_try(std::vector<std::pair<Subs::Array1D<double>, double> >& params, Subs::Array1D<double>& psum, 
		    Afunc& func, int ihi, double fac);

  void amoeba_get_psum(const std::vector<std::pair<Subs::Array1D<double>, double> >& params, Subs::Array1D<double>& psum);
}

/** Simplex minimisation routine. Start from N+1 cornered object where N is number
 * of variables.
 * \param params n+1 pairs, each with an Subs::Array1D defining the corner and the function value at that corner
 * \param ftol fractional toleranc
xe on minimum. On exit all corners should have values within fraction ftol of each other 
 * \param func function object which returns chi**2 when called as func(vector<double>& vec) where vec are the variable 
 * parameter values. This should be defined by inheritance from the Vfunc class ('Vector function'). See subs.h
 * \param nfunc number of calls
 */
void Subs::amoeba(std::vector<std::pair<Subs::Array1D<double>, double> >& params, double ftol, int nmax, Afunc& func, int& nfunc){
    
  const double TINY = 1.e-10;
  
  for(size_t i=0; i<params.size(); i++)
    if(int(params.size()) != params[i].first.size() + 1)
      throw Subs_Error("amoeba: number of parameters is not one less than the number of parameter vectors on input");
  
  int ndim = params[0].first.size();
  Subs::Array1D<double> psum(ndim);
  nfunc = 0;
  
  amoeba_get_psum(params, psum);
  
  for(;;){
    int ilo = 0;
    int inhi;
    int ihi = params[0].second > params[1].second ? (inhi=2,1) : (inhi=1,2);
    for(size_t i=0; i<params.size(); i++){
      if(params[i].second <= params[ilo].second) ilo = i;
      if(params[i].second > params[ihi].second) {
	inhi = ihi;
	ihi  = i;
      }else if(params[i].second > params[inhi].second && int(i) != ihi) {
	inhi = i;
      }
    }
    double rtol = 2.*fabs(params[ihi].second-params[ilo].second)/(fabs(params[ihi].second)+fabs(params[ilo].second)+TINY);
    
    if(rtol < ftol){
      std::swap(params[0].second,params[ilo].second);
      for(int i=0; i<ndim; i++)
	std::swap(params[0].first[i], params[ilo].first[i]);
      break;
    }
    
    if(nfunc >= nmax){
      std::cerr << "amoeba: nmax = " << nmax << " exceeded." << std::endl;
      break;
    }
    
    nfunc += 2;
    
    double ytry = amoeba_try(params, psum, func, ihi, -1.0);
    if(ytry <= params[ilo].second){
      ytry = amoeba_try(params, psum, func, ihi, 2.0);
    }else if(ytry >= params[inhi].second){
      double ysave = params[ihi].second;
      ytry = amoeba_try(params, psum, func, ihi, 0.5);
      if(ytry >= ysave){
	for(size_t i=0; i<params.size(); i++){
	  if(int(i) != ilo){
	    for(int j=0; j<ndim; j++)
	      params[i].first[j] = psum[j] = 0.5*(params[i].first[j]+params[ilo].first[j]);
	    params[i].second = func(psum);
	  }
	}
	nfunc += ndim;
	amoeba_get_psum(params, psum);
      }
    }else{
      nfunc--;
    }
  }
}

namespace Subs {

    //! Helper routine for amoeba
    void amoeba_get_psum(const std::vector<std::pair<Subs::Array1D<double>, double> >& params, 
			 Subs::Array1D<double>& psum){
	for(int j=0; j<psum.size(); j++){
	    double sum = 0.;
	    for(size_t i=0; i<params.size(); i++)
		sum += params[i].first[j];
	    psum[j] = sum;
	}
    }

    //! Helper routine for amoeba
    double amoeba_try(std::vector<std::pair<Subs::Array1D<double>, double> >& params, 
		      Subs::Array1D<double>& psum, Subs::Afunc& func, int ihi, double fac){
	int ndim = params[0].first.size();
	Subs::Array1D<double> ptry(ndim);
	double fac1 = (1.-fac)/ndim;
	double fac2 = fac1-fac;
	for(int j=0; j<ndim; j++)
	    ptry[j] = psum[j]*fac1-params[ihi].first[j]*fac2;
	double ytry = func(ptry);
	if(ytry < params[ihi].second){
	    params[ihi].second = ytry;
	    for(int j=0; j<ndim; j++){
		psum[j] += ptry[j] - params[ihi].first[j];
		params[ihi].first[j] = ptry[j];
	    }
	}
	return ytry;
    }
}
    


