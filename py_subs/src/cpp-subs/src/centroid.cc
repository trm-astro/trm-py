#include "trm/subs.h"
#include "trm/constants.h"

// Class to compute the derivative and second-derivative of the cross-correlation (for rtsafe)
 
class CCFfd : public Subs::RTfunc {
  
  const float* data;
  const int p1, p2, plo, phi;
  const double sigma;
  
public:
  
  /* Constructor storing fixed data
   * Note that the peak position must be good because it will be used to calculate
   * the range over which to calculate the CCF which will be held fixed.
   */
  CCFfd(const float data[], int p1, int p2, float fwhm, int peak) : 
    data(data), p1(p1), p2(p2), 
    plo(std::max(p1,peak-int(3.*fwhm+1))), phi(std::min(p2,peak+int(3.*fwhm+1))),
    sigma(fwhm/Constants::EFAC)  {}
					   
  // Function operator
  void operator()(double p, double& f, double& d) const { 
    double x, rat, yw;
    f = d = 0.;
    for(int i=plo; i<=phi; i++){
      x    = double(i)-p;
      rat  = Subs::sqr(x/sigma);
      yw   = exp(-0.5*rat)*data[i];
      f   += x*yw;
      d   += (rat-1.)*yw;
    }
  }
};


// Class to compute the value only of the cross-correlation (for dbrent)
class CCFf : public Subs::Sfunc {
  
  const float* data;
  const int p1, p2, plo, phi;
  const double sigma;
  const bool emission;
  
public:
  
  /* Constructor storing fixed data
   * Note that the peak position must be good because it will be used to calculate
   * the range over which to calculate the CCF which will be held fixed.
   */
  CCFf(const float data[], int p1, int p2, float fwhm, int peak, bool emission) : 
    data(data), p1(p1), p2(p2), 
    plo(std::max(p1,peak-int(3.*fwhm+1))), phi(std::min(p2,peak+int(3.*fwhm+1))),
    sigma(fwhm/Constants::EFAC), emission(emission)  {}
					   
  // Function operator
  double operator()(double p) { 
    double x, f = 0.;
    for(int i=plo; i<=phi; i++){
      x    = double(i)-p;
      f   += exp(-0.5*Subs::sqr(x/sigma))*data[i];
    }
    return (emission ? -f : f);
  }
};

// Class to compute the derivative only of the cross-correlation (for dbrent)
class CCFd : public Subs::Sfunc {
  
  const float* data;
  const int p1, p2, plo, phi;
  const double sigma;
  const bool emission;
  
public:
  
  /* Constructor storing fixed data
   * Note that the peak position must be good because it will be used to calculate
   * the range over which to calculate the CCF which will be held fixed.
   */
  CCFd(const float data[], int p1, int p2, float fwhm, int peak, bool emission) : 
    data(data), p1(p1), p2(p2), 
    plo(std::max(p1,peak-int(3.*fwhm+1))), phi(std::min(p2,peak+int(3.*fwhm+1))),
    sigma(fwhm/Constants::EFAC), emission(emission)  {}
					   
  // Function operator
  double operator()(double p) { 
    double x, rat;
    double d = 0.;
    for(int i=plo; i<=phi; i++){
      x    = double(i)-p;
      rat  = x/sigma;
      d   += rat*exp(-0.5*Subs::sqr(rat))*data[i]/sigma;
    }
    return (emission ? -d : d);
  }
};

/**
 * Gaussian centroiding routine. This measures the position of a peak
 * by cross-correlation with a gaussian of fixed width. This is fairly
 * robust in the presence of noise. It first searches for a maximum in the
 * cross-correlation, then runs a Newton-Raphson-like refinement to tweak the value.
 * Take care near the ends of the array: the routine effectively assumes that beyond the
 * ends of the array it carries on = 0. Thus it may be best to remove a background value
 * if you want to avoid possible problems associated with this assumption.
 *
 * \param data     the data array
 * \param var      variances on the data
 * \param p1       the first pixel of the array to consider. Array starts at 0
 * \param p2       the last pixel of the array to consider. Must be at least 2 more than p1.
 * \param fwhm     the FWHM of gaussian. This will be limited to be at least 2 since below this value the routine 
 * locks onto pixel centres.
 * \param start    the starting position
 * \param emission true if the line is a peak rather than a trough
 * \param pos      the position (returned)
 * \param epos     the uncertainty in the returned position (returned)
 * \exception The routine throws exceptions if it fails to find a local maximum of the 
 * cross-correlation. The likelihood of this can be reduced bu subtraction of a constant if the mean
 * level is significantly different from 0.
 */

void Subs::centroid(const float data[], const float var[], int p1, int p2, float fwhm, float start, 
		    bool emission, double& pos, float& epos){
  
  fwhm = std::max(2.f, fwhm);
  const int NCNT = 13;
  double cnt[NCNT];
  const float SIG     = fwhm/Constants::EFAC;
  const int WIDTH     = int(3.*fwhm+1);
  const int SEARCH    = std::max(2, std::min(WIDTH, NCNT/2));

  if(p1 < 0 || p2-p1 < 2)
    throw Subs_Error("Subs::centroid: invalid range " + Subs::str(p1) + " to " + Subs::str(p2));

  // Hunt for a local maximum in cross-correlation. 
  int peak, ichamp = int(start+0.5);
  int ix, cycle = 0;
  double sum, champ;

  epos = -1.;

  do{
    // clamp the presumed maximum to lie within range.
    ichamp = std::max(p1, std::min(p2, ichamp));
    peak   = ichamp;
    cycle++;

    // Give up if too many cycles or if current best pixel is stuck near the end limits
    if(cycle > 10 || (cycle > 1 && (ichamp == p1 || ichamp == p2)))
      throw Subs_Error("Subs::centroid: failed to find a local maximum of the cross-correlation after " + Subs::str(cycle) +
		       " cycles with ichamp = " + Subs::str(ichamp) + " cf range " + Subs::str(p1) + " to " + Subs::str(p2));

    // Evaluate X-corr in range around current best position
    for(int ip=peak-SEARCH; ip<=peak+SEARCH; ip++){
      ix  = ip-peak+SEARCH;
      sum = 0.;
      for(int iq=std::max(p1,ip-WIDTH); iq<=std::min(p2,ip+WIDTH); iq++)
	sum  += exp(-Subs::sqr(float(iq-ip)/SIG)/2.)*data[iq];
      cnt[ix] = sum;
    }

    // Work out point of maximum x-corr (accounting for emission or absorption)
    champ  = cnt[0];
    ichamp = peak-SEARCH;
    for(ix = 1; ix<2*SEARCH+1; ix++){
      if((emission  && cnt[ix] > champ) || (!emission && cnt[ix] < champ)){
	champ  = cnt[ix];
	ichamp = peak - SEARCH + ix;
      }
    }

    // Try again if the maximum is too close to the end of the search range.

  }while(abs(ichamp-peak) == SEARCH);

  // OK have found a maximum/minimum within range. Use 'rtsafe' to refine it.
  const CCFfd ccffd(data, p1, p2, fwhm, ichamp);

  try{
    pos = rtsafe(ccffd, ichamp-1, ichamp+1, 1.e-7);
  }
  catch(const Subs::Subs_Error& err){

    // rtsafe can fail when the gradient at either end of the range has the same sign.
    // try to get round this with dbrent (experiments show this to be much slower but 
    // more robust)

    CCFf ccff(data, p1, p2, fwhm, ichamp, emission);
    CCFd ccfd(data, p1, p2, fwhm, ichamp, emission);
    dbrent(ichamp-1, ichamp, ichamp+1, ccff, ccfd, 1.e-7, false, 0., pos);

  }

  // Calculate error in bin position from variance on each point
  double x1 = 0., x2 = 0., x, rat, epf;
  int plo = std::max(p1,int(pos+0.5)-WIDTH);
  int phi = std::min(p2,int(pos+0.5)+WIDTH);
  for(int i=plo; i<=phi; i++){
    x    = i - pos;
    rat  = Subs::sqr(x/SIG);
    epf  = exp(-rat*0.5);
    x1  += Subs::sqr(x*epf)*var[i];
    x2  += data[i]*(rat-1.)*epf;
  }
  if(x2*x2 == 0.) 
    throw Subs_Error("Subs::centroid: failed to calculate an uncertainty on the position");
  epos = sqrt(x1/(x2*x2));
}

