#include <cmath>
#include <cfloat>
#include "trm/subs.h"

namespace Subs {
    double dmod(double a, double b);
    void spread(double y, double *yy, int n, double x, int m);
}

/** Press & Rybicki's fast method for Lomb-Scargle periodogram adapted to handle
 * uncertainties.
 * \param x X data
 * \param y Y data
 * \param e uncertainties. Set equal to sample RMS for standard normalisation, <= 0 to ignore
 * \param n the number of points
 * \param ofac oversampling factor (typically 4)
 * \param hifac maximum frequency as multiple of 'average Nyquist'
 * \param freq returned frequencies
 * \param pgram returned periodogram
 */
void Subs::fasper(double *x, float *y, float *e,  int n, double ofac, double hifac, Subs::Buffer1D<double>& freq, Subs::Buffer1D<double>& pgram){

    const int MACC = 4;
    const int nfreq  = int(0.5*ofac*hifac*n);
    const int nfreqt = int(ofac*hifac*n*MACC);
    int nf     = 64;

    while(nf < nfreqt) nf <<= 1;
    int ndim = nf << 1;

    try{
        freq.resize(ndim);
        pgram.resize(ndim);
    }
    catch(const Subs_Error& e){
        throw Subs_Error("Error in Subs::fasper: " + e);
    }

    double xmin = DBL_MAX, xmax = -DBL_MAX;
    for(int i=0; i<n; i++){
        if(e[i] > 0.){
            if(x[i] < xmin) xmin = x[i];
            if(x[i] > xmax) xmax = x[i];
        }
    }
    const double xdif = xmax-xmin;  
    freq  = 0;
    pgram = 0;
    const double fac   = ndim/(xdif*ofac);
    const double dndim = double(ndim);

    // "Extirpolate" into large array 
    double wsum = 0., ck, ckk, w;
    for(int i=0; i<n; i++){
        if(e[i] > 0.){
            w     = 1./(e[i]*e[i]);
            wsum += w;
            ck    = dmod((x[i]-xmin)*fac,dndim);
            ckk   = dmod(2.*ck,dndim);
            spread(w*y[i],freq,ndim,ck,MACC);
            spread(w,pgram,ndim,ckk,MACC);
        }
    }

    fftr(freq,ndim,1);
    fftr(pgram,ndim,1);

    double df = 1./(xdif*ofac);
    double cwt, swt, hypo, hc2wt, hs2wt, den, cterm, sterm;
    for(int i=0,k=2; i<nfreq; i++,k+=2){
        hypo     = sqrt(sqr(pgram[k])+sqr(pgram[k+1]));
        hc2wt    = 0.5*pgram[k]/hypo;
        hs2wt    = 0.5*pgram[k+1]/hypo;
        cwt      = sqrt(0.5+hc2wt);
        swt      = sign(sqrt(0.5-hc2wt),hs2wt);
        den      = 0.5*wsum+hc2wt*pgram[k]+hs2wt*pgram[k+1];
        cterm    = sqr(cwt*freq[k]+swt*freq[k+1])/den;
        sterm    = sqr(cwt*freq[k+1]-swt*freq[k])/(wsum-den);
        freq[i]  = df*(i+1);
        pgram[i] = (cterm+sterm)/2.;
    }
    freq.set_npix(nfreq);
    pgram.set_npix(nfreq);
}

/** Press & Rybicki's fast method for Lomb-Scargle periodogram adapted to handle
 * uncertainties. This implementation is to allow a fixed set of frequencies to be used
 * but requires more work from the user to work out suitable inputs.
 * \param x X data
 * \param y Y data
 * \param e uncertainties. Set equal to sample RMS for standard normalisation, <= 0 to ignore
 * \param n the number of points
 * \param fmax maximum frequency, cycles/unit x
 * \param nfreq number of frequencies. 
 * \param freq returned frequencies. Will be returned with nfreq frequencies.
 * \param pgram returned periodogram
 */
void Subs::fasper(double *x, float *y, float *e, int n, double fmax, int nfreq, Subs::Buffer1D<double>& freq, Subs::Buffer1D<double>& pgram){
  
    // Compute valid X range
    double xmin = DBL_MAX, xmax = -DBL_MAX;
    for(int i=0; i<n; i++){
        if(e[i] > 0.){
            if(x[i] < xmin) xmin = x[i];
            if(x[i] > xmax) xmax = x[i];
        }
    }
    const double xdif  = xmax-xmin;
    const double ofac  = nfreq/(xdif*fmax);
    if(ofac  < 1.)
        throw Subs::Subs_Error("void Subs::fasper(double, float, float, int, double, int, Subs::Buffer1D<double>&, Subs::Buffer1D<double>&): " //
                               "oversampling factor = " + Subs::str(ofac) + " and is < 1. Need more frequencies or a smaller maximum frequency.");

    const int MACC   = 4;
    const int nfreqt = 2*MACC*nfreq;

    // Find first power of 2 > 2*nfreqt
    int nf = 64;
    while(nf < nfreqt) nf <<= 1;
    const int ndim = nf << 1;

    try{
        freq.resize(ndim);
        pgram.resize(ndim);
    }
    catch(const Subs_Error& e){
        throw Subs_Error("Error in Subs::fasper: " + e);
    }

    freq  = 0;
    pgram = 0;

    const double fac   = ndim*(fmax/nfreq);
    const double dndim = double(ndim);

    // "Extirpolate" into large array 
    double wsum = 0., ck, ckk, w;
    for(int i=0; i<n; i++){
        if(e[i] > 0.){
            w     = 1./(e[i]*e[i]);
            wsum += w;
            ck    = dmod(fac*(x[i]-xmin), dndim);
            ckk   = dmod(2.*ck, dndim);
            spread(w*y[i],freq,ndim,ck,MACC);
            spread(w,pgram,ndim,ckk,MACC);
        }
    }

    // Carry out FFTs
    fftr(freq,ndim,1);
    fftr(pgram,ndim,1);

    double df = fmax/nfreq;
    double cwt, swt, hypo, hc2wt, hs2wt, den, cterm, sterm;
    for(int i=0,k=2; i<nfreq; i++,k+=2){
        hypo     = sqrt(sqr(pgram[k])+sqr(pgram[k+1]));
        hc2wt    = 0.5*pgram[k]/hypo;
        hs2wt    = 0.5*pgram[k+1]/hypo;
        cwt      = sqrt(0.5+hc2wt);
        swt      = sign(sqrt(0.5-hc2wt),hs2wt);
        den      = 0.5*wsum+hc2wt*pgram[k]+hs2wt*pgram[k+1];
        cterm    = sqr(cwt*freq[k]+swt*freq[k+1])/den;
        sterm    = sqr(cwt*freq[k+1]-swt*freq[k])/(wsum-den);
        freq[i]  = df*(i+1);
        pgram[i] = (cterm+sterm)/2.;
    }
    freq.set_npix(nfreq);
    pgram.set_npix(nfreq);
}

/** Press & Rybicki's fast method to calculate amplitude spectra resulting from a
 * a least-squares fit of a sinusoid at each frequency. 
 * \param x X data
 * \param y Y data
 * \param e uncertainties. <= 0 to ignore
 * \param n the number of points
 * \param fmax maximum frequency, cycles/unit x
 * \param nfreq number of frequencies. 
 * \param freq returned frequencies
 * \param amps returned amplitudes
 */
void Subs::famp(double *x, float *y, float *e,  int n, double fmax, int nfreq, Subs::Buffer1D<double>& freq, Subs::Buffer1D<double>& amps){
  
    // Compute valid X range
    double xmin = DBL_MAX, xmax = -DBL_MAX;
    for(int i=0; i<n; i++){
        if(e[i] > 0.){
            if(x[i] < xmin) xmin = x[i];
            if(x[i] > xmax) xmax = x[i];
        }
    }
    const double xdif  = xmax-xmin;
    const double ofac  = nfreq/(xdif*fmax);
    if(ofac  < 1.)
        throw Subs::Subs_Error("void Subs::famp(double, float, float, int, double, int, Subs::Buffer1D<double>&, Subs::Buffer1D<double>&): " //
                               "oversampling factor = " + Subs::str(ofac) + " and is < 1. Need more frequencies or a smaller maximum frequency.");
  
    const int MACC = 4;
    const int nfreqt = 2*MACC*nfreq;
    int nf     = 64;

    while(nf < nfreqt) nf <<= 1;
    int ndim = nf << 1;

    try{
        freq.resize(ndim);
        amps.resize(ndim);
    }
    catch(const Subs_Error& e){
        throw Subs_Error("Error in Subs::famp: " + e);
    }

    freq = 0;
    amps = 0;
    const double fac   = ndim/(xdif*ofac);
    const double dndim = double(ndim);

    // "Extirpolate" into large array 
    double wsum = 0., ck, ckk, w;
    for(int i=0; i<n; i++){
        if(e[i] > 0.){
            w     = 1./(e[i]*e[i]);
            wsum += w;
            ck    = dmod((x[i]-xmin)*fac,dndim);
            ckk   = dmod(2.*ck,dndim);
            spread(w*y[i],freq,ndim,ck,MACC);
            spread(w,amps,ndim,ckk,MACC);
        }
    }

    fftr(freq,ndim,1);
    fftr(amps,ndim,1);

    double df = 1./(xdif*ofac);
    double hypo, hc2wt, hs2wt;
    for(int i=0,k=2; i<nfreq; i++,k+=2){

        hypo     = sqrt(sqr(amps[k])+sqr(amps[k+1]));
        hc2wt    = 0.5*amps[k]/hypo;
        hs2wt    = 0.5*amps[k+1]/hypo;

        double ac = (0.5*wsum-hc2wt)*freq[k] - hs2wt*freq[k+1];
        double as = -hs2wt*freq[k] + (0.5*wsum+hc2wt)*freq[k+1];

        freq[i]  = df*(i+1);
        amps[i]  = sqrt(ac*ac + as*as)/(sqr(0.5*wsum) - sqr(hc2wt) - sqr(hs2wt));

    }
    freq.set_npix(nfreq);
    amps.set_npix(nfreq);
}

namespace Subs {
    void spread(double y, double *yy, int n, double x, int m){
	
        const int NFAC[] = {1,1,2,6,24,120,720,5040,40320,362880};
	
        if(m >= 9)
            throw Subs::Subs_Error("factorial table too small in spread inside Subs::fasper");
	
        int ix = int(x);
	
        if(x == double(ix))
            yy[ix] += y;
        else{
            int ilo  = std::min(std::max(int(x-0.5*m),0),int(n-m));
            int ihi  = ilo + m - 1;
            int nden = NFAC[m-1];
            double fac  = x - ilo;
            for(int j=ilo+1; j<=ihi; j++) fac *= (x-j);
            yy[ihi] += y*fac/(nden*(x-ihi));
            for(int j=ihi-1; j>=ilo; j--){
                nden = (nden/(j+1-ilo))*(j-ihi);
                yy[j] += y*fac/(nden*(x-j));
            }
        }
    }
    
    /** Removes b from a until a is less than b
     */
    
    double dmod(double a, double b){
        while(a >= b) a -= b;
        return a;
    }
}

