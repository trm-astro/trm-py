// Fits multiple sines and cosines of fixed frequency to a
// data vector. 
//
// double x[ndata]  -- x array (input)
// float  y[ndata]  -- y array (input)
// float  e[ndata]  -- error array (-ve to ignore) (input)
// long unsigned int ndata -- number of points (input)
// double freq[nfreq] -- array of frequencies (input)
// long unsigned int nfreq -- number of frequencies (input)
// double a[nfreq]    -- cosine amplitudes (returned)
// double b[nfreq]    -- sine amplitudes (returned)
// double ea[nfreq]   -- uncertainties on cosine amplitudes (returned)
// double eb[nfreq]   -- uncertainties on sine amplitudes (returned)

#include <math.h>
#include "trm/subs.h"

void Subs::sincos(double *x, float *y, float *e, long unsigned int ndata, 
		 double *freq, long unsigned int nfreq, double *a, double *b,
		 double *ea, double *eb, double *rab){
  
  const double twopi = 8.*atan(1.);
  long unsigned int i, j, n, ndim;

  for(i=0;i<nfreq;i++){
    a[i]   = 0.;
    b[i]   = 0.;
    ea[i]  = -1.;
    eb[i]  = -1.;
    rab[i] =  0.;
  }
  ndim  = 2*nfreq;

  // Get workspace

  double v[ndim];
  double **M = new double *[ndim];
  for(i=0;i<ndim;i++){
    v[i] = 0.;
    M[i] = new double [ndim];
    for(j=0;j<ndim;j++)
      M[i][j] = 0.;
  }

  double ysq = 0., xx, temp;
  float ww;
  long unsigned int np = 0;
  long unsigned int i1, i2, j1, j2; 
  double cni, sni, cnj, snj;
  
  for(n=0;n<ndata;n++){
    if(e[n] > 0.){
      np++;
      temp = twopi*x[n];
      ww   = 1./(e[n]*e[n]);
      ysq += ww*y[n]*y[n];
      for(i=0;i<nfreq;i++){
	xx = temp*freq[i];
	cni = cos(xx);
	sni = sin(xx);

	// Accumulate sums

	i1     = 2*i;
	i2     = i1+1;
	v[i1] += ww*cni*y[n];
	v[i2] += ww*sni*y[n];
	for(j=i;j<nfreq;j++){
	  xx  = temp*freq[j];
	  cnj = cos(xx);
	  snj = sin(xx);
	  j1  = 2*j;
	  j2  = j1 + 1;
	  M[i1][j1] += ww*cni*cnj;
          M[i1][j2] += ww*cni*snj;
	  M[i2][j2] += ww*sni*snj;
	}
      }
    }
  }
  if(np<2)
    throw Subs_Error("Too few points to fit inside sincos");

  // Symmetrise matrix M

  for(i=1;i<ndim;i++)
    for(j=0;j<i;j++)
      M[i][j] = M[j][i];

  // Solve simultaneous equations to determine fit coefficients

  unsigned int indx[ndim];
  double det;
  double unit[ndim];
  double IM[ndim][ndim];

  if(ludcmp(M,ndim,indx,det)){
    lubksb(M,ndim,indx,v);

    // Store fit coefficients

    for(i=0;i<nfreq;i++){
      a[i] = v[2*i];
      b[i] = v[2*i+1];
    }

    // Invert matrix to get uncertainties

    for(i=0;i<ndim;i++){
      for(j=0;j<ndim;j++)
	unit[j] = 0.;
      unit[i] = 1.;
      lubksb(M,ndim,indx,unit);
      for(j=0;j<ndim;j++)
	IM[i][j] = unit[j];
    }
      
    // Store errors on A and B and correlation coefficients

    for(i=0;i<nfreq;i++){
      i1 = 2*i;
      i2 = i1 + 1;
      ea[i]  = sqrt(IM[i1][i1]);
      eb[i]  = sqrt(IM[i2][i2]);
      rab[i] = IM[i2][i1]/sqrt(IM[i2][i2]*IM[i1][i1]);
    }
  }
  delete M;
}
