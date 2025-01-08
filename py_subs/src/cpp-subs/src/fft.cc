#include <cstdlib>
#include <cmath>
#include <iostream>
#include "trm/subs.h"
#include "trm/constants.h"

/**  FFT routine. 
 *  
 * \param data the data array. This must have alternate real and imaginary components (see NR for more
 *  details).
 * \param nump the number of points. This must be a power of 2; zero pad if necessary
 * \param flag flag  *  If flag = 1, data[nump] is replaced by it's discrete fourier transform.
 *  If flag =-1, data[nump] is replaced by it's inverse FT.
 */

void Subs::fft(float* data, long unsigned int nump, int flag) {

  unsigned long m, mmax, i, j, istep;
  double wtemp, wr, wpr, wpi, wi, theta;
  float rtemp, itemp;

  // A few stupidity checks...
  if((flag != 1) && (flag != -1)) throw Subs_Error("Input flag to Subs::fft must be +/-1!");

  // Keep shifting i until the msb goes off the end.
  // If there is anything left, this wasn't 2^(int)!        

  i = nump << 1;
  j = nump;
  while(i > j) {
    j = i;
    i <<= 1;
  }
  if(i != 0) throw Subs_Error(Subs::str(nump) + " is not an integer power of 2!");

  // because we are manipulating the data alot, it is
  // more efficient to do in a normal array.
  // Bit reversal (part of FFT algorithm)

  j=1;
  for(i=1; i<nump; i+=2) {
    if(j > i) {
      std::swap(data[j-1], data[i-1]);
      std::swap(data[j],   data[i]);
    }
    m = nump >> 1;
    while((m >= 2) && (j > m)) {
      j -= m;
      m >>= 1;
    }
    j += m;
  }

  // Blame Daniel and Lanczos for this bit...
  // Some kind of trigonometric recurrence (must be done
  // in double precision) 

  mmax = 2;
  while(nump > mmax) {
    istep = mmax << 1;
    theta = flag * (Constants::TWOPI/(double)mmax); /* Inverse or forward? */
    wtemp = sin(0.5 * theta);
    wpr = -2.0 * sqr(wtemp);
    wpi = sin(theta);
    wr = 1.0;
    wi = 0.0;
    for(m=1; m<mmax; m+=2) {
      for(i=m; i<=nump; i+=istep) {
	j = i+mmax;
	rtemp = wr*data[j-1] - wi*data[j];
	itemp = wr*data[j] + wi*data[j-1];
	data[j-1] = data[i-1] - rtemp;
	data[j] = data[i] - itemp;
	data[i-1] += rtemp;
	data[i] += itemp;
      }
      wtemp = wr;
      wr = wr*wpr - wi*wpi + wr;
      wi = wi*wpr + wtemp*wpi + wi;
    }
    mmax = istep;
  }
} 

/** FFT routine for real data. Output data is in real, imag, real, imag form as per 
 * normal except for data[1] which is f_(N/2) since both this and f_0 are real, not complex.
 *
 * \param data real array of data
 * \param nump the number of points, a power of 2
 * \param flag flag=1, data[nump] is replaced by it's discrete fourier transform.
 * If flag=-1, data[nump] is replaced by it's inverse FT .(Multiply by 2/n to get
 * back original data)
 */

void Subs::fftr(float* data, unsigned long int nump, int flag) {

  unsigned long int i,i1,i2,i3,i4,np1;
  float c1=0.5,c2,h1r,h1i,h2r,h2i;
  double wr,wi,wpr,wpi,wtemp,theta;
 
  theta = Constants::PI/(double) (nump>>1);
  if (flag == 1) {
    c2 = -0.5;
    fft(data,nump,1);
  } else {
    c2=0.5;
    theta = -theta;
  }
  wtemp=sin(0.5*theta);
  wpr = -2.0*wtemp*wtemp;
  wpi=sin(theta);
  wr=1.0+wpr;
  wi=wpi;
  np1=nump+1;
  for (i=1;i<(nump>>2);i++) {
    i4=1+(i3=np1-(i2=1+(i1=i+i)));
    h1r=c1*(data[i1]+data[i3]);
    h1i=c1*(data[i2]-data[i4]);
    h2r = -c2*(data[i2]+data[i4]);
    h2i=c2*(data[i1]-data[i3]);
    data[i1]=h1r+wr*h2r-wi*h2i;
    data[i2]=h1i+wr*h2i+wi*h2r;
    data[i3]=h1r-wr*h2r+wi*h2i;
    data[i4] = -h1i+wr*h2i+wi*h2r;
    wr=(wtemp=wr)*wpr-wi*wpi+wr;
    wi=wi*wpr+wtemp*wpi+wi;
  }
  if (flag == 1) {
    data[0] = (h1r=data[0])+data[1];
    data[1] = h1r-data[1];
  } else {
    data[0]=c1*((h1r=data[0])+data[1]);
    data[1]=c1*(h1r-data[1]);
    fft(data,nump,-1);
  }
} 

/**  FFT routine. 
 *  
 * \param data the data array. This must have alternate real and imaginary components (see NR for more
 *  details).
 * \param nump the number of points. This must be a power of 2; zero pad if necessary
 * \param flag flag  *  If flag = 1, data[nump] is replaced by it's discrete fourier transform.
 *  If flag =-1, data[nump] is replaced by it's inverse FT. 
 */

void Subs::fft(double* data, long unsigned int nump, int flag) {

  unsigned long m, mmax, i, j, istep;
  double wtemp, wr, wpr, wpi, wi, theta;
  double rtemp, itemp;

  // A few stupidity checks...

  if((flag != 1) && (flag != -1)) throw Subs_Error("Input flag to Subs::fft must be +/-1!");

  // Keep shifting i until the msb goes off the end.
  // If there is anything left, this wasn't 2^(int)!        

  i = nump << 1;
  j = nump;
  while(i > j) {
    j = i;
    i <<= 1;
  }
  if(i != 0) throw Subs_Error(Subs::str(nump) + " is not an integer power of 2!");

  // because we are manipulating the data a lot, it is
  // more efficient to do in a normal array.
  // Bit reversal (part of FFT algorithm)

  j=1;
  for(i=1; i<nump; i+=2) {
    if(j > i) {
      std::swap(data[j-1], data[i-1]);
      std::swap(data[j],   data[i]);
    }
    m = nump >> 1;
    while((m >= 2) && (j > m)) {
      j -= m;
      m >>= 1;
    }
    j += m;
  }

  // Blame Daniel and Lanczos for this bit...
  // Some kind of trigonometric recurrence (must be done
  // in double precision) 

  mmax = 2;
  while(nump > mmax) {
    istep = mmax << 1;
    theta = flag * (Constants::TWOPI/(double)mmax); /* Inverse or forward? */
    wtemp = sin(0.5 * theta);
    wpr = -2.0 * sqr(wtemp);
    wpi = sin(theta);
    wr = 1.0;
    wi = 0.0;
    for(m=1; m<mmax; m+=2) {
      for(i=m; i<=nump; i+=istep) {
	j = i+mmax;
	rtemp = wr*data[j-1] - wi*data[j];
	itemp = wr*data[j] + wi*data[j-1];
	data[j-1] = data[i-1] - rtemp;
	data[j] = data[i] - itemp;
	data[i-1] += rtemp;
	data[i]   += itemp;
      }
      wtemp = wr;
      wr = wr*wpr - wi*wpi + wr;
      wi = wi*wpr + wtemp*wpi + wi;
    }
    mmax = istep;
  }
} 

/** FFT routine for real data. Output data is in real, imag, real, imag form as per 
 * normal except for data[1] which is f_(N/2) since both this and f_0 are real, not complex.
 *
 * \param data real array of data
 * \param nump the number of points, a power of 2
 * \param flag flag=1, data[nump] is replaced by it's discrete fourier transform.
 * If flag=-1, data[nump] is replaced by it's inverse FT. (Multiply by 2/n to get
 * back original data)
 */

void Subs::fftr(double* data, unsigned long int nump, int flag) {

  unsigned long int i,i1,i2,i3,i4,np1;
  double c1=0.5,c2,h1r,h1i,h2r,h2i;
  double wr,wi,wpr,wpi,wtemp,theta;
 
  theta = Constants::PI/(double) (nump>>1);
  if (flag == 1) {
    c2 = -0.5;
    fft(data,nump,1);
  } else {
    c2=0.5;
    theta = -theta;
  }
  wtemp=sin(0.5*theta);
  wpr = -2.0*wtemp*wtemp;
  wpi = sin(theta);
  wr  = 1.0+wpr;
  wi  = wpi;
  np1=nump+1;
  for (i=1;i<(nump>>2);i++) {
    i4=1+(i3=np1-(i2=1+(i1=i+i)));
    h1r =  c1*(data[i1]+data[i3]);
    h1i =  c1*(data[i2]-data[i4]);
    h2r = -c2*(data[i2]+data[i4]);
    h2i =  c2*(data[i1]-data[i3]);
    data[i1] = h1r+wr*h2r-wi*h2i;
    data[i2] = h1i+wr*h2i+wi*h2r;
    data[i3] = h1r-wr*h2r+wi*h2i;
    data[i4] = -h1i+wr*h2i+wi*h2r;
    wr=(wtemp=wr)*wpr-wi*wpi+wr;
    wi=wi*wpr+wtemp*wpi+wi;
  }
  if (flag == 1) {
    data[0] = (h1r=data[0])+data[1];
    data[1] = h1r-data[1];
  } else {
    data[0]=c1*((h1r=data[0])+data[1]);
    data[1]=c1*(h1r-data[1]);
    fft(data,nump,-1);
  }
} 


/** FFT of 2 real functions simultaneously. This is useful for convolutions.
 *
 * \param data1 first array of data, of length n
 * \param data2 second array of data, of length n
 * \param fft1  first FFT, of length 2*n
 * \param fft2  second FFT, of length 2*n
 * \param n     the number of points, a power of 2
 */

void Subs::twofft(double data1[], double data2[], double fft1[], double fft2[], unsigned long int n) {

  // Pack two arrays into one complex array
  for(unsigned long int j=0, jj=1; j<n; j++, jj+=2){
    fft1[jj-1] = data1[j];
    fft2[jj]   = data2[j];
  }

  fft(fft1, n, 1);

  fft2[0] = fft1[1];
  fft1[1] = fft2[1] = 0.;

  unsigned long int nn2 = n + n;
  unsigned long int nn3 = 1 + nn2;

  for(unsigned long int j=2; j<=n; j+=2){

    double rep = 0.5*(fft1[j]+fft1[nn2-j]);
    double rem = 0.5*(fft1[j]-fft1[nn2-j]);
    double aip = 0.5*(fft1[j+1]+fft1[nn3-j]);
    double aim = 0.5*(fft1[j+1]-fft1[nn3-j]);

    fft1[j]     =  rep;
    fft1[j+1]   =  aim;
    fft1[nn2-j] =  rep;
    fft1[nn3-j] = -aim;

    fft2[j]     =  aip;
    fft2[j+1]   = -rem;
    fft2[nn2-j] =  aip;
    fft2[nn3-j] =  rem;

  }
}


