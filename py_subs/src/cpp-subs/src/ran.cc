#include "trm/subs.h"

namespace Subs {
    void psdes(UINT4& lword, UINT4& rword);
}

/** Random number generator 'ran1' from Numerical Recipes
 * \param seed seed integer, -ve to initialise.
 * \return uniform random number between 0 and 1
 */
double Subs::ran1(INT4 &seed){
  const INT4 IA    = 16807;
  const INT4 IM    = 2147483647;
  const double    AM   = 1.0/IM;
  const INT4 IQ    = 127773;
  const INT4 IR    = 2836;
  const INT4 NTAB  = 32;
  const INT4 NDIV  = 1+(IM-1)/NTAB;
  const double    EPS  = 1.2e-7;
  const double    RNMX = 1.0-EPS;

  int j;
  INT4 k;
  static bool first = true;
  static INT4 iy=0;
  static INT4 iv[NTAB];
  double temp;

  if(seed <= 0 || first){
    first = false;
    if(seed == 0){
      seed = 1;
    }else if(seed < 0){
      seed = -seed;
    }
    for(j=NTAB+7;j>=0;j--){
      k = seed/IQ;
      seed = IA*(seed-k*IQ)-IR*k;
      if(seed < 0) seed += IM;
      if(j < NTAB) iv[j] = seed;
    }
    iy = iv[0];
  }
  k = seed/IQ;
  seed = IA*(seed-k*IQ)-IR*k;
  if(seed < 0) seed += IM;
  j  = iy/NDIV;
  iy = iv[j];
  iv[j] = seed;
  if((temp=AM*iy)> RNMX) return RNMX;
  else return temp;
}

/** Random number generator 'ran2' from Numerical Recipes
 * \param seed seed integer, -ve to initialise.
 * \return uniform random number between 0 and 1
 */
double Subs::ran2(INT4 &seed){
  const INT4 IM1   = 2147483563;
  const INT4 IM2   = 2147483399;
  const double    AM    = 1.0/IM1;
  const INT4 IMM1  = IM1-1;
  const INT4 IA1   = 40014;
  const INT4 IA2   = 40692;
  const INT4 IQ1   = 53668;
  const INT4 IQ2   = 52774;
  const INT4 IR1   = 12211;
  const INT4 IR2   = 3791;
  const INT4 NTAB  = 32;
  const INT4 NDIV  = 1+IMM1/NTAB;
  const double    EPS   = 1.2e-7;
  const double    RNMX  = 1.0-EPS;

  int j;
  long k;
  static bool first = true;
  static long seed2=123456789;
  static long iy=0;
  static long iv[NTAB];
  double temp;

  if(seed <= 0 || first){
    first = false;
    if(seed == 0){
      seed = 1;
    }else if(seed < 0){
      seed = -seed;
    }
    seed2 = seed;
    for(j=NTAB+7; j>=0; j--){
      k    = seed/IQ1;
      seed = IA1*(seed-k*IQ1)-k*IR1;
      if(seed < 0) seed += IM1;
      if(j < NTAB) iv[j] = seed;
    }
    iy = iv[0];
  }
  k     = seed/IQ1;
  seed  = IA1*(seed-k*IQ1)-k*IR1;
  if(seed < 0) seed += IM1;
  k     = seed2/IQ2;
  seed2 = IA2*(seed2-k*IQ2)-k*IR2;
  if(seed2 < 0) seed2 += IM2;
  j = iy/NDIV;
  iy = iv[j] - seed2;
  iv[j] = seed;
  if(iy < 1) iy += IMM1;
  if((temp=AM*iy) > RNMX){
    return RNMX;
  }else{
    return temp;
  }
}
  
/** Random number generator 'ran3' from Numerical Recipes
 * \param seed seed integer, -ve to initialise.
 * \return uniform random number between 0 and 1
 */
double Subs::ran3(INT4& seed){
  const INT4 MBIG  = 1000000000L;
  const INT4 MSEED = 161803398L;

  static int inext, inextp;
  static INT4 ma[56];
  static bool first = true;
  INT4 mj, mk;

  if(seed < 0 || first){
    first = false;
    mj = Subs::abs(MSEED-Subs::abs(seed));
    mj %= MBIG;
    ma[55] = mj;
    mk = 1;
    int ii;
    for(int i=1; i<=54; i++){
      ii = (21*i) % 55;
      ma[ii] = mk;
      mk = mj - mk;
      if(mk < 0) mk += MBIG;
      mj = ma[ii];
    }
    for(int k=1; k<=4; k++){
      for(int i=1; i<=55; i++){
	ma[i] -= ma[1+(i+30) % 55];
	if(ma[i] < 0) ma[i] += MBIG;
      }
    }
    inext = 0;
    inextp = 31;
    seed = 1;
  }

  if(++inext == 56)  inext = 1;
  if(++inextp == 56) inextp = 1;
  mj = ma[inext] - ma[inextp];
  if(mj < 0) mj += MBIG;
  ma[inext] = mj;
  return mj/double(MBIG);
}

/** Random number generator 'ran4' from Numerical Recipes;
 *  MUST be 32-bit integers (INT4).
 * \param seed seed integer, -ve to initialise.
 * \return uniform random number between 0 and 1
 */

double Subs::ran4(INT4& seed){

  static long idums = 0;
  UINT4 rword, lword;

  // For 32-bit integers only!!!
  const double NORMALISE  = pow(2.,-32);

  if(seed < 0){
    idums = -seed;
    seed  = 1;
  }
  
  rword = seed;
  lword = idums;
  
  psdes(lword, rword);
  seed++;

  // NR do some fancy machine dependent masking here for reasons I
  // don't quite follow; I could see no obvious difference in speed 
  // and the price is portability. I think this should suffice
  return NORMALISE*rword;

}

namespace Subs {

    void psdes(Subs::UINT4& lword, Subs::UINT4& rword){
	unsigned long i, ia, ib, iswap, itmph=0, itmpl=0;
	const Subs::UINT4 NITER = 4;
	static Subs::UINT4 c1[NITER] = {
	    0xbaa96887L, 0x1e17d32cL, 0x03bcdc3cL, 0x0f33d1b2L};
	static Subs::UINT4 c2[NITER] = {
	    0x4b0f3b58L, 0xe874f0c3L, 0x6955c5a6L, 0x55a7ca46L};
	
	for(i=0;i<NITER;i++){
	    ia     = (iswap = rword) ^ c1[i];
	    itmpl  = ia & 0xffff;
	    itmph  = ia >> 16;
	    ib     = itmpl*itmpl + ~(itmph*itmph);
	    rword  = lword ^ (((ia = (ib >> 16) | ((ib & 0xffff) << 16)) ^ c2[i]) + itmpl*itmph);
	    lword  = iswap;
	}
    }   
}
