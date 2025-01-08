#include <math.h>
#include "trm/subs.h"

/** Generates gaussian random deviates, mean 0, rms 1
 * Uses ran1 as the generator.
 * \param seed seed integer
 * \return Gaussian random number
 */
double Subs::gauss1(INT4 &seed){

  static int iset=0;
  static double gset;
  double fac, rsq, vv1, vv2;

  if(seed < 0) iset = 0;
  if(iset == 0){
    do{
      vv1 = 2.0f*ran1(seed)-1.0f;
      vv2 = 2.0f*ran1(seed)-1.0f;
      rsq = vv1*vv1+vv2*vv2;
    }while(rsq >= 1.0f || rsq == 0.0f);
    fac = sqrt(-2.0f*log(rsq)/rsq);
    gset = vv1*fac;
    iset = 1;
    return vv2*fac;
  } else{
    iset = 0;
    return gset;
  }
}

/** Generates gaussian random deviates, mean 0, rms 1
 * Uses ran2 as the generator.
 * \param seed seed integer
 * \return Gaussian random number
 */
double Subs::gauss2(INT4 &seed){

  static int iset=0;
  static double gset;
  double fac, rsq, vv1, vv2;

  if(seed < 0) iset = 0;
  if(iset == 0){
    do{
      vv1 = 2.0*ran2(seed)-1.0;
      vv2 = 2.0*ran2(seed)-1.0;
      rsq = vv1*vv1+vv2*vv2;
    }while(rsq >= 1.0 || rsq == 0.0);
    fac = sqrt(-2.0*log(rsq)/rsq);
    gset = vv1*fac;
    iset = 1;
    return vv2*fac;
  } else{
    iset = 0;
    return gset;
  }
}

/** Generates gaussian random deviates, mean 0, rms 1
 * Uses ran3 as the generator.
 * \param seed seed integer
 * \return Gaussian random number
 */
double Subs::gauss3(INT4 &seed){

  static int iset=0;
  static double gset;
  double fac, rsq, vv1, vv2;

  if(seed < 0) iset = 0;
  if(iset == 0){
    do{
      vv1 = 2.0*ran3(seed)-1.0;
      vv2 = 2.0*ran3(seed)-1.0;
      rsq = vv1*vv1+vv2*vv2;
    }while(rsq >= 1.0 || rsq == 0.0);
    fac = sqrt(-2.0*log(rsq)/rsq);
    gset = vv1*fac;
    iset = 1;
    return vv2*fac;
  } else{
    iset = 0;
    return gset;
  }
}

/** Generates gaussian random deviates, mean 0, rms 1
 * Uses ran4 as the generator.
 * \param seed seed integer
 * \return Gaussian random number
 */
double Subs::gauss4(INT4 &seed){

  static int iset=0;
  static double gset;
  double fac, rsq, vv1, vv2;

  if(seed < 0) iset = 0;
  if(iset == 0){
    do{
      vv1 = 2.0*ran4(seed)-1.0;
      vv2 = 2.0*ran4(seed)-1.0;
      rsq = vv1*vv1+vv2*vv2;
    }while(rsq >= 1.0 || rsq == 0.0);
    fac = sqrt(-2.0*log(rsq)/rsq);
    gset = vv1*fac;
    iset = 1;
    return vv2*fac;
  } else{
    iset = 0;
    return gset;
  }
}




