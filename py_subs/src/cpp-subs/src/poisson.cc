#include <cmath>
#include "trm/constants.h"
#include "trm/subs.h"

/** Generate numbers distributed as Poisson random numbers
 * of mean mu using the rejection method. I have verified
 * that this routine works. Uses ran1.
 *
 * \param mu mean value
 * \param seed seed integer, set negative to initialise.
 */
float Subs::poisson1(float mu, INT4& seed){

  static float sq, alxm,g,oldm=-1.;
  float em,t,y;

  if(mu < 12.){ 
    if(mu != oldm){
      oldm = mu;
      g = exp(-mu);
    }
    em = -1.;
    t = 1.0;
    do {
      ++em;
      t *= ran1(seed);
    } while(t > g);
  }else{
    if(mu != oldm){
      oldm = mu;
      sq=sqrt(2.*mu);
      alxm=log(mu);
      g=mu*alxm-gammln(mu+1.);
    }
    do{
      do{
	y=tan(Constants::PI*ran1(seed));
	em=sq*y+mu;
      } while(em < 0.);
      em = floor(em);
      t=0.9*(1.+y*y)*exp(em*alxm-gammln(em+1.)-g);
    } while(ran1(seed) > t);
  }
  return em;
}

/** Generate numbers distributed as Poisson random numbers
 * of mean mu using the rejection method. I have verified
 * that this routine works. Uses ran2.
 *
 * \param mu mean value
 * \param seed seed integer, set negative to initialise.
 */
float Subs::poisson2(float mu, INT4& seed){

  static float sq, alxm,g,oldm=-1.;
  float em,t,y;

  if(mu < 12.){ 
    if(mu != oldm){
      oldm = mu;
      g = exp(-mu);
    }
    em = -1.;
    t = 1.0;
    do {
      ++em;
      t *= ran2(seed);
    } while(t > g);
  }else{
    if(mu != oldm){
      oldm = mu;
      sq=sqrt(2.*mu);
      alxm=log(mu);
      g=mu*alxm-gammln(mu+1.);
    }
    do{
      do{
	y=tan(Constants::PI*ran2(seed));
	em=sq*y+mu;
      } while(em < 0.);
      em = floor(em);
      t=0.9*(1.+y*y)*exp(em*alxm-gammln(em+1.)-g);
    } while(ran2(seed) > t);
  }
  return em;
}

/** Generate numbers distributed as Poisson random numbers
 * of mean mu using the rejection method. I have verified
 * that this routine works. Uses ran3.
 *
 * \param mu mean value
 * \param seed seed integer, set negative to initialise.
 */
float Subs::poisson3(float mu, INT4& seed){

  static float sq, alxm,g,oldm=-1.;
  float em,t,y;

  if(mu < 12.){ 
    if(mu != oldm){
      oldm = mu;
      g = exp(-mu);
    }
    em = -1.;
    t = 1.0;
    do {
      ++em;
      t *= ran3(seed);
    } while(t > g);
  }else{
    if(mu != oldm){
      oldm = mu;
      sq=sqrt(2.*mu);
      alxm=log(mu);
      g=mu*alxm-gammln(mu+1.);
    }
    do{
      do{
	y=tan(Constants::PI*ran3(seed));
	em=sq*y+mu;
      } while(em < 0.);
      em = floor(em);
      t=0.9*(1.+y*y)*exp(em*alxm-gammln(em+1.)-g);
    } while(ran3(seed) > t);
  }
  return em;
}

/** Generate numbers distributed as Poisson random numbers
 * of mean mu using the rejection method. I have verified
 * that this routine works. Uses ran4.
 *
 * \param mu mean value
 * \param seed seed integer, set negative to initialise.
 */
float Subs::poisson4(float mu, INT4& seed){

  static float sq, alxm,g,oldm=-1.;
  float em,t,y;

  if(mu < 12.){ 
    if(mu != oldm){
      oldm = mu;
      g = exp(-mu);
    }
    em = -1.;
    t = 1.0;
    do {
      ++em;
      t *= ran4(seed);
    } while(t > g);
  }else{
    if(mu != oldm){
      oldm = mu;
      sq=sqrt(2.*mu);
      alxm=log(mu);
      g=mu*alxm-gammln(mu+1.);
    }
    do{
      do{
	y=tan(Constants::PI*ran3(seed));
	em=sq*y+mu;
      } while(em < 0.);
      em = floor(em);
      t=0.9*(1.+y*y)*exp(em*alxm-gammln(em+1.)-g);
    } while(ran3(seed) > t);
  }
  return em;
}
