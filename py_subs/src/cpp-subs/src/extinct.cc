#include "trm/subs.h"

/** Returns extinction in magnitudes as a function of wavelength scaled
 * to the extinction in magnitudes at V according to Cardelli, Clayton &
 * Mathis, 1989, ApJ, 345, 245.
 *
 * \param lambda wavelength, microns. The function returns something for 
 * lambda in the range 0.1 to 3.5 microns but remember that there is evidence for
 * deviations from the mean law used here, especially below 0.16 microns.  The function will return
 * the value of the closest limit if specified outside this range. Refer to the above paper 
 * if you are really concerned about accuracy.
 * \param ratio this is the usual A(V)/E(B-V) ratio. If in doubt use 3.1, but note
 * that the shape of the extinction does depend upon this parameter.
 */

double Subs::extinct(double lambda, double ratio){
  
  double x = 1./lambda;
  double a = 0, b = 0;

  if(x <= 1.1){

    // Infrared
    a =  0.574*pow(std::max(x,0.3), 1.61);
    b = -0.527*pow(std::max(x,0.3), 1.61);

  }else if(x > 1.1 && x <= 3.3){

    // Optical
    double y = x - 1.82;
    a = 1  + y*(0.17699 + y*(-0.50447 + y*(-0.02427 + y*(0.72085 + y*(0.01979 + y*(-0.77530 + y*0.32999))))));
    b = y*(1.41338 + y*(2.28305 + y*(1.07233 + y*(-5.38434  + y*(-0.62251 + y*(5.30260 - y*2.09002))))));

  }else if(x > 3.3 && x <= 8){

    // Ultraviolet
    a =  1.752 - 0.316*x - 0.104/(sqr(x-4.67) + 0.341);
    b = -3.090 + 1.825*x + 1.206/(sqr(x-4.62) + 0.263);
    if(x > 5.9){
      a +=  -0.04473*sqr(x-5.9) - 0.009779*pow(x-5.9, 3);
      b +=  +0.2130*sqr(x-5.9)  + 0.1207*pow(x-5.9, 3);
    }

  }else if(x > 8){

    // Far ultraviolet
    double y = std::min(x-8.,2.);

    a = -1.073 + y*(-0.628 + y*(0.137 - 0.070*y)); 
    b = 13.670 + y*(4.257 + y*(-0.420 + 0.374*y));

  }

  return a + b/ratio;
}
