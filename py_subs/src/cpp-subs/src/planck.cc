#include "trm/subs.h"
#include "trm/constants.h"

/** Computes the Planck function Bnu = (2 h \nu^3/c^2)/(exp(h \nu/kT) - 1)
 *  as a function of wavelength and temperature. Output units are W/m**2/Hz/sr.
 *
 * \param wave wavelength in nanometres
 * \param temp temperature in K
 */

double Subs::planck(double wave, double temp){
  
    const double FAC1 = 2.e27*Constants::H*Constants::C;
    const double FAC2 = 1.e9*Constants::H*Constants::C/Constants::K;
    
    double efac = FAC2/(wave*temp);
    if(efac > 40.){
	return FAC1*exp(-efac)/(wave*sqr(wave));
    }else{
	return FAC1/(exp(efac)-1.)/(wave*sqr(wave));
    }
}

/** Computes the logarithmic derivative of the Planck function 
 * Bnu wrt wavelength (i.e. d ln(Bnu) / d ln(lambda)) as a function of wavelength and temperature
 * \param wave wavelength in nanometres
 * \param temp temperature in K
 */

double Subs::dplanck(double wave, double temp){
    
    const double FAC2 = 1.e9*Constants::H*Constants::C/Constants::K;
    
    double efac = FAC2/(wave*temp);
    return efac/(1.-exp(-efac)) - 3.;
}

/** Computes the logarithmic derivative of the Planck function 
 * Bnu wrt T (i.e. d ln(Bnu) / d ln(T)) as a function of wavelength and temperature
 * \param wave wavelength in nanometres
 * \param temp temperature in K
 */

double Subs::dlpdlt(double wave, double temp){
    
    const double FAC2 = 1.e9*Constants::H*Constants::C/Constants::K;
    
    double efac = FAC2/(wave*temp);
    return efac/(1.-exp(-efac));
}

