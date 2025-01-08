#include <cstdlib>
#include <iostream>
#include "trm/subs.h"
#include "trm/constants.h"
#include "trm/vec3.h"
#include "trm/roche.h"

/** 
 * set_earth computes the earth vector given an inclination and orbital phase. 
 *
 * \param  iangle orbital inclination
 * \param  phase  orbital phase (1 = one orbit)
 *
 * \return a Subs::Vec3 as the earth vector (unit vector)
 */

Subs::Vec3 Roche::set_earth(double iangle, double phase){

    double cosi, sini, ri = Subs::deg2rad(iangle);
    cosi      = cos(ri);
    sini      = sin(ri);
    
    double sinp, cosp, rp = Constants::TWOPI*phase;
    cosp      = cos(rp);
    sinp      = sin(rp);
    
    return Subs::Vec3(sini*cosp, -sini*sinp, cosi);
}

/** 
 * set_earth computes the earth vector given an inclination and orbital phase. 
 *
 * \param  cosi cosine of orbital inclination
 * \param  sini sine of orbital inclination
 * \param  phase  orbital phase (1 = one orbit)
 *
 * \return a Subs::Vec3 as the earth vector (unit vector)
 */

Subs::Vec3 Roche::set_earth(double cosi, double sini, double phase){

    double sinp, cosp, rp = Constants::TWOPI*phase;
    cosp      = cos(rp);
    sinp      = sin(rp);
    
    return Subs::Vec3(sini*cosp, -sini*sinp, cosi);
}


