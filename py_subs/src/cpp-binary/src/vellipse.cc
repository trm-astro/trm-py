#include <cstdlib>
#include <iostream>
#include "trm/constants.h"
#include "trm/subs.h"
#include "trm/binary.h"

/** vellipse compute the x and y velocity of a particle executing
 * an elliptical orbit about a star located at the origin. The long axis of
 * the ellipse is parallel to the x axis and it extends in the negative x direction, i.e.
 * the periastron is on the negative x-axis. The particle orbits anti-clockwise.
 * \param m mass of star (solar)
 * \param e eccentricity of orbit.
 * \param a semi-major axis (solar radii)
 * \param angle from periastron, radians
 * \return Returns a structure with all the information.
 */ 

Binary::Einfo Binary::vellipse(double m, double e, double a, double theta){

    // Compute polar coords
    double sint  = sin(theta);
    double cost  = cos(theta);

    // semi latus rectum
    double l = a*(1-e*e);

    // useful number
    double fac = 1 + e*cost;

    // orbital frequency
    double omega = sqrt(Constants::G*Constants::MSUN*m/pow(Constants::RSUN*a,3));

    // eccentric anomaly
    double cosE = (e + cost)/fac;
    double E    = acos(cosE);
    E = (sint >= 0) ? E : Constants::TWOPI - E;
    double sinE = sin(E);

    // rate of change of eccentric anom.
    double Edot = omega/(1-e*cosE);

    // rate of change of theta
    double zeta = (1+e)/(1-e);
    double cosE2 = cos(E/2);
    double sinE2 = sin(E/2);
    double fac0  = Subs::sqr(cosE2) + zeta*Subs::sqr(sinE2);
    double thetadot = sqrt(zeta)*Edot/fac0;

    // velocity in radial direction
    double vr = e*l*sint*thetadot/(fac*fac);

    // theta velocity
    double vt = l*thetadot/fac;

    // x, y velocity in solar radii/sec
    double vx = vr*cost-vt*sint;
    double vy = vr*sint + vt*cost;

    // Now for the partial derivatives. _x means
    // wrt x etc.

    double r       = l/fac;
    double theta_x = -sint/r;
    double theta_y = +cost/r;
    double r_x     = cost;
    double r_y     = sint;
    double l_x     = r_x*fac - e*r*sint*theta_x;
    double l_y     = r_y*fac - e*r*sint*theta_y;
    double omega_x = -3*omega/2.*l_x/l;
    double omega_y = -3*omega/2.*l_y/l;
    double E_x, E_y;
    if(sint == 0){
        E_x = 0.;
        E_y = 0.;
    }else{
        E_x = (1-e*e)/Subs::sqr(fac)*sint/sinE*theta_x;
        E_y = (1-e*e)/Subs::sqr(fac)*sint/sinE*theta_y;
    }

    double fac1   = 1 - e*cosE;
    double Edot_x = (omega_x - e*omega*sinE*E_x/fac1)/fac1;
    double Edot_y = (omega_y - e*omega*sinE*E_y/fac1)/fac1;

    double thetadot_x = sqrt(zeta)/fac0*(Edot_x-(zeta-1)*sinE*E_x*Edot/2./fac0);
    double thetadot_y = sqrt(zeta)/fac0*(Edot_y-(zeta-1)*sinE*E_y*Edot/2./fac0);

    double vr_x = e/(fac*fac)*(l*sint*thetadot_x + l*cost*theta_x*thetadot + l_x*sint*thetadot + 2*e/fac*l*sint*sint*theta_x*thetadot);
    double vr_y = e/(fac*fac)*(l*sint*thetadot_y + l*cost*theta_y*thetadot + l_y*sint*thetadot + 2*e/fac*l*sint*sint*theta_y*thetadot);
    double vt_x = 1/fac*(l*thetadot_x + l_x*thetadot + e/fac*l*sint*theta_x*thetadot);
    double vt_y = 1/fac*(l*thetadot_y + l_y*thetadot + e/fac*l*sint*theta_y*thetadot);

    double vx_x = vr_x*cost - vt_x*sint - (vr*sint + vt*cost)*theta_x;
    double vx_y = vr_y*cost - vt_y*sint - (vr*sint + vt*cost)*theta_y;
    double vy_x = vr_x*sint + vt_x*cost + (vr*cost - vt*sint)*theta_x;
    double vy_y = vr_y*sint + vt_y*cost + (vr*cost - vt*sint)*theta_y;

    // Return data, converting velocities to km/s
    Einfo temp = {1.e-3*Constants::RSUN*vx, 1.e-3*Constants::RSUN*vy, vx_x, vx_y, vy_x, vy_y, r, omega, E};
    return temp;
}


