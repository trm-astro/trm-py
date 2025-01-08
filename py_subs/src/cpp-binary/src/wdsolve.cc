/*

!!begin
!!title   wdsolve -- solves for parameters of a white-dwarf binary
!!author  T.R. Marsh
!!created 19 May 2006
!!descr   solves for the parameters of a white-dwarf binary
!!root    wdsolve
!!index   wdsolve
!!class   Programs
!!css     style.css
!!head1   wdsolve - solves for the parameters of a white dwarf binary

Given a mass ratio, inclination, orbital period and radius of white dwarf 
scaled by separation (presumed to be the primary star), this program computes 
the masses of the two stars and the radius of the white dwarf assuming a 
mass-radius relation for the white dwarf. This is the sort of information 
one has for an eclipsing white dwarf in a binary.

!!table
!!arg{q}{Mass ratio = M2/M1}
!!arg{iangle}{Orbital inclination angle, degrees}
!!arg{period}{Orbital period, days}
!!arg{r1}{Radius of the white dwarf, units of separation}
!!arg{over}{Oversize factor. The factor by which the accretor exceeds the zero-temperature radius from Eggleton. 
To estimate effect of temperature expansion. over=1 gives standard zero temp result.}
!!table

!!end

*/

#include <cstdlib>
#include <iostream>
#include "trm/subs.h"
#include "trm/format.h"
#include "trm/input.h"
#include "trm/binary.h"

int main (int argc, char *argv[]){

  try{

    // Construct Input object

    Subs::Input input(argc, argv, Binary::BINARY_ENV, Binary::BINARY_DIR);

    // Define inputs

    input.sign_in("q",       Subs::Input::LOCAL, Subs::Input::PROMPT);
    input.sign_in("iangle",  Subs::Input::LOCAL, Subs::Input::PROMPT);
    input.sign_in("period",  Subs::Input::LOCAL, Subs::Input::PROMPT);
    input.sign_in("r1",      Subs::Input::LOCAL, Subs::Input::PROMPT);
    input.sign_in("over",    Subs::Input::LOCAL, Subs::Input::PROMPT);

    double q;
    input.get_value("q",  q, 0.5, 1.e-5, 1.e5, "mass ratio = M2/M1");
    double iangle;
    input.get_value("iangle", iangle, 80., 0., 90., "orbital inclination angle (degrees)");
    double period;
    input.get_value("period", period, 0.1, 1.e-3, 1.e3, "orbital period (days)");
    double r1;
    input.get_value("r1", r1, 0.01, 1.e-5, 0.99, "radius of white dwarf (units of separation)");
    double over;
    input.get_value("over", over, 1., 0.1, 10., "oversize factor relative to zero temp radius");


    period *= Constants::DAY;

    double mlo = 0.05, mhi = 1.39;
    if(over*Binary::mr_wd_eggleton(mlo) < r1*Binary::orbital_separation(mlo, mlo*q, period) ||
       over*Binary::mr_wd_eggleton(mhi) > r1*Binary::orbital_separation(mhi, mhi*q, period))
      throw Binary::Binary_Error("Cannot find solution in range " + Subs::str(mlo) + " to " + Subs::str(mhi));

    double m = (mlo+mhi)/2.;
    while(abs(mhi-mlo) > 1.e-7){

      if(over*Binary::mr_wd_eggleton(m) < r1*Binary::orbital_separation(m, m*q, period)) 
	mhi = m;
      else
	mlo = m;

      m = (mlo+mhi)/2.;

    }
    double a = Binary::orbital_separation(m, m*q, period);
    Subs::Format form(8);
    std::cout << std::endl;
    std::cout << "a  = " << form(a) << " solar radii"    << std::endl;
    std::cout << "M1 = " << form(m) << " solar masses"   << std::endl;
    std::cout << "M2 = " << form(q*m) << " solar masses" << std::endl;
    double vscale = Constants::TWOPI/period*Constants::RSUN/1e3*a*sin(Constants::TWOPI*iangle/360.);
    std::cout << "K1 = " << form(vscale*q/(1+q)) << " km/s" << std::endl;
    std::cout << "K2 = " << form(vscale/(1+q))   << " km/s" << std::endl;
    std::cout << "log(g) of white dwarf = " << form(log10(Constants::G*Constants::MSUN*m/pow(Constants::RSUN*r1*a,2))+2.) << std::endl;
    double jdotgr = Binary::jdotgr(m, q*m, a);
    double b = 3*jdotgr/2*period*period;
    std::cout << "Quadratic coeff of ephemeris = " << form(b)   << " sec" << std::endl;

  }

  catch(const Binary::Binary_Error& err){
    std::cerr << "Binary::Binary_Error exception:" << std::endl;
    std::cerr << err << std::endl;
    exit(EXIT_FAILURE);
  }
  catch(const std::string& err){
    std::cerr << "string exception:" << std::endl;
    std::cerr << err << std::endl;
    exit(EXIT_FAILURE);
  }
}









