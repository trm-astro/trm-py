#include <cfloat>
#include <cstdlib>
#include <iostream>
#include "trm/subs.h"
#include "trm/format.h"
#include "trm/constants.h"
#include "trm/input.h"
#include "trm/vec3.h"
#include "trm/binary.h"
#include "trm/roche.h"

/*

!!begin
!!title impact -- computes impact stuff
!!author T.R. Marsh
!!revised 24 April 2006
!!descr computes impact stuff to do with V407 Vul
!!root  impact
!!index impact
!!class Programs
!!css   style.css
!!head1 impact - direct impact computation.

This program computes whether a direct impact between two white dwarfs
occurs and where it does so if it does. It reports X, Y and the radius scaled by 
the binary separation, as well as the L1 distance plus the separation in solar radii
and orbital period in seconds.

!!head2 Arguments

!!table
!!arg{m1}{Mass of accretor, solar masses. Its radius will be computed using Eggleton's M-R relation.}
!!arg{m2}{Mass of donor, solar masses.}
!!arg{period}{Orbital period, seconds. If <= 0, it will be computed assuming that the secondary also 
follows Eggleton's M-R relation}
!!table

!!end

*/

int main (int argc, char *argv[]){

  try{

    // Construct Input object

    Subs::Input input(argc, argv, Roche::ROCHE_ENV, Roche::ROCHE_DIR);

    // Define inputs

    input.sign_in("m1",      Subs::Input::LOCAL, Subs::Input::PROMPT);
    input.sign_in("m2",      Subs::Input::LOCAL, Subs::Input::PROMPT);
    input.sign_in("period",  Subs::Input::LOCAL, Subs::Input::PROMPT);

    double m1;
    input.get_value("m1",  m1, 0.6, 1.e-3, 1.4, "primary mass (solar)");
    double m2;
    input.get_value("m2", m2, 0.2, 1.e-4, 1.4, "secondary mass (solar)");
    double period;
    input.get_value("period", period, 500., -1e10, 1e5, "orbital period (seconds, <=0 automatic)");

    const double q    = m2/m1;
    double r1         = Binary::mr_wd_eggleton(m1);
    double r2, a, rl2 = Roche::rlobe_eggleton(q);

    if(period <= 0){
      r2  = Binary::mr_wd_eggleton(m2);
      a   = r2/rl2;
      period   = Binary::orbital_period(m1, m2, a);
    }else{
      a = Binary::orbital_separation(m1, m2, period);
    }

    Subs::Vec3 p;
    p.x()   = r1/a/sqrt(2.);
    p.y()   = r1/a/sqrt(2.);
    p.z()   = 0.;
    double pot   = Roche::rpot(q,p), x, y;
    Subs::Format form(8);

    if(Roche::hits(q,pot,x,y)){
      std::cout << "Stream hits at X = " << form(x) << ", Y = " << form(y) << " (units of a)" << std::endl;
      std::cout << "Radius = " << form(sqrt(x*x+y*y)) << ", RL1 = " << form(Roche::xl1(q)) << " (units of a)" <<std::endl;
      std::cout << "Angle  = " << form(360.*atan2(y,x)/Constants::TWOPI) << " degrees ahead of secondary star" <<std::endl;
      
    }else{
      std::cout << "Stream does not hit. RL1 = " << Roche::xl1(q) << std::endl;
    }
    std::cout << "Separation = " << form(a) << " solar radii, period = " << form(period) << " seconds." << std::endl;
  }

  catch(const Roche::Roche_Error& err){
    std::cerr << "Roche::Roche_Error exception:" << std::endl;
    std::cerr << err.what() << std::endl;
    exit(EXIT_FAILURE);
  }
  catch(const std::string& err){
    std::cerr << "string exception:" << std::endl;
    std::cerr << err << std::endl;
    exit(EXIT_FAILURE);
  }
}









