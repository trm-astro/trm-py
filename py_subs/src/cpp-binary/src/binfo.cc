/*

!!begin
!!title   binfo -- gives lots of stuff about binaries
!!author  T.R. Marsh
!!created 27 Sep 2006
!!descr   gives lots of stuff about binaries
!!root    binfo
!!index   binfo
!!class   Programs
!!css     style.css
!!head1   binfo -- gives lots of stuff about binaries

Given masses and a period, this reports the separation and gravitational
wave angular momentum loss rate.

!!table
!!arg{m1}{Mass of star 1}
!!arg{m2}{Mass of star 2}
!!arg{period}{Orbital period, days}
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

    input.sign_in("m1",       Subs::Input::LOCAL, Subs::Input::PROMPT);
    input.sign_in("m2",       Subs::Input::LOCAL, Subs::Input::PROMPT);
    input.sign_in("period",  Subs::Input::LOCAL, Subs::Input::PROMPT);

    double m1;
    input.get_value("m1", m1, 0.6, 1.e-5, 1.e5, "Mass of primary");
    double m2;
    input.get_value("m2", m2, 0.6, 1.e-5, 1.e5, "Mass of secondary");
    double period;
    input.get_value("period", period, 0.1, 1.e-3, 1.e3, "orbital period (days)");

    period *= Constants::DAY;

    double a = Binary::orbital_separation(m1, m2, period);
    double jdotgr = Binary::jdotgr(m1, m2, a);

    Subs::Format form(8);
    std::cout << "Orbital separation = " << form(a) << " solar radii" << std::endl;
    std::cout << "Jdot/J(GR)         = " << form(jdotgr) << "/second" << std::endl;

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









