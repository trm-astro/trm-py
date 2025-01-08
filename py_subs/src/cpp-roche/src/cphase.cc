/*

!!begin
!!title   cphase -- contact phjases
!!author  T.R. Marsh
!!created 20 Aug 2002
!!descr   computes contact phases in Roche geometry
!!root    cphase
!!index   cphase
!!class   Programs
!!css     style.css
!!head1   cphase -- contact phase computation

This program computes the phases of third and fourth contacts of a white dwarf in a CV.

!!table
!!arg{q}{Mass ratio = M1/M2}
!!arg{i}{Inclination, degrees}
!!arg{m1}{Mass or primary, solar masses}
!!arg{period}{Orbital period, days}
!!arg{ntheta}{Number of angles around the star being eclipsed}
!!arg{acc}{Accuracy for eclipse computation}
!!table

!!end

*/

#include <cstdlib>
#include <iostream>
#include "trm/subs.h"
#include "trm/format.h"
#include "trm/input.h"
#include "trm/binary.h"
#include "trm/roche.h"

using Subs::operator+;

int main (int argc, char *argv[]){

  try{

    Subs::Format form(8);

    // Construct Input object

    Subs::Input input(argc, argv, Roche::ROCHE_ENV, Roche::ROCHE_DIR);

    // Define inputs

    input.sign_in("q",       Subs::Input::LOCAL, Subs::Input::PROMPT);
    input.sign_in("i",       Subs::Input::LOCAL, Subs::Input::PROMPT);
    input.sign_in("m1",      Subs::Input::LOCAL, Subs::Input::PROMPT);
    input.sign_in("period",  Subs::Input::LOCAL, Subs::Input::PROMPT);
    input.sign_in("ntheta",  Subs::Input::LOCAL, Subs::Input::PROMPT);
    input.sign_in("acc",     Subs::Input::LOCAL, Subs::Input::PROMPT);

    double q;
    input.get_value("q",  q, 0.2, 1.e-5, 1.e5, "mass ratio = M2/M1");
    double iangle;
    input.get_value("i",  iangle, 88., 1., 90., "orbital inclination (degrees)");
    double m1;
    input.get_value("m1",  m1, 0.9, 0.01, 1.38, "mass of white dwarf (solar)");
    double m2 = q*m1;
    double r1 = Binary::mr_wd_eggleton(m1);
    std::cout << "Radius of white dwarf = " << form(r1) << " solar radii" << std::endl;
    double period;
    input.get_value("period", period, 0.15, 0.01, 10., "orbital period (days)");
    double separation = Binary::orbital_separation(m1, m2, 86400*period);
    r1 /= separation;
    std::cout << "Radius of white dwarf/separation = " << form(r1) << std::endl;
    int ntheta;
    input.get_value("ntheta",  ntheta, 100, 2, 100000, "number of angles around star being eclipsed");
    double acc;
    input.get_value("acc",  acc, 0.005, 0., 0.2, "accuracy for eclipse computation by 'blink'");

    double third;
    bool yes_third = Roche::third_contact(q, iangle, r1, ntheta, acc, third);
    if(yes_third){
      std::cout << " Third contact phase  = " << form(third) << std::endl;
    }else{
      std::cout << "There is no third contact" << std::endl;
    }

    double fourth;
    bool yes_fourth = Roche::fourth_contact(q, iangle, r1, ntheta, acc, fourth);
    if(yes_fourth){
      std::cout << "Fourth contact phase  = " << form(fourth) << std::endl;
    }else{
      std::cout << "There is no fourth contact" << std::endl;
    }

    if(yes_third && yes_fourth)
      std::cout << "Phase width of white dwarf eclipse (mid-ingress to mid-egress) = " << form(third+fourth)
	   << ", phase width of ingress and egress = " << form(fourth-third) << std::endl;
      std::cout << "Width of white dwarf eclipse (mid-ingress to mid-egress) = " << form(86400*period*(third+fourth))
	   << " seconds, width of ingress and egress = " << form(86400*period*(fourth-third)) << " seconds" << std::endl;
      std::cout << "Width of white dwarf eclipse (mid-ingress to mid-egress) = " << form(period*(third+fourth))
	   << " days, width of ingress and egress = " << form(period*(fourth-third)) << " days" << std::endl;
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









