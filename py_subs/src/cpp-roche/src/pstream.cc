/*

!!begin
!!title   pstream -- plots impact and Roche lobes
!!author  T.R. Marsh
!!created 19 March 2002
!!descr   plots Roche lobes and stream/disc with impact
!!root    pstream
!!index   pstream
!!class   Programs
!!css     style.css
!!head1   pstream - direct impact stuff

!!emph{pstream} outputs the stream in a double white dwarf binary.


!!table

!!arg{q}{Mass ratiop, = m2/m1}
!!arg{stype}{Stream type, v for velocity, p for position}
!!arg{step}{step size in units of K1+K2 for velocity, or 'a' for position}
!!arg{npoints}{Number of points}
!!table

!!end

*/

#include <cstdlib>
#include <iostream>
#include <string>
#include "trm/subs.h"
#include "trm/input.h"
#include "trm/binary.h"
#include "trm/roche.h"

using Subs::operator+;

int main (int argc, char *argv[]){

  try{

    // Construct Input object

    Subs::Input input(argc, argv, Roche::ROCHE_ENV, Roche::ROCHE_DIR);

    // Define inputs

    input.sign_in("q",       Subs::Input::LOCAL, Subs::Input::PROMPT);
    input.sign_in("stype",   Subs::Input::LOCAL, Subs::Input::PROMPT);
    input.sign_in("step",    Subs::Input::LOCAL, Subs::Input::PROMPT);
    input.sign_in("npoints", Subs::Input::LOCAL, Subs::Input::PROMPT);

    double q;
    input.get_value("q", q, 0.5, 1.e-5, 1.e5, "mass ratio (=m2/m1");
    char stype;
    input.get_value("stype", stype, 'V', "vVpP", "what type of stream? V(elocity) or P(osition)");
    stype = toupper(stype);
    double step;
    input.get_value("step", step, 1.e-3, 1.e-10, 1., "step size between points");
    int npoints;
    input.get_value("npoints", npoints, 500, 2, 10000000, "number of points");

    float *x = new float[npoints], *y = new float[npoints], *r = new float[npoints];

    std::cout << "# " << std::endl;
    std::cout << "# q     = " << q     << std::endl;
    std::cout << "# stype = " << stype << std::endl;
    std::cout << "# " << std::endl;
    if(stype == 'V'){
      Roche::vstream(q, step, x, y, r, npoints, 1);
      for(int i=0; i<npoints; i++)
	std::cout << x[i] << " " << y[i] << " " << r[i] << "\n";
    }else{
      Roche::stream(q, step, x, y, npoints);
      for(int i=0; i<npoints; i++)
	std::cout << x[i] << " " << y[i] << "\n";
    }
    delete[] r;
    delete[] y;
    delete[] x;
  }
  catch(const std::string& err){
    std::cerr << err << std::endl;
  }
}









