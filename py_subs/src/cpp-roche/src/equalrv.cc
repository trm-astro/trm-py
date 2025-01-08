/*

!!begin
!!title   equalrv -- shows radial velocity lines
!!author  T.R. Marsh
!!created Mar 2004
!!descr   shows radial velocity pattern on disc in binary
!!root    equalrv
!!index   equalrv
!!class   Programs
!!css     style.css
!!head1   equalrv - direct impact computation.

This program plots the Roche lobes and lines of equal radial velocities
over it.

!!table
!!arg{q}{Mass ratio = m2/m1}
!!arg{vscale}{Velocity scale = k1+k2}
!!arg{vstep}{Velocity step between equal RV lines}
!!arg{nstep}{Number of steps both to negative and positive velocities}
!!arg{r1}{Radius of the accretor (units of Rl1)}
!!arg{rdisc}{Disc radius (units of Rl1)}
!!arg{phase}{Orbital phase}
!!arg{device}{Plot device}
!!arg{width}{Plot width, inches, 0 for default}
!!arg{x1}{Left X limit for plot}
!!arg{x2}{Right X limit for plot}
!!arg{y1}{Lower Y limit for plot}
!!arg{y2}{Upper Y limit for plot}
!!arg{lwidth}{Line width, multiple of default}
!!arg{colours}{Do you want colours?}
!!table

!!end

*/

#include <cstdlib>
#include <cfloat>
#include <iostream>
#include "trm/subs.h"
#include "trm/input.h"
#include "trm/plot.h"
#include "trm/roche.h"

using Subs::operator+;

int main (int argc, char *argv[]){

  try{

    // Construct Input object
    Subs::Input input(argc, argv, Roche::ROCHE_ENV, Roche::ROCHE_DIR);

    // Define inputs
    input.sign_in("q",      Subs::Input::LOCAL, Subs::Input::PROMPT);
    input.sign_in("vscale", Subs::Input::LOCAL, Subs::Input::PROMPT);
    input.sign_in("vstep",  Subs::Input::LOCAL, Subs::Input::PROMPT);
    input.sign_in("nstep",  Subs::Input::LOCAL, Subs::Input::PROMPT);
    input.sign_in("r1",     Subs::Input::LOCAL, Subs::Input::PROMPT);
    input.sign_in("rdisc",  Subs::Input::LOCAL, Subs::Input::PROMPT);
    input.sign_in("phase",  Subs::Input::LOCAL, Subs::Input::PROMPT);
    input.sign_in("device", Subs::Input::GLOBAL, Subs::Input::PROMPT);
    input.sign_in("width",  Subs::Input::GLOBAL, Subs::Input::PROMPT);
    input.sign_in("x1",     Subs::Input::LOCAL, Subs::Input::PROMPT);
    input.sign_in("x2",     Subs::Input::LOCAL, Subs::Input::PROMPT);
    input.sign_in("y1",     Subs::Input::LOCAL, Subs::Input::PROMPT);
    input.sign_in("y2",     Subs::Input::LOCAL, Subs::Input::PROMPT);
    input.sign_in("lwidth", Subs::Input::LOCAL, Subs::Input::PROMPT);
    input.sign_in("colours",Subs::Input::LOCAL, Subs::Input::PROMPT);

    double q;
    input.get_value("q", q, 0.2, 1.e-3, 1000., "mass ratio (m1/m2)");
    double vscale;
    input.get_value("vscale", vscale, 500., 0.1, 10000., "velocity scale (k1+k2)");
    double vstep;
    input.get_value("vstep", vstep, 100., 0.001, 5000., "velocity step");
    int nstep;
    input.get_value("nstep", nstep, 20, 0, 100, "number of steps");
    double r1;
    input.get_value("r1", r1, 0.02, 0., 1., "accretor radius (units of Rl1)");
    double rdisc;
    input.get_value("rdisc", rdisc, 0.7, 0.001, 1., "disc radius (units of Rl1)");
    double phase;
    input.get_value("phase", phase, 0., -1., 1., "orbital phase");
    std::string device;
    input.get_value("device", device, "/xs", "plot device");
    float width;
    input.get_value("width", width, 0.f, 0.f, 100.f, "width of plot in inches");
    float x1;
    input.get_value("x1", x1, -1.f, -FLT_MAX, FLT_MAX, "left X limit for plot");
    float x2;
    input.get_value("x2", x2,  1.f, -FLT_MAX, FLT_MAX, "right X limit for plot");
    float y1;
    input.get_value("y1", y1, -1.f, -FLT_MAX, FLT_MAX, "lower X limit for plot");
    float y2;
    input.get_value("y2", y2,  1.f, -FLT_MAX, FLT_MAX, "upper X limit for plot");
    int lwidth;
    input.get_value("lwidth", lwidth,  1, 1, 100, "line width");
    bool colours;
    input.get_value("colours", colours,  true, "do you want colours?");

    float r11  = cos(Constants::TWOPI*phase);
    float r12  = -sin(Constants::TWOPI*phase);
    float r21  = -r12;
    float r22  = r11;
    float cofm = q/(1.+q);

    Subs::Plot plot(device);
    cpgpap(width,1.);
    cpgslw(lwidth);
    cpgsvp(0.,1.,0.,1.);
    cpgwnad(x1, x2, y1, y2);
    const int NLOBE = 200;
    const int NRADV = 50;
    const int NSTRM = 50;
    const int NPLOT = std::max(NLOBE, std::max(NRADV, NSTRM));
    Subs::Buffer1D<float> x(NPLOT), y(NPLOT);

    Roche::lobe1(q, x, y, NPLOT);
    float xt, yt;
    for(int i=0; i<NPLOT; i++){
      xt   = r11*(x[i]-cofm) + r12*y[i];
      yt   = r21*(x[i]-cofm) + r22*y[i];
      x[i] = xt;
      y[i] = yt;
    }
    cpgsls(2);
    cpgline(NPLOT, x, y);

    if(colours)
      cpgsci(2);
    Roche::lobe2(q, x, y, NPLOT);
    for(int i=0; i<NPLOT; i++){
      xt   = r11*(x[i]-cofm) + r12*y[i];
      yt   = r21*(x[i]-cofm) + r22*y[i];
      x[i] = xt;
      y[i] = yt;
    }
    cpgsls(1);
    cpgsfs(3);
    cpgpoly(NPLOT, x, y);
    cpgline(NPLOT, x, y);
    float xa = -cofm*r11;
    float ya = -cofm*r21;
    float xl1 = Roche::xl1(q);
    rdisc *= xl1;
    cpgsfs(1);
    cpgsls(1);
    
    for(int i=-nstep; i<=nstep; i++){
      float vrad = vstep*i;
      float c = sqrt(1+q)*(vrad/vscale+q*r21/(1+q));
      float rmax = std::min(rdisc, 0.9999999/(c*c));
      for(int j=0; j<NRADV; j++){
	float r = rmax*j/(NRADV-1);
	float sint = c*sqrt(r);
	float cost = sqrt(1.-sint*sint);
	x[j] = xa + r*cost;
	y[j] = ya + r*sint;
      }
      if(i < 0){
	if(colours)
	  cpgsci(4);
	else
	  cpgsls(2);
      }else if(i > 0){
	if(colours)
	  cpgsci(2);
	else
	  cpgsls(1);
      }else{
	if(colours)
	  cpgsci(1);
      }
      cpgline(NRADV, x, y);
      for(int j=0; j<NRADV; j++) x[j] = xa-(x[j]-xa);
      cpgline(NRADV, x, y);
    }
    if(colours)
      cpgsci(1);
    cpgsfs(1);
    cpgcirc(xa, ya, r1);
    cpgsfs(2);
    cpgsls(1);
    cpgcirc(xa, ya, rdisc);

    cpgsls(1);
    if(colours)
      cpgsci(2);
    Roche::streamr(q, 0.9*rdisc, x, y, NSTRM);
    for(int i=0; i<NPLOT; i++){
      xt   = r11*(x[i]-cofm) + r12*y[i];
      yt   = r21*(x[i]-cofm) + r22*y[i];
      x[i] = xt;
      y[i] = yt;
    }
    cpgline(NSTRM, x, y);
    if(colours)
      cpgsci(3);
    cpgsch(3);
    cpgpt1(x[NSTRM-1],y[NSTRM-1],18);
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









