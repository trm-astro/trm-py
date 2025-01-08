/*

!!begin
!!title   vspace -- velocity space plot
!!author  T.R. Marsh
!!created Mar 2004
!!descr   velocity space plot
!!root    vspace
!!index   vspace
!!class   Programs
!!css     style.css
!!head1   vspace - velocity space plot

This program plots the Roche lobes and lines of equal radial velocities
over it in velocity space.

!!table
!!arg{q}{Mass ratio = m2/m1}
!!arg{vscale}{Velocity scale = k1+k2}
!!arg{vstep}{Velocity step between equal RV lines}
!!arg{nstep}{Number of steps both to negative and positive velocities}
!!arg{r1}{Radius of the accretor (units of Rl1)}
!!arg{rdisc}{Disc radius (units of Rl1)}
!!arg{phase}{Orbital phase}
!!arg{device}{Plot device}
!!arg{width}{Plot width, inches}
!!arg{x1}{Left X limit for plot, velocity units}
!!arg{x2}{Right X limit for plot, velocity units}
!!arg{y1}{Lower Y limit for plot, velocity units}
!!arg{y2}{Upper Y limit for plot, velocity units}
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
    input.get_value("width", width, 0.f, 0.f, 100.f, "left X limit for plot");
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

    Subs::Plot plot(device);
    cpgpap(width,1.);
    cpgslw(lwidth);
    cpgsvp(0.,1.,0.,1.);
    cpgwnad(x1, x2, y1, y2);
    const int NLOBE = 200;
    const int NSTRM = 200;
    const int NPLOT = std::max(NLOBE, NSTRM);
    Subs::Buffer1D<float> vx(NPLOT), vy(NPLOT), rad(NSTRM);

    Roche::vlobe1(q, vx, vy, NPLOT);
    for(int i=0; i<NPLOT; i++){
      vx[i] *= vscale;
      vy[i] *= vscale;
    }
    cpgsls(2);
    cpgline(NPLOT, vx, vy);

    Roche::vlobe2(q, vx, vy, NPLOT);
    for(int i=0; i<NPLOT; i++){
      vx[i] *= vscale;
      vy[i] *= vscale;
    }
    cpgsls(1);
    cpgsfs(3);
    if(colours)
      cpgsci(2);
    cpgpoly(NPLOT, vx, vy);
    cpgline(NPLOT, vx, vy);

    float kw    = vscale*q/(1.+q);
    float xl1   = Roche::xl1(q);
    float vdisc = vscale/sqrt((1+q)*xl1*rdisc);

    cpgsfs(2);
    cpgsls(1);
    if(colours)
      cpgsci(1);
    cpgcirc(0.f, -kw, vdisc);

    float vaccrete = vscale/sqrt((1+q)*xl1*r1);
    cpgcirc(0.f, -kw, vaccrete);

    Roche::vstream(q, 0.01, vx, vy, rad, NSTRM, 1);
    int nlast = 0;
    for(int i=0; i<NSTRM; i++){
      vx[i] *= vscale;
      vy[i] *= vscale;
      if(rad[i] > xl1*rdisc) nlast++;
    }

    if(colours)
      cpgsci(2);
    cpgline(nlast, vx, vy);    
    if(colours)
      cpgsci(3);
    cpgsch(3);
    cpgpt1( vx[nlast-1], vy[nlast-1], 18);

    Roche::vstream(q, 0.01, vx, vy, rad, NSTRM, 2);
    nlast = 0;
    for(int i=0; i<NSTRM; i++){
      vx[i] *= vscale;
      vy[i] *= vscale;
      if(rad[i] > xl1*rdisc) nlast++;
    }

    if(colours)
      cpgsci(2);
    cpgline(nlast, vx, vy);    
    if(colours)
      cpgsci(3);
    cpgsch(3);
    cpgpt1( vx[nlast-1], vy[nlast-1], 18);

    if(colours)
      cpgsci(1);

    float cosp = cos(Constants::TWOPI*phase);
    float sinp = sin(Constants::TWOPI*phase);
    float xp1, yp1, xp2, yp2;
    bool ok_to_plot;
    for(int i=-nstep; i<=nstep; i++){
      float vrad = vstep*i;
      if(i < 0){
	if(colours)
	  cpgsci(2);
	else
	  cpgsls(1);
      }else if(i > 0){
	if(colours)
	  cpgsci(4);
	else
	  cpgsls(2);
      }else{
	cpgsci(1);
      }
      if(cosp == 0.f){
	xp1 = x1;
	yp1 = -vrad;
	xp2 = x2;
	yp2 = -vrad;
	ok_to_plot = (yp1 > y1 && yp1 < y2);

      }else if(sinp == 0.f){
	yp1 = y1;
	xp1 = vrad;
	yp2 = y2;
	xp2 = vrad;
	ok_to_plot = (xp1 > x1 && xp1 < x2);

      }else{
	xp1 = x1;
	yp1 = (xp1*cosp - vrad)/sinp;
	xp2 = x2;
	yp2 = (xp2*cosp - vrad)/sinp;
	if(yp1 < y1 && yp2 > y1){
	  yp1 = y1;
	  xp1 = (yp1*sinp + vrad)/cosp;
	}
	if(yp1 > y2 && yp2 < y2){
	  yp1 = y2;
	  xp1 = (yp1*sinp + vrad)/cosp;
	}
	if(yp2 < y1 && yp1 > y1){
	  yp2 = y1;
	  xp2 = (yp2*sinp + vrad)/cosp;
	}
	if(yp2 > y2 && yp1 < y2){
	  yp2 = y2;
	  xp2 = (yp2*sinp + vrad)/cosp;
	}
	ok_to_plot = (yp1 >= y1 && yp1 <= y2 && yp2 >= y1 && yp2 <= y2);
      }

      if(ok_to_plot){
	cpgmove(xp1, yp1);
	cpgdraw(xp2, yp2);
      }
    }

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









