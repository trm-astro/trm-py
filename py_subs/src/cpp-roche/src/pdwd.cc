/*

!!begin
!!title   pdwd -- plots impact and Roche lobes
!!author  T.R. Marsh
!!created 19 March 2002
!!revised 31 March 2005
!!descr   plots Roche lobes and stream/disc with impact
!!root    pdwd
!!index   pdwd
!!class   Programs
!!css     style.css
!!head1   pdwd - direct impact stuff

!!emph{pdwd} plots the stream in a double white dwarf binary.


!!table

!!arg{m1}{Mass of accretor, solar masses}

!!arg{m2}{Mass of donor, solar masses}

!!arg{type}{Type of system: 'cv' ==> a low mass MS donor with
R = M in solar units. 'am' ==> white dwarf donor.}

!!arg{device}{Plot device}

!!arg{x1}{lower x limit of plot (units of a)}

!!arg{x2}{upper x limit of plot}

!!arg{y1}{lower y limit of plot}

!!arg{y2}{upper y limit of plot}

!!arg{width}{width of plot inches}

!!arg{aspect}{aspect ratio}

!!arg{lwidth}{width of plot lines, as a multiple of the default}

!!table

!!end

*/

#include <cstdlib>
#include <iostream>
#include <string>
#include "cpgplot.h"
#include "trm/subs.h"
#include "trm/format.h"
#include "trm/plot.h"
#include "trm/input.h"
#include "trm/vec3.h"
#include "trm/binary.h"
#include "trm/roche.h"

int main (int argc, char *argv[]){

  try{

    // Construct Input object

    Subs::Input input(argc, argv, Roche::ROCHE_ENV, Roche::ROCHE_DIR);

    // Define inputs
    input.sign_in("m1",      Subs::Input::LOCAL,  Subs::Input::PROMPT);
    input.sign_in("m2",      Subs::Input::LOCAL,  Subs::Input::PROMPT);
    input.sign_in("type",    Subs::Input::LOCAL,  Subs::Input::PROMPT);
    input.sign_in("device",  Subs::Input::GLOBAL, Subs::Input::NOPROMPT);
    input.sign_in("x1",      Subs::Input::LOCAL,  Subs::Input::PROMPT);
    input.sign_in("x2",      Subs::Input::LOCAL,  Subs::Input::PROMPT);
    input.sign_in("y1",      Subs::Input::LOCAL,  Subs::Input::PROMPT);
    input.sign_in("y2",      Subs::Input::LOCAL,  Subs::Input::PROMPT);
    input.sign_in("width",   Subs::Input::LOCAL,  Subs::Input::PROMPT);
    input.sign_in("aspect",  Subs::Input::LOCAL,  Subs::Input::PROMPT);
    input.sign_in("lwidth",  Subs::Input::LOCAL,  Subs::Input::PROMPT);
    input.sign_in("reverse", Subs::Input::LOCAL,  Subs::Input::PROMPT);

    float m1;
    input.get_value("m1", m1, 0.6f, 1.e-3f, 1.4f, "primary mass (solar)");
    float m2;
    input.get_value("m2", m2, min(0.1f, m1), 1.e-3f, m1, "secondary mass (solar)");
    const double q  = m2/m1;
    std::string type; 
    input.get_value("type", type, "cv", "system type ('cv' or 'am')");
    if(Subs::toupper(type) != "CV" && Subs::toupper(type) != "AM")
      throw Roche::Roche_Error("type must be either 'cv' or 'am'");
    std::string device; 
    input.get_value("device", device, "/xs", "plot device");
    float x1, x2, y1, y2;
    input.get_value("x1", x1, -1.f, -20.f, 20.f, "left-hand X limit");
    input.get_value("x2", x2, +1.f, -20.f, 20.f, "right-hand X limit");
    input.get_value("y1", y1, -1.f, -20.f, 20.f, "lower Y limit");
    input.get_value("y2", y2, +1.f, -20.f, 20.f, "upper Y limit");
    float width;
    input.get_value("width", width, 6.f, 0.f, 100.f, "plot width (inches)");
    float aspect;
    input.get_value("aspect", aspect, 0.6f, 0.01f, 100.f, "aspect ratio (height/width)");
    int  lwidth;
    input.get_value("lwidth", lwidth, 2, 1, 100, "line width");
    bool reverse;
    input.get_value("reverse", reverse, false, "reverse foreground/background colours");

    const int NPLOT = 200;
    float x[NPLOT], y[NPLOT];

    Subs::Plot plot(device);
    cpgpap(width, aspect);
    if(reverse){
      float bred, bgreen, bblue, fred, fgreen, fblue;
      cpgqcr(0, &bred, &bgreen, &bblue);
      cpgqcr(1, &fred, &fgreen, &fblue);
      cpgscr(1, bred,  bgreen, bblue);
      cpgscr(0, fred,  fgreen, fblue);
    }
    cpgsvp(0., 1., 0., 1.);
    cpgwnad(x1,x2,y1,y2);
    cpgslw(lwidth);
    cpgsch(1.5);
    cpgscf(2);
    cpgsci(4);

    // Secondary star
    cpgsci(2);
    cpgsfs(3);
    Roche::lobe2(q, x, y, NPLOT);
    cpgpoly(NPLOT, x, y);
    cpgline(NPLOT, x, y);

    // Primary lobe
    cpgsci(1);
    cpgsls(2);
    Roche::lobe1(q, x, y, NPLOT);
    cpgline(NPLOT, x, y);

    const float r1     = Binary::mr_wd_eggleton(m1);
    float r2;
    if(Subs::toupper(type) == "CV")
      r2 = Binary::mr_main_sequence(m2);
    else
      r2 = Binary::mr_wd_eggleton(m2);

    const float rl2    = Roche::rlobe_eggleton(q);
    const float a      = r2/rl2;
    const float period = Binary::orbital_period(m1,m2,a);
    Subs::Format form(10);

    std::cout << "Period = " << form(period) << " seconds." << std::endl;
    int iperiod = int(period+0.5);

    // Primary star. Find potential of a point 45 degrees around
    Subs::Vec3 p;
    p.x()   = r1/a/sqrt(2.);
    p.y()   = r1/a/sqrt(2.);
    p.z()   = 0.;
    double pot = Roche::rpot(q,p);
    Roche::flobe1(q, x, y, NPLOT, pot);
    cpgsfs(3);
    cpgsls(1);
    cpgpoly(NPLOT, x, y);
    cpgline(NPLOT, x, y);

    std::cout << "r1/a = " << r1/a << std::endl;
    // Stream
    double xh, yh;

    if(Roche::hits(q,pot,xh,yh)){

      Roche::streamr(q,sqrt(xh*xh+yh*yh),x,y,NPLOT);
      cpgline(NPLOT, x, y);
      cpgsci(3);
      cpgsch(2);
      cpgpt(1, x+NPLOT-1, y+NPLOT-1, 18);

      float xspot = x[NPLOT-1];
      float yspot = y[NPLOT-1];

      const int NARC = 200;
      float xarc[NARC], yarc[200];
      double angle1 = atan2(yspot, xspot), angle;
      double rad = 1.15*sqrt(xspot*xspot+yspot*yspot);
      const double aadd[3] = {1.4, 0.9, 0.4};
      const int ci[3] = {2, 3, 4};
      const double DELTA = 0.02;
      double ca, sa;
      cpgsfs(1);
      for(int j=0; j<3; j++){
	for(int i=0; i<NARC/2; i++){
	  angle = angle1 + aadd[j]*i/(NARC/2-1);
	  ca = cos(angle);
	  sa = sin(angle);
	  xarc[i]        = (rad+DELTA)*ca;
	  xarc[NARC-1-i] = rad*ca;
	  yarc[i]        = (rad+DELTA)*sa;
	  yarc[NARC-1-i] = rad*sa;
	}
	cpgsci(ci[j]);
	cpgpoly(NARC, xarc, yarc);
	rad += 1.5*DELTA;
      }
      rad *= 1.2;
      cpgsch(1.2);
      for(int j=0; j<3; j++){
	angle = angle1 + aadd[j]/2.;
	ca    = rad*cos(angle);
	sa    = rad*sin(angle);
	cpgsci(ci[j]);
	cpgarro(ca, sa, 1.5*ca, 1.5*sa);
      }

    }else{
    
      float rd = 0.65*Roche::xl1(q);
      cpgsfs(2);
      cpgcirc(0.,0.,rd);
      Roche::streamr(q,rd,x,y,NPLOT);
      cpgline(NPLOT, x, y);
      cpgsci(3);
      cpgsch(2);
      cpgpt(1, x+NPLOT-1, y+NPLOT-1, 18);

    }

    // Top label
    cpgsci(2);
    std::string title = "M1 = " + Subs::str(m1) + ", M2 = " + Subs::str(m2) +
      ", P = " + Subs::str(iperiod) + " sec";

    cpglab(" ", " ", title.c_str());

    
  }

  catch(const std::string& err){
    std::cerr << "string exception:" << std::endl;
    std::cerr << err << std::endl;
  }
  catch(...){
    std::cerr << "unknown error occurred" << std::endl;
  }
}









