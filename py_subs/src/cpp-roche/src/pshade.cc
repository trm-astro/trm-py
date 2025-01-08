/*

!!begin
!!title   pshade -- plots Roche lobes and shadows
!!author  T.R. Marsh
!!created 30 March 2006
!!descr   plots Roche lobes and shadows
!!root    pshade
!!index   pshade
!!class   Programs
!!css     style.css
!!head1   pshade - direct impact stuff

!!emph{pdwd} plots the stream in a double white dwarf binary.


!!table

!!arg{q}{Mass ratio, M2/M1}

!!arg{iangle}{Inclination angle}

!!arg{phase}{Orbital phase}

!!arg{r1}{Radius of primary star, units of separation}

!!arg{rdisc}{Radius of disc, units of separation}

!!arg{acc}{Accuracy of position determination}

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
#include "trm/plot.h"
#include "trm/input.h"
#include "trm/vec3.h"
#include "trm/binary.h"
#include "trm/roche.h"

using Subs::operator+;

int main (int argc, char *argv[]){

  try{

    // Construct Input object

    Subs::Input input(argc, argv, Roche::ROCHE_ENV, Roche::ROCHE_DIR);

    // Define inputs
    input.sign_in("q",       Subs::Input::LOCAL,  Subs::Input::PROMPT);
    input.sign_in("iangle",  Subs::Input::LOCAL,  Subs::Input::PROMPT);
    input.sign_in("phase",   Subs::Input::LOCAL,  Subs::Input::PROMPT);
    input.sign_in("r1",      Subs::Input::LOCAL,  Subs::Input::PROMPT);
    input.sign_in("rdisc",   Subs::Input::LOCAL,  Subs::Input::PROMPT);
    input.sign_in("acc",     Subs::Input::LOCAL,  Subs::Input::PROMPT);
    input.sign_in("device",  Subs::Input::GLOBAL, Subs::Input::NOPROMPT);
    input.sign_in("x1",      Subs::Input::LOCAL,  Subs::Input::PROMPT);
    input.sign_in("x2",      Subs::Input::LOCAL,  Subs::Input::PROMPT);
    input.sign_in("y1",      Subs::Input::LOCAL,  Subs::Input::PROMPT);
    input.sign_in("y2",      Subs::Input::LOCAL,  Subs::Input::PROMPT);
    input.sign_in("width",   Subs::Input::LOCAL,  Subs::Input::PROMPT);
    input.sign_in("aspect",  Subs::Input::LOCAL,  Subs::Input::PROMPT);
    input.sign_in("lwidth",  Subs::Input::LOCAL,  Subs::Input::PROMPT);
    input.sign_in("reverse", Subs::Input::LOCAL,  Subs::Input::PROMPT);

    double q;
    input.get_value("q", q, 0.2, 0.001, 10., "mass ratio, (M2/M1)");
    double iangle;
    input.get_value("iangle", iangle, 82., 0., 89.9999, "inclination angle (deg)");
    double phase;
    input.get_value("phase", phase, 0., -10000., 10000., "orbital phase (deg)");
    double r1;
    input.get_value("r1", r1, 0.01, 0., 1., "radius of primary (units of a)");
    double rdisc;
    input.get_value("rdisc", rdisc, 0.3, 0., 1., "radius of disc (units of a)");
    double acc;
    input.get_value("acc", acc, 1.e-4, 1.e-10, 0.01, "accuracy of position determinations (units of a)");
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

    const int NPLOT = 400;
    float x[NPLOT], y[NPLOT];
    bool shade[NPLOT];

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

    cpgsls(1);

    // Primary star
    cpgsfs(1);
    cpgsci(4);
    cpgcirc(0.,0.,r1);

    // Shadow of primary
    double fac = 1./cos(M_PI*iangle/180);
    double ct = cos(2*M_PI*phase);
    double st = sin(2*M_PI*phase);
    for(int i=0; i<NPLOT; i++){
      double theta = M_PI/2 + M_PI*i/(NPLOT-1);
      double xt = r1*fac*cos(theta);
      double yt = r1*sin(theta);
      x[i] =   ct*xt + st*yt;
      y[i] =  -st*xt + ct*yt;
    }
    cpgsls(4);
    cpgline(NPLOT, x, y);

    // Stream
    Roche::streamr(q,rdisc,x,y,NPLOT);
    cpgsls(1);
    cpgsci(10);
    cpgline(NPLOT, x, y);
    cpgsch(2.5);
    cpgpt(1, x+NPLOT-1, y+NPLOT-1, 18);

    // Disc
    cpgsfs(2);
    cpgcirc(0.,0.,rdisc);
    
    // Shadow
    cpgsls(2);
    cpgsci(2);
    Roche::roche_shadow(q, iangle, phase, 10., acc, x, y, shade, NPLOT);
    bool draw = false;
    for(int i=0; i<NPLOT; i++){
      if(shade[i]){
	if(draw){ 
	  cpgdraw(x[i], y[i]);
	}else{
	  draw = true;
	  cpgmove(x[i],y[i]);
	}
      }else{
	draw = false;
      }
    }
  }

  catch(const std::string& err){
    std::cerr << "string exception:" << std::endl;
    std::cerr << err << std::endl;
  }
  catch(...){
    std::cerr << "unknown error occurred" << std::endl;
  }
}









