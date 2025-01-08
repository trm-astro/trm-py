/*

!!begin
!!title   algol -- computes algol stuff
!!author  T.R. Marsh
!!created Nov 2001
!!descr   computes algol stuff to do with V403 Vul
!!root    algol
!!index   algol
!!class   Programs
!!css     style.css
!!head1   algol - direct impact computation.

This program computes critical parameters for direct impact between
two white dwarfs. For each primary mass, it prints out the critical donor
mass for impact to occur and the critical donor mass at which impact
occurs and the donor just sees the impact spot.

!!table
!!arg{m1low}{Mass of primary, lower limit, solar masses}
!!arg{m1high}{Mass of primary, upper limit, solar masses}
!!arg{nmass}{Number of masses}
!!table

!!end

*/

#include <climits>
#include <cstdlib>
#include <iostream>
#include "trm/subs.h"
#include "trm/input.h"
#include "trm/binary.h"
#include "trm/roche.h"

int main (int argc, char *argv[]){

  try{

    // Construct Input object

    Subs::Input input(argc, argv, Roche::ROCHE_ENV, Roche::ROCHE_DIR);

    // Define inputs

    input.sign_in("m1low",  Subs::Input::LOCAL, Subs::Input::PROMPT);
    input.sign_in("m1high", Subs::Input::LOCAL, Subs::Input::PROMPT);
    input.sign_in("nmass",  Subs::Input::LOCAL, Subs::Input::PROMPT);

    float m1lo;
    input.get_value("m1low",  m1lo, 0.1f, 1.e-3f, 1.4f, "primary mass, lower limit (solar)");
    float m1hi;
    input.get_value("m1high", m1hi, max(1.2f, m1lo), m1lo, 1.4f, "primary mass, upper limit (solar)");
    int nm;
    input.get_value("nmass", nm, 100, 1, INT_MAX, "number of masses");
    
    double q, qlo, qhi, qc1, rl2, a;
    double x, y, m1, m2, r1, r2, pot, xh, yh;
    
    // Primary star. Find potential of a point 45 degrees around

    Subs::Vec3 p;
    for(int j=0; j<nm; j++){

      qlo = 0.0001;
      qhi = 2.0;

      m1  = m1lo + (m1hi-m1lo)*j/(nm-1);
      r1  = Binary::mr_wd_eggleton(m1);

      // First check lower limit
      m2  = m1/qlo;
      r2  = Binary::mr_wd_eggleton(m2);
      rl2 = Roche::rlobe_eggleton(qlo);
      a   = r2/rl2;

      p.x()   = r1/a/sqrt(2.);
      p.y()   = r1/a/sqrt(2.);
      p.z()   = 0.;
      pot   = Roche::rpot(qlo,p);

      if(Roche::hits(qlo,pot,x,y))
	throw Roche::Roche_Error("For q = " + Subs::str(qlo) + ", stream unexpectedly hits!");

      // then upper limit
      m2  = m1/qhi;
      r2  = Binary::mr_wd_eggleton(m2);
      rl2 = Roche::rlobe_eggleton(qhi);
      a   = r2/rl2;

      p.x()   = r1/a/sqrt(2.);
      p.y()   = r1/a/sqrt(2.);
      pot   = Roche::rpot(qhi,p);

      if(!Roche::hits(qhi,pot,xh,yh))
	throw Roche::Roche_Error("For q = " + Subs::str(qhi) + ", stream unexpectedly does not hit!");

      while(qhi-qlo > 1.e-4){
	q   = (qlo+qhi)/2.;
	m2  = m1/q;
	r2  = Binary::mr_wd_eggleton(m2);
	rl2 = Roche::rlobe_eggleton(qhi);
	a   = r2/rl2;
	p.x() = r1/a/sqrt(2.);
	p.y() = r1/a/sqrt(2.);
	pot = Roche::rpot(q,p);
	if(Roche::hits(q,pot,x,y)){
	  qhi = q;
	}else{
	  qlo = q;
	}
      }
      qc1 = qhi;
      if(Roche::irrad(qc1,x,y)){
	std::cout << m1 << " " << m1*qc1 << " " << 0. << std::endl;
      }else{
	qlo = qc1;
	qhi = 2.0;
	if(!Roche::irrad(qhi,xh,yh)){
	  std::cout << m1 << " " << m1*qc1 << " " << 1. << std::endl;
	}else{
	  while(qhi-qlo > 1.e-4){
	    q = (qlo+qhi)/2.;
	    m2  = m1/q;
	    r2  = Binary::mr_wd_eggleton(m2);
	    rl2 = Roche::rlobe_eggleton(qhi);
	    a   = r2/rl2;
	    p.x() = r1/a/sqrt(2.);
	    p.y() = r1/a/sqrt(2.);
	    pot = Roche::rpot(q,p);
	    Roche::hits(q,pot,x,y);
	    if(Roche::irrad(q,x,y)){
	      qhi = q;
	    }else{
	      qlo = q;
	    }
	  }
	  std::cout << m1 << " " << m1*qc1 << " " << m1*qlo << std::endl;
	}
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









