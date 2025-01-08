/*

!!begin
!!title   pphi - Roche potential plot
!!author  T.R. Marsh
!!created Aug 2003
!!descr   plots Roche potential versus lambda
!!root    pphi
!!index   pphi
!!class   Programs
!!css     style.css
!!head1   pphi - plots Roche potential

Given a mass ratio, orbital inclination and filling factor this program
computes and plots the Roche potential along the line of sight for a given
orbital phase, IF it crosses the reference sphere around the roche distorted
star.

!!end

*/

#include <iostream>
#include "cpgplot.h"
#include "trm/subs.h"
#include "trm/constants.h"
#include "trm/plot.h"
#include "trm/vec3.h"
#include "trm/roche.h"

int main(){
    Subs::Vec3 r;
    double q, iangle, ffac, lam1, lam2;
    std::cout << "Enter q, i and filling factor: ";
    std::cin >> q >> iangle >> ffac;
    double rsphere, pref, phase;
    Roche::ref_sphere(q, Roche::SECONDARY, 1., ffac, rsphere, pref);
    const Subs::Vec3 cofm2(1.,0.,0.);
    double ri = Subs::deg2rad(iangle);
    
    for(;;){
	
	std::cout << "Enter triplet of numbers for point position and phase: ";
	std::cin >> r >> phase;
	
	Subs::Vec3 earth = Roche::set_earth(iangle, phase);
	
	try{
	    if(Roche::sphere_eclipse(earth, r, cofm2, rsphere, lam1, lam2)){
		std::cout << "Eclipsed: " << lam1 << " " << lam2 << std::endl;
		const int NX=200;
		float x[NX], y[NX], pmin = 1.e30, pmax = -1.e30;
		for(int ix=0; ix<NX; ix++){
		    double lam = lam1 + (lam2-lam1)*ix/double(NX-1);
		    Subs::Vec3 p = r + lam*earth;
		    
		    x[ix] = lam;
		    y[ix] = Roche::rpot(q,p);
		    pmin   = std::min(y[ix], pmin);
		    pmax   = std::max(y[ix], pmax);
		}
		float range = pmax - pmin;
		pmin -= 0.05*range;
		pmax += 0.05*range;
	
		Subs::Plot plot("/xs");
		cpgslw(2);
		cpgsch(1.5);
		cpgscf(2);
		cpgsci(4);
		cpgenv(lam1,lam2,pmin,pmax, 0,0);
		cpgsci(2);
		cpgline(NX, x, y);
		cpgsls(2);
		cpgmove(lam1, pref);
		cpgdraw(lam2, pref);
		
	    }else{
		std::cout << "no eclipse" << std::endl;
	    }

	}
	catch(Roche::Roche_Error err){
	    std::cerr << "Roche::Roche_Error\n" << err.what() << std::endl;
	}
    }
}
