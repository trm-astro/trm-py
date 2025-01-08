/*

!!begin
!!title   pcont - Roche potential contours
!!author  T.R. Marsh
!!created Aug 2003
!!descr   plots contours of Roche potential
!!root    pcont
!!index   pcont
!!class   Programs
!!css     style.css
!!head1   pcont - contours of Roche potential

Given a mass ratio, orbital inclination and filling factor this program
computes phases and vector multpliers within which, if the point is ever
eclipsed, the line of sight from it dips below the stellar photosphere.
It then computes a grid between these values which it plots as a greyscale
and contours with a green contour for the critical Roche potential. It keeps
prompting for new points until you give up.

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
    double q, iangle, ffac, phi1, phi2, lam1, lam2;
    std::cout << "Enter q, i and filling factor: ";
    std::cin >> q >> iangle >> ffac;
    double rref, pref;
    Roche::ref_sphere(q, Roche::SECONDARY, 1., ffac, rref, pref);
    const Subs::Vec3 cofm2(1.,0.,0.);
    double ri = Subs::deg2rad(iangle);
    double cosi = cos(ri), sini = sin(ri);
    for(;;){
    
	std::cout << "Enter triplet of numbers for point position: ";
	std::cin >> r;

	try{
	    if(Roche::sphere_eclipse(cosi, sini, r, cofm2, rref, phi1, phi2, lam1, lam2)){
		std::cout << "Eclipsed by sphere: " << phi1 << " " << phi2 << " " << lam1 << " " << lam2 << std::endl;
		const int NX=100, NY=100;
		float arr[NX*NY], pmin = 1.e30, pmax = -1.e30;
		for(int iy=0, i=0; iy<NY; iy++){
		    double lam = lam1 + (lam2-lam1)*iy/double(NY-1);
		    for(int ix=0; ix<NX; ix++, i++){ 
			double phi = phi1 + (phi2-phi1)*ix/double(NX-1);
			Subs::Vec3 earth = Roche::set_earth(cosi, sini, phi);
			Subs::Vec3 x = r + lam*earth; 
			arr[i] = Roche::rpot(q, x);
			pmin   = std::min(arr[i], pmin);
			pmax   = std::max(arr[i], pmax);
		    }
		}
		Subs::Plot plot("/xs");
		cpgslw(2);
		cpgsch(1.5);
		cpgscf(2);
		cpgsci(4);
		cpgenv( phi1, phi2, lam1, lam2, 0, 0);
		cpglab("Phi","Lambda","");
		cpgsci(2);
		float tr[6];
		tr[1] = (phi2-phi1)/double(NX-1);
		tr[0] = phi1 - tr[1];
		tr[2] = 0;
		tr[5] = (lam2-lam1)/double(NY-1);
		tr[3] = lam1 - tr[5];
		tr[4] = 0.;
		cpggray(arr, NX, NY, 1, NX, 1, NY, pmax, pmin,  tr);
		const int NCONT = 20;
		float cont[NCONT];
		for(int i=0; i<NCONT; i++) cont[i] = pmin + (pmax-pmin)*(i+1)/float(NCONT+1);
		cpgsci(2);
		cpgcont(arr, NX, NY, 1, NX, 1, NY, cont, NCONT,  tr);
		cont[0] = pref;
		cpgsci(3);
		cpgcont(arr, NX, NY, 1, NX, 1, NY, cont, 1,  tr);
		
		double phi, lam;
		Roche::pot_min(q, cosi, sini, Roche::SECONDARY, 1., r, phi1, phi2, lam1, lam2, rref, pref, 1.e-3, phi, lam);
		cpgpt1(phi, lam, 5);
		
		double ingress, egress;
		if(Roche::ingress_egress(q, Roche::SECONDARY, 1., ffac, iangle, 1.e-6, r, ingress, egress)){
		    std::cout << "point is eclipsed by star from " << ingress << " to " << egress << std::endl;
		}else{
		    std::cout << "point is not eclipsed by star" << std::endl;
		}
	    }else{
		std::cout << "no eclipse by sphere" << std::endl;
	    }
	}
	catch(Roche::Roche_Error err){
	    std::cerr << "Roche::Roche_Error\n" << err.what() << std::endl;
	}
    }
}
