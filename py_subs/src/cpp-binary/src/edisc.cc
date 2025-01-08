/*

!!begin
!!title   edisc -- computes a line profile from an elliptical disc
!!author  T.R. Marsh
!!created 14 Aug 2006
!!descr   computes a line profile from an elliptical disc
!!root    edisc
!!index   edisc
!!class   Programs
!!css     style.css
!!head1   edisc - computes a line profile from an elliptical disc

!!emph{edisc} computes the emission line profile from a disc consisting of a series
of concentric elliptical orbits of the same eccentricity and scaled only by semi-major
axis. 

!!table
!!arg{m}{Mass of central object, solar masses.}
!!arg{e}{eccentricity of orbits}
!!arg{alpha}{Orientation of the ellipses to the line of sight. If alpha = 0, the observer is looking along the
major axis, with the periastron furthest from him/her. As alpha becomes greater than zero, the periastron region
rotates so that it becomes blue-shifted.}
!!arg{a1}{Inner semi-major axis, solar radii.}
!!arg{a2}{Outer semi-major axis, solar radii.}
!!arg{dr}{Radial width of annuli}
!!arg{beta}{Emmissivity scales at r**beta}
!!arg{data}{Data file to load, 'none' if none. Velocity, flux, error column data otherwise. Velocity in km/s}
!!arg{vmax}{If no data file loaded. Maximum velocity to store (km/s)}
!!arg{vstep}{If no data file loaded. Velocity step (determines number of pixels) (km/s)}
!!arg{nsplit}{Splitting factor for velocity pixels}
!!arg{fwhm}{FWHM resolution (km/s)}
!!arg{q}{Factor controlling the importance of shear broadening in the disc, roughly equivalent to
Q*sin(i)*tan(i) from Horne & Marsh 1986}
!!arg{az}{Amount of azimuthal asymmetry to impose. Multiplies emissivity by (1+az*cos(theta))
where theta is angle from periastron. az can lie from -1 to +1.}
!!arg{plot}{Yes or no plot to /xs}
!!table

!!end

*/

#include <cstdlib>
#include <iostream>
#include "cpgplot.h"
#include "trm/subs.h"
#include "trm/plot.h"
#include "trm/array1d.h"
#include "trm/format.h"
#include "trm/input.h"
#include "trm/binary.h"

int main (int argc, char *argv[]){

    try{

        // Construct Input object

        Subs::Input input(argc, argv, Binary::BINARY_ENV, Binary::BINARY_DIR);

        // Define inputs

        input.sign_in("m",       Subs::Input::LOCAL, Subs::Input::PROMPT);
        input.sign_in("e",       Subs::Input::LOCAL, Subs::Input::PROMPT);
        input.sign_in("alpha",   Subs::Input::LOCAL, Subs::Input::PROMPT);
        input.sign_in("a1",      Subs::Input::LOCAL, Subs::Input::PROMPT);
        input.sign_in("a2",      Subs::Input::LOCAL, Subs::Input::PROMPT);
        input.sign_in("dr",      Subs::Input::LOCAL, Subs::Input::PROMPT);
        input.sign_in("beta",    Subs::Input::LOCAL, Subs::Input::PROMPT);
        input.sign_in("data",    Subs::Input::LOCAL, Subs::Input::PROMPT);
        input.sign_in("vmax",    Subs::Input::LOCAL, Subs::Input::PROMPT);
        input.sign_in("vstep",   Subs::Input::LOCAL, Subs::Input::PROMPT);
        input.sign_in("nsplit",  Subs::Input::LOCAL, Subs::Input::PROMPT);
        input.sign_in("fwhm",    Subs::Input::LOCAL, Subs::Input::PROMPT);
        input.sign_in("q",       Subs::Input::LOCAL, Subs::Input::PROMPT);
        input.sign_in("az",      Subs::Input::LOCAL, Subs::Input::PROMPT);
        input.sign_in("plot",    Subs::Input::LOCAL, Subs::Input::PROMPT);

        Subs::Format form(10);

        double m;
        input.get_value("m",  m, 0.5, 1.e-5, 1.e5, "mass of central object (solar masses)");
        double e;
        input.get_value("e",  e, 0.5, 0., 0.99999999, "eccentricity of orbits");
        double alpha;
        input.get_value("alpha", alpha, 80., -10000., 10000., "argument of pericentre (degrees)");
        double a1;
        input.get_value("a1",  a1, 0.1, 1.e-5, 1.e5, "inner semi-major axis");
        double a2;
        input.get_value("a2",  a2, std::max(a1, 100.), a1, 1.e8, "outer semi-major axis");
        double dr;
        input.get_value("dr",  dr, 0.01, 1.e-4, 10., "radial width of annuli");
        double beta;
        input.get_value("beta", beta, -2., -100., 100., "emissivity power-law exponent");
        std::string sdata;
        double vstart, vstep;
        int npix;
        input.get_value("data", sdata, "none", "data file to store (km/s)");
        Subs::Array1D<Subs::xyz<double,float,float> > data;
        if(sdata == "none"){
            double vmax;
            input.get_value("vmax", vmax, 2000., 1., 100000., "maximum velocity to store (km/s)");
            input.get_value("vstep", vstep, 100., 1., 10000., "velocity step size (km/s)");
            npix   = int(2*vmax/vstep) + 1;
            vstart = -vstep*npix/2;
        }else{
            data.load_ascii(sdata);
            vstart = data[0].x;
            npix   = data.size();
            double vend = data[npix-1].x;
            vstep = (vend-vstart)/(npix-1);
            double dmax = 0.;
            for(int i=0; i<npix; i++){
                double dev = Subs::abs((data[i].x - (vstart + vstep*i))/vstep);
                if(dev > dmax) dmax = dev;
            }
            std::cerr << "Maximum deviation from linearity = " << form(dmax) << " pixels" << std::endl;
            vstart -= vstep/2.;
            std::cerr << "Velocity step size = " << form(vstep) << " km/s, number of pixels = " << npix << std::endl;
            std::cerr << "Velocity range = " << form(vstart) << " to " << form(vstart + npix*vstep) << std::endl;
        }

        int nsplit;
        input.get_value("nsplit", nsplit, 5, 1, 1000, "sub-division factor");
        const int NFINE = nsplit*npix;

        double fwhm;
        input.get_value("fwhm", fwhm, 200., 1., 10000., "FWHM resolution (km/s)");
        double q;
        input.get_value("q", q, 0., 0., 1000., "importance of shear broadening");
        double az;
        input.get_value("az", az, 0., -1., 1., "azimuthal asymmetry");
        bool mplot;
        input.get_value("plot", mplot, false, "do you want a plot?");

        Subs::Array1D<double> fprof(NFINE);

        fprof = 0.;

        const int NA = int((a2-a1)/dr) + 1;
        Binary::Einfo einfo;
        for(int n=0; n<NA; n++){
            double a    = a1 + (a2-a1)*(n+0.5)/NA;
            int ntheta = int(Constants::TWOPI*a/dr);
            ntheta = ntheta > 32 ? ntheta : 32;
            for(int nt=0; nt<ntheta; nt++){

                double theta = Constants::TWOPI*nt/ntheta;

                // Compute RV information
                einfo = Binary::vellipse(m, e, a, theta);

                double area  = Constants::TWOPI/ntheta*dr/a*Subs::sqr(einfo.r);
                double emiss = area*pow(einfo.r,beta)*(1+az*cos(theta));

                double cosa  = cos(Constants::TWOPI*alpha/360);
                double sina  = sin(Constants::TWOPI*alpha/360);
                double rv    = einfo.vx*cosa + einfo.vy*sina;
                double shear = (einfo.vx_x*cosa*cosa + (einfo.vx_y+einfo.vy_x)*cosa*sina + einfo.vy_y*sina*sina)/einfo.omega;

                int nadd = int(nsplit*(rv-vstart)/vstep);

                if(nadd >=0 && nadd < NFINE) fprof[nadd] += sqrt(1+Subs::sqr(q*shear))*emiss;
            }
        }

        // create blurr array
        const int NBLURR = int(4.*nsplit*fwhm/vstep + 1);
        Subs::Array1D<double> blurr(2*NBLURR+1);
        double norm = 0.;
        for(int nb=0; nb<2*NBLURR+1; nb++){
            blurr[nb] = exp(-Subs::sqr(vstep*(nb-NBLURR)/nsplit/(fwhm/2.3548))/2.);
            norm += blurr[nb];
        }
        blurr /= norm;

        // blurr array
        Subs::Array1D<double> tfprof(NFINE);
        tfprof = 0.;
        for(int nf=0; nf<NFINE; nf++){
            for(int nb=0, nadd=nf-NBLURR; nb<2*NBLURR+1; nb++, nadd++){
                if(nadd >=0 && nadd < NFINE)
                    tfprof[nadd] += blurr[nb]*fprof[nf];
            }
        }

        // Bin and print out
        Subs::Array1D<float>  x(npix);
        Subs::Array1D<double> y(npix);
        y = 0.;
        double sum1 = 0., sum2 = 0.;
        for(int np=0; np<npix; np++){
            x[np] = vstart + vstep*(np+0.5);
            for(int nfp=nsplit*np; nfp<nsplit*(np+1); nfp++)
                y[np] += tfprof[nfp];

            y[np] /= nsplit;

            // derive optimum scaling factor
            if(sdata != "none"){
                sum1 += data[np].y*y[np]/Subs::sqr(data[np].z);
                sum2 += Subs::sqr(y[np]/data[np].z);
            }

        }

        // apply scaling
        if(sdata != "none") y *= sum1/sum2;
      
        Subs::Array1D<double> dy(npix);
        if(sdata != "none"){
            double chisq = 0.;
            for(int np=0; np<npix; np++){
                dy[np] = data[np].y;
                chisq += Subs::sqr((dy[np]-y[np])/data[np].z);
            }
            std::cerr << "Chisq = " << form(chisq) << std::endl;
        }

        // plot
        float ymax = 1.1*y.max();
        if(sdata != "none") ymax = std::max(ymax, float(1.1*dy.max()));
        if(mplot){
            Subs::Plot plot("/xs");
            cpgsch(1.5);
            cpgslw(3);
            cpgenv(vstart, vstart + npix*vstep, -0.1*ymax, ymax, 0, 0);
            cpglab("Velocity (km/s)", "Flux density", " " );
            if(sdata != "none"){
                cpgsci(2);
                for(int np=0; np<npix; np++){
                    cpgmove(x[np], dy[np] - data[np].z);
                    cpgdraw(x[np], dy[np] + data[np].z);
                }
                cpgsci(1);
                pgpt(x, dy, 17);
                cpgsls(4);
                cpgsci(2);
            }
	    
            // plot model
            pgline(x, y);
        }   
 
        // Output the fit.
        for(int np=0; np<npix; np++)
            std::cout << form(x[np]) << " " << y[np] << std::endl; 
	    
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









