#include <cmath>
#include <cstdlib>
#include <iostream>
#include "trm/subs.h"
#include "trm/constants.h"
#include "trm/vec3.h"
#include "trm/roche.h"

/**
 * The line of sight to any fixed point in a binary sweeps out a cone at the
 * binary rotates. Positions on the cone can be parameterised by the orbital
 * phase phi and the multiplier ('lambda') needed to get from the fixed point.
 * The question pot_min tries to solve is "does the cone intersect a surface
 * of fixed Roche potential lying within a Roche lobe?". It does so by
 * minimisation over a region of phi and lambda. It stops as soon as any
 * potential below a critical value is found. The initial range of phi and
 * lambda can be determined using sphere_eclipse which calculates them for a
 * sphere.
 *
 * \param q      mass ratio = M2/M1
 * \param cosi   cosine orbital inclination
 * \param sini   sine orbital inclination
 * \param star   which star (needed for asynchronous case)
 * \param spin   ratio of spin/orbital
 * \param p      point of origin
 * \param phi1   minimum phase within which eclipse may occur (0 - 1)
 * \param phi2   maximum phase within which an eclipse may occur (> phi1)
 * \param lam1   minimum multiplier line of sight crosses eclipsing star
 * \param lam2   maximum multiplier line of sight crosses eclipsing star
 * \param rref   reference radius
 * \param pref   reference potential.
 * \param acc    absolute accuracy in position to go for
 * \param phi    phi at minimum potential. Ingress occurs between phi1 and phi
                 if there is an eclipse. Egress occurs between phi and phi2.
 * \param lam    lambda at minimum potential
 * \return true if minimum potential is below the reference
 */

bool Roche::pot_min(double q, double cosi, double sini, STAR star, double spin,
                    const Subs::Vec3& p, double phi1, double phi2, double lam1,
                    double lam2, double rref, double pref, double acc,
                    double& phi, double& lam){

    if(q <= 0.) throw Roche_Error("q = " + Subs::str(q) + "(<= 0.) in pot_min");

    double dphi, dlam, gdphi, gdlam, hdphi, hdlam;

    // Start in the middle
    phi = (phi1+phi2)/2.;
    lam = (lam1+lam2)/2.;

    double rp = Constants::TWOPI*phi;
    double cosp = cos(rp), sinp = sin(rp);
    Subs::Vec3 earth(sini*cosp, -sini*sinp, cosi);
    double pot;
    rpot_val_grad(q, star, spin, earth, p, lam, pot, dphi, dlam);
    if(pot <= pref) return true;

    gdphi = -dphi;
    gdlam = -dlam;
    dphi  = hdphi = gdphi;
    dlam  = hdlam = gdlam;

    const int ITMAX = 200;

    // Accuracy in position corresponds to an accuracy in potential
    const double DELPHI = q/(1.+q)*Subs::sqr(acc/rref)/2.;
    double pmin, gam, dgg, gg;
    bool jammed;
    for(int its=0; its<ITMAX; its++){

        linmin(q, star, spin, cosi, sini, p, phi, lam, dphi, dlam,
               phi1, phi2, lam1, lam2, pref, acc, pmin, jammed);

        // Various reasons for stopping
        if(pmin <= pref) return true;
        if(jammed || fabs(pmin-pot) < DELPHI) return false;

        pot  = pmin;
        rp   = Constants::TWOPI*phi;
        cosp = cos(rp);
        sinp = sin(rp);
        earth.set(sini*cosp, -sini*sinp, cosi);
        rpot_grad(q, star, spin, earth, p, lam, dphi, dlam);

        gg  = gdphi*gdphi + gdlam*gdlam;
        if(gg == 0.) return false;

        dgg = (dphi+gdphi)*dphi + (dlam+gdlam)*dlam;
        gam = dgg/gg;

        gdphi = -dphi;
        gdlam = -dlam;
        dphi  = hdphi = gdphi + gam*hdphi;
        dlam  = hdlam = gdlam + gam*hdlam;

    }

    throw Roche_Error("Roche::pot_min: too many iterations.");

}

/*
 * linmin minimises along a line in phase, lambda space. It returns at the
 * minimum or as soon as the potential drops below a reference value. It is
 * assumed that the potential is dropping with x at the starting point.
 *
 * \param q      mass ratio = M2/M1
 * \param star   star in question
 * \param spin   spin to orbital frequency ratio
 * \param cosi   cosine of orbital inclination (both passed to speed computations)
 * \param sini   sine of orbital inclination
 * \param p      point of origin
 * \param phi    value of phase at start (0 - 1) and at end (returned)
 * \param lam    value of lambda at start and at end (returned)
 * \param dphi   rate at which phi changes
 * \param dlam   rate at which lambda changes
 * \param phi1   minimum value of phi
 * \param phi2   maximum value of phi
 * \param lam1   minimum lambda
 * \param lam2   maximum lambda
 * \param pref   reference potential.
 * \param acc    accuracy in position
 * \param pmin   value of function at minimum
 * \param jammed true if minimum is on a boundary
 */

void Roche::linmin(double q, STAR star, double spin, double cosi, double sini,
                   const Subs::Vec3& p, double& phi, double& lam, double dphi,
                   double dlam, double phi1, double phi2, double lam1, double lam2,
                   double pref, double acc, double& pmin, bool& jammed){

    // Create functions for 1D minimisation.
    Rpot   func(q, star, spin, cosi, sini, p, phi, dphi, lam, dlam);
    Drpot dfunc(q, star, spin, cosi, sini, p, phi, dphi, lam, dlam);

    // Current point equivalent to x=0. Compute maximum before hitting a
    // boundary, and the boundary in question.
    double xmax = 1.e30;
    double xend;
    int nbound = 0;

    if(dphi != 0.){
        xend = (phi1-phi)/dphi;
        if(xend > 0.&& xend < xmax){
            nbound = 1;
            xmax   = xend;
        }
        xend = (phi2-phi)/dphi;
        if(xend > 0. && xend < xmax){
            nbound = 2;
            xmax   = xend;
        }
    }

    if(dlam != 0.){
        xend = (lam1-lam)/dlam;
        if(xend > 0.&& xend < xmax){
            nbound = 3;
            xmax   = xend;
        }
        xend = (lam2-lam)/dlam;
        if(xend > 0.&& xend < xmax){
            nbound = 4;
            xmax   = xend;
        }
    }

    // Now the aim is to bracket the minimum, while accounting for the maximum
    // possible step so that we can then apply dbrent.
    double xa = 0.,         fa = func(xa);
    double xb = 1.e-8*xmax, fb = func(xb);
    int nten = 0;
    const int NTEN = 7;
    while((fb >= fa || xa == xb) && nten < NTEN){
        nten++;
        xb *= 10.;
        fb  = func(xb);
    }
    if(fb <= pref){
        phi += dphi*xb;
        lam += dlam*xb;
        pmin = fb;
        return;
    }
    if(fb >= fa){
        // Let's hope that we have not stepped past the minimum without
        // knowing it
        pmin = fa;
        return;
    }

    // OK, so fb < fa so we are heading downhill at least. Now try
    // to find other side starting from xb, looking for a point when
    // we go up or dip below the critical potential
    bool bracketted = false;
    double xc, fc, dc, xbold = xb, fbold=fb;
    const int NTRY = 5;
    xmax -= xb;
    for(int n=1; n<=NTRY; n++){
        xc = xbold + xmax*n/NTRY;
        fc = func(xc);
        if(fc <= pref){
            phi += dphi*xc;
            lam += dlam*xc;
            pmin = fc;
            return;
        }
        if(fc < fb){
            xb = xc;
            fb = fc;
        }else{
            bracketted = true;
            break;
        }
    }

    jammed = false;

    if(!bracketted){

        dc = dfunc(xc);
        if(dc > 0.){

            // We have crashed into the end stop without crossing the minimum
            // but the derivative says that we are going up. Damn!  Go back to
            // old xb and check that derivative was going down there it really
            // ought to have been ...
            xb = xbold;
            double db = dfunc(xb);
            if(db < 0.){

                // OK, let's try to zero in on the point at which the derivative
                // switches sign.
                double xm = (xb+xc)/2.;
                while((xc-xb) > 1.e-8*xc){
                    xm = (xb+xc)/2.;
                    if(dfunc(xm) > 0.){
                        xc = xm;
                    }else{
                        xb = xm;
                    }
                }
                double fm = func(xm);
                if(fm <= pref){
                    phi += dphi*xm;
                    lam += dlam*xm;
                    pmin = fm;
                    return;
                }
                if(fm < fc && fm < fbold){
                    xb = xm;
                    bracketted = true;
                }else{
                    std::cerr << "Diagnostics" << std::endl;
                    std::cerr << "xb,xc,db,dc= " << xb << "," << xc << "," << db << ","
                              << dc << std::endl;
                    std::cerr << "fbold,fm-fbold,fm-fc= " << fbold << ", " << fm-fbold
                              << ", " << fm-fc << std::endl;
                    std::cerr << "q,ci,si,p = " << q << ", " << cosi << ", " << sini
                              << ", " << p << std::endl;
                    std::cerr << "phi, lam, dphi, dlam = " << phi << ", " << lam
                              << ", " << dphi << ", " << dlam << std::endl;
                    std::cerr << "phi1,phi2,lam1,lam2 = " << phi1 << ", " << phi2
                              << ", " << lam1 << ", " << lam2 << std::endl;
                    std::cerr << "pref, acc = " << pref << ", " << acc << std::endl;
                    for(int i=0; i<=50; i++){
                        xb = xmax*i/50.;
                        std::cerr << i << " " << xb << " " << func(xb) << " "
                                  << dfunc(xb) << std::endl;
                    }
                    throw Roche_Error("Roche::linmin: failed to bracket minimum, error 3");
                }

            }else{

                std::cerr << "Diagnostics" << std::endl;
                std::cerr << "xb,xc,db,dc= " << xb << "," << xc << "," << db
                          << "," << dc << std::endl;
                std::cerr << "q,ci,si,p = " << q << ", " << cosi << ", " << sini
                          << ", " << p << std::endl;
                std::cerr << "phi, lam, dphi, dlam = " << phi << ", " << lam
                          << ", " << dphi << ", " << dlam << std::endl;
                std::cerr << "phi1,phi2,lam1,lam2 = " << phi1 << ", " << phi2
                          << ", " << lam1 << ", " << lam2 << std::endl;
                std::cerr << "pref, acc = " << pref << ", " << acc << std::endl;
                for(int i=0; i<=50; i++){
                    xb = xmax*i/50.;
                    std::cerr << i << " " << xb << " " << func(xb) << " "
                              << dfunc(xb) << std::endl;
                }
                throw Roche_Error("Roche::linmin: failed to bracket minimum, error 1");
            }

        }else{

            // We are trapped on a boundary; re-define line minimisation functions.
            jammed = true;
            phi  += dphi*xmax;
            lam  += dlam*xmax;
            xmax = 1.;
            rpot_grad(q, star, spin, set_earth(cosi, sini, phi), p, lam, dphi, dlam);
            switch(nbound){
            case 1:
                phi  = phi1;
                dphi = 0.;
                if(dlam > 0.){
                    dlam = lam1 - lam;
                }else{
                    dlam = lam2 - lam;
                }
                break;
            case 2:
                phi  = phi2;
                dphi = 0.;
                if(dlam > 0.){
                    dlam = lam1 - lam;
                }else{
                    dlam = lam2 - lam;
                }
                break;
            case 3:
                lam  = lam1;
                dlam = 0.;
                if(dphi > 0.){
                    dphi = phi1 - phi;
                }else{
                    dphi = phi2 - phi;
                }
                break;
            case 4:
                lam  = lam2;
                dlam = 0.;
                if(dphi > 0.){
                    dphi = phi1 - phi;
                }else{
                    dphi = phi2 - phi;
                }
            }
            func.reset(q, star, spin, cosi, sini, p,  phi, dphi, lam, dlam);
            dfunc.reset(q, star, spin, cosi, sini, p, phi, dphi, lam, dlam);

            // Again try to bracket minimum
            nten = 0;
            xb   = 1.e-6;
            fb   = func(xb);
            while(fb >= fa && nten < NTEN){
                xb *= 10.;
                nten++;
                fb = func(xb);
            }
            if(fb <= pref){
                phi += dphi*xb;
                lam += dlam*xb;
                pmin = fb;
                return;
            }
            if(fb >= fa){
                pmin = fa;
                return;
            }

            // Now the bracketting steps.
            bracketted = false;
            for(int n=1; n<=NTRY; n++){
                xc = xmax*n/NTRY;
                fc = func(xc);
                if(fc <= pref){
                    phi += dphi*xc;
                    lam += dlam*xc;
                    pmin = fc;
                    return;
                }
                if(fc < fb){
                    xb = xc;
                    fb = fc;
                }else{
                    bracketted = true;
                    break;
                }
            }
            if(!bracketted){

                if(dfunc(xc) > 0.){

                    std::cerr << "Diagnostics" << std::endl;
                    std::cerr << "xb,xc,fb,fc= " << xb << "," << xc << ","
                              << "," << fb << "," << fc << std::endl;
                    std::cerr << "q,ci,si,p = " << q << ", " << cosi << ", " << sini
                              << ", " << p << std::endl;
                    std::cerr << "phi, lam, dphi, dlam = " << phi << ", " << lam
                              << ", " << dphi << ", " << dlam << std::endl;
                    std::cerr << "phi1,phi2,lam1,lam2 = " << phi1 << ", " << phi2
                              << ", " << lam1 << ", " << lam2 << std::endl;
                    std::cerr << "pref, acc = " << pref << ", " << acc << std::endl;
                    for(int i=0; i<=50; i++){
                        xb = xmax*i/50.;
                        std::cerr << i << " " << xb << " " << func(xb) << " "
                                  << dfunc(xb) << std::endl;
                    }

                    throw Roche_Error("Roche::linmin failed to bracket minimum, error 2");
                }else{
                    phi  += dphi*xmax;
                    lam  += dlam*xmax;
                    pmin = fc;
                    return;
                }
            }
        }
    }

    // OK, finally, minimum bracketted, so refine it.
    double xmin;
    double xacc = acc/sqrt(Subs::sqr(Constants::TWOPI*dphi)+Subs::sqr(dlam));
    pmin = dbrent(xa, xb, xc, func, dfunc, xacc, true, pref, xmin);

    phi += dphi*xmin;
    lam += dlam*xmin;

}

