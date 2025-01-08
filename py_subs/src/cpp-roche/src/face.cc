#include <cstdlib>
#include <cmath>
#include <iostream>
#include <string>
#include "trm/subs.h"
#include "trm/vec3.h"
#include "trm/roche.h"

/**
 * 'face' computes the position and orientation of a face on either star in a binary assuming Roche geometry given
 * a direction, a reference radius and a potential.
 *
 * \param q    the mass ratio = M2/M1.
 * \param star specifies which star, primary or secondary is under consideration.
 * \param spin ratio of star in questions spin to the orbital frequency
 * \param dirn the direction (unit) vector from the centre of mass of the secondary to the face in question.
 * \param rref reference radius. This is a radius large enough to guarantee crossing of the reference potential. See ref_sphere
 * \param pref reference potential. This defines the precise location of the face.
 * \param acc  location accuracy (units of separation)
 * \param pvec position vector of centre of face (position vector in standard binary coordinates), returned
 * \param dvec orientation vector perpendicular to face, returned
 * \param r    distance from centre of mass of star, returned
 * \param g    magnitude of gravity at face, returned
 * \exception The routine throws exceptions if it cannot bracket the reference potential. This can occur if the reference radius fails to enclose
 * the face in question, or if the face is so deep in the potential that the initial search fails to reach it. Finally if acc is set too low an
 * exception may be thrown if too many binary chops occur. The behaviour at the L1 point is undefined so do not try to call it there.
 */

void Roche::face(double q, STAR star, double spin, const Subs::Vec3& dirn, double rref, double pref, double acc,
                 Subs::Vec3& pvec, Subs::Vec3& dvec, double& r, double& g){

    // centre of mass in question
    const Subs::Vec3 cofm = (star == PRIMARY) ? Subs::Vec3(0.,0.,0.) : Subs::Vec3(1.,0.,0.);

    // Pointer to Roche potential and potential derivative functions
    double (*rp)(double q, double spin, const Subs::Vec3& p) = (star == PRIMARY) ? &rpot1 : &rpot2;
    Subs::Vec3 (*drp)(double q, double spin, const Subs::Vec3& p) = (star == PRIMARY) ? &drpot1 : &drpot2;

    // A check on the reference radius & potential
    double tref = rp(q, spin, cofm + rref*dirn);
    if(tref < pref)
        throw Roche_Error("Roche::face error: point at reference radius = " + Subs::str(rref) +
                          " appears to have a lower potential = " + Subs::str(tref) +
                          " than the reference = " + Subs::str(pref) +
                          "\nOther params: q = " + Subs::str(q) +
                          ", direction = (" + Subs::str(dirn.x()) +
                          "," + Subs::str(dirn.y()) + "," +
                          Subs::str(dirn.z()) + ")" );

    // Find r1 r2 such that r1 is below reference potential and r2 is above.
    double r1 = rref/2, r2 = rref;
    tref = pref + 1;

    const int MAXSEARCH = 30;
    for(int i=0; i<MAXSEARCH && tref > pref; i++){
        r1   = r2/2;
        tref = rp(q, spin, cofm + r1*dirn);
        if(tref > pref)
            r2 = r1;
    }
    if(tref > pref)
        throw Roche_Error("Roche::face error: could not find a radius with a potential below the reference potential; probably bad inputs.");

    // OK now refine with a binary chop. Crude but robust.
    const int MAXCHOP = 100;
    int nchop = 0;
    while(r2-r1 > acc && nchop < MAXCHOP){
        r = (r1+r2)/2.;
        pvec = cofm + r*dirn;
        if(rp(q, spin, pvec) < pref){
            r1 = r;
        }else{
            r2 = r;
        }
        nchop++;
    }
    if(nchop == MAXCHOP)
        throw Roche_Error("Roche::face error: reached maximum number of binary chops = " + Subs::str(MAXCHOP) );

    r     = (r1+r2)/2.;
    pvec  = cofm + r*dirn;
    dvec  = drp(q, spin, pvec);
    g     = dvec.length();
    dvec /= g;
}




