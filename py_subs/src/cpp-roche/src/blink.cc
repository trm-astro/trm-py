#include <cstdlib>
#include <cmath>
#include <iostream>
#include <string>
#include "trm/subs.h"
#include "trm/vec3.h"
#include "trm/roche.h"

/**
 * blink tests for the occultation of a point in a semi-detached
 * binary. blink = true when point is eclipsed, otherwise = false. It works 
 * by first trying to eliminate as many paths that go nowhere near the 
 * Roche lobe before getting down to the hard work of stepping
 * along the path.
 *
 * NB Note that this version of !!emph{blink} works in units of binary separation
 * NOT the inner lagrangian point distance.
 *
 * \param q the mass ratio = M2/M1.
 * \param r the position vector (in units of binary separation)
 * \param e unit vector pointing towards Earth. Standardly this is (sini(i)*cos(phi),-sin(i)*sin(phi),cos(i))
 * \param acc step size parameter. acc specifies the size of steps taken when trying to see if the 
 * photon path goes inside the Roche lobe. The step size is roughly acc times the radius of the lobe filling star. 
 * This means that the photon path could mistakenly not be occulted if it passed less than about (acc**2)/8 of 
 * the radius below the surface of the Roche lobe. acc of order 0.1 should therefore do the job.
 * \return false = not eclipsed; true = eclipsed.
 */

bool Roche::blink(double q, const Subs::Vec3& r, const Subs::Vec3& e, double acc){

  int i, nstep;
  double rl1, rsphere, par1, par2, fac, a, b, c, x1, x2, x3;
  double yy, xc, r2, p, p1, p2, dp, r1, par, rr, rs1, rs2;
  double xm, cmax, xt;

  // For speed several variables which only need computing if q is changed
  // are preserved between invocations

  static double qq = -100., xcm, c1, c2, step, pp, crit;

  // Compute q dependent quantities

  if(q != qq){
    if(q <= 0.)
      throw Roche::Roche_Error("Invalid mass ratio (=" + Subs::str(q) + ") in blink.");

    if(acc <= 0.)
      throw Roche::Roche_Error("Invalid accuracy parameter (=" + Subs::str(acc) + ") in blink.");

    qq   = q;
    xcm  = 1./(1.+q);
    c1   = 2.*xcm;
    xcm *= q;
    c2   = 2.*xcm;

    // Locate the inner Lagrangian point (L1)

    rl1 = xl1(q);

    // Evaluate Roche potential at L1 point.
    
    r1   = rl1;
    r2   = 1.-rl1;
    xc   = rl1 - xcm;
    crit = c1/r1 + c2/r2 + xc*xc;

    // The red star lies entirely within the sphere centred on its
    // centre of mass and reaching the inner Lagrangian point.
    
    rsphere = 1.-rl1;
    pp      = rsphere*rsphere;
    step    = rsphere*acc;
  }

  // From now on computations are done every call. Main point is
  // to try to bail out as soon as possible to save time 

  // evaluate closest approach distance to sphere.

  b = e.x()*(xt = r.x()-1.) + e.y()*r.y() + e.z()*r.z();
  c = xt*xt + r.y()*r.y() + r.z()*r.z() - pp;

  // Photon path crosses sphere at two points given by quadratic
  // equation l*l + 2*b*l + c = 0. First check that this has
  // real roots. If not, there is no eclipse. This test should
  // eliminate most cases. 

  if((fac = b*b - c) <= 0.) return false;

  fac  = sqrt(fac);

  // If the larger root is negative, the photon starts after the sphere.
  // This should get rid of another whole stack

  if((par2 = -b + fac) <= 0.) return false;

  //  Now for the hard work. The photon's path does go
  //  inside the sphere ... groan. First evaluate smaller root 
  //  but limit to >= 0 to ensure that only photon path
  //  and not its extrapolation backwards in included.

  par1 =  - b - fac;
  par1 = par1 > 0.? par1 : 0.;

  // Now follow the photon's path in finite steps to see if 
  // it intercepts the red star. We start off at closest
  // approach point to increase chance of only needing one 
  // computation of the roche potential

  par =  b < 0. ? - b : 0.;
  x1  = r.x() + par*e.x();
  x2  = r.y() + par*e.y();
  x3  = r.z() + par*e.z();

  // Test roche potential for an occultation 
  xm = x1 - 1.;
  rs2 = xm*xm + (rr = (yy = x2*x2) + x3*x3);

  // Point at c of m of red star. Definitely eclipsed, avoids
  // division by 0 later

  if(rs2 <= 0.) return true;

  r1 = sqrt(rs1 = x1*x1 + rr);
  r2 = sqrt(rs2);

  // Deeper in well than inner Lagrangian, therefore eclipsed.

  xc = x1 - xcm;
  if((c  = c1/r1 + c2/r2 + xc*xc + yy) > crit) return true;

  // Now we need to step. Determine step direction by 
  // evaluating the first derivative

  a = x1*e.x() + x2*e.y();
  b = a + x3*e.z();
  if(-c1*b/(rs1*r1)-c2*(b-e.x())/(rs2*r2)+2.*(a-xcm*e.x()) > 0.){
    p1 = par;
    p2 = par2;
  }else{
    p1 = par;
    p2 = par1;
  }

  // Loop while the Roche potential increases in depth

  nstep = int(floor(fabs(p2-p1)/step+0.5));
  nstep = nstep > 2 ? nstep : 2;
  dp    = (p2-p1)/nstep;
  cmax  = c - 1.;
  i     = 0;
  while(c > cmax && i < nstep){
    i++;
    cmax = c;
    p  = p1 + dp*i;
    x1 = r.x() + p*e.x();
    x2 = r.y() + p*e.y();
    x3 = r.z() + p*e.z();

    // test roche potential for an occultation, again guarding
    // against division by 0
    
    xm = x1 - 1.;
    r2 = sqrt(xm*xm + (rr = (yy = x2*x2) + x3*x3));
    if(r2 <= 0.) return true;
    r1 = sqrt(x1*x1 + rr); 
    xc = x1 - xcm;
    if((c =  c1/r1 + c2/r2 + xc*xc + yy) > crit) return true;
  }
  return false;
}




