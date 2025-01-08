#include <cstdlib>
#include <cmath>
#include <iostream>
#include <string>
#include "trm/subs.h"
#include "trm/constants.h"
#include "trm/roche.h"

// The routine here works by considering the 'line-of-sight cone' (LOSC) formed by the line of sight
// towards a point when considering all possible phases. Since the disc is modelled as cylindrically
// symmetric with the z-axis as the axis, then at various times one needs to know how the LOSC
// intersects with circles with parallel axes. The way this is worked out is by considering the
// the interesection of the LOSC with the plane containing the circle of interest. The intersection
// is itself a circle, so that it boils down to considering the intersection of two circles.

namespace Roche {

  // This enumerates the 5 possible outcomes of the LOSC intersection with a circle.
  enum CIRCLE {

    // Line of sight cone starts at or above the circle of interest
    ABOVE,

    // Line of sight circle is everywhere inside circle of interest
    INSIDE,    
    
    // Line of sight circle is everywhere outside circle of interest
    OUTSIDE,   

    // Line of sight circle is separated from the circle of interest
    SEPARATE,
  
    // Line of sight circle cone intersects the circle of interest
    CROSSING   
  };
};

double cut_phase(double rxy, double rcone, double radius);

Roche::CIRCLE circle_eclipse(double rxy, double z, double zcirc, double radius, double tani, double& phase);  

/**
 * disc_eclipse works out phase ranges during which a cylindrically symmetric, flared disc 
 * running between a pair of radii eclipses a given point. 
 *
 * \param iangle  the orbital inclination, degrees. 90 = edge on.
 * \param r       the position vector of the point in question (units of binary separation)
 * \param rdisc1 inner disc  radius, units of separation
 * \param rdisc2 outer disc radius, units of separation
 * \param beta exponent of flaring, so that the height scales as r**beta. beta should be >= 1
 * \param height disc height at unit radius in disc (even if it does not exist)
 * \return a vector of ingress and egress phase pairs during which the point in question is eclipsed.
 * The ingress phase will always be between 0 and 1 while the egress phase will be larger than this, but
 * by no more than 1 cycle. If the vector is null, no eclipse takes place.
 */

std::vector<std::pair<double,double> > Roche::disc_eclipse(double iangle, double rdisc1, double rdisc2, double beta, double height, const Subs::Vec3& r){

  if(beta < 1) throw Roche_Error("Roche::disc_eclipse: beta must be >= 1");
  
  // Compute and store cosine and sine of inclination if need be.
  static double iangle_old = -1.e30, sini, cosi;
  if(iangle != iangle_old){
    iangle_old = iangle;
    sini = sin(Constants::PI/180.*iangle);
    cosi = cos(Constants::PI/180.*iangle);
  }
  
  std::vector<std::pair<double,double> > temp;
  
  // Compute height of disc at outer boundary
  const double HOUT = height*pow(rdisc2, beta);
  
  // Deal with points too high ever to be eclipsed whatever the inclination
  if(r.z() >= HOUT) return temp;
  
  // Special case of exactly edge-on, only curved outer edge matters.
  if(cosi == 0.){
    if(fabs(r.z()) < HOUT){
      double rxy = sqrt(Subs::sqr(r.x()) + Subs::sqr(r.y()));
      if(rxy <= rdisc2){
	temp.push_back(std::make_pair(0.,1.));
      }else{
	double subtend = asin(rdisc2/rxy)/Constants::TWOPI;
	double pcen    = atan2(r.y(), -r.x())/Constants::TWOPI;
	double ingress = pcen - subtend;
	ingress = ingress - floor(ingress);
	double egress  = ingress + 2.*subtend;
	temp.push_back(std::make_pair(ingress,egress));
      }
    }
    return temp;
  }
  
  // Work out distance from axis
  const double RXY   = sqrt(Subs::sqr(r.x()) + Subs::sqr(r.y()));
  
  if(rdisc1 < RXY && RXY < rdisc2 && fabs(r.z()) < height*pow(RXY, beta)){
    // Point is inside disc and so is eclipsed
    temp.push_back(std::make_pair(0., 1.1));
    return temp;
  }
  
  const double TANI = sini/cosi;  
  Roche::CIRCLE result;
  double phase, ingress, egress;
  
  if(RXY < rdisc2 && r.z() >= height*pow(std::max(rdisc1, RXY), beta)){

    // Point is in approximately conical region above the disc. Just need to check whether
    // it is not occulted by the edge of the disc
    result = circle_eclipse(RXY, r.z(), HOUT, rdisc2, TANI, phase);
    
    if(result == Roche::OUTSIDE){
      // point will be occulted by the disc edge at all phases
      temp.push_back(std::make_pair(0., 1.1));
    }else if(result == Roche::CROSSING){
      // point partially occulted by disc edge; work out phases
      double phi0 = atan2(r.y(), -r.x())/Constants::TWOPI;
      ingress     = phi0 + phase;
      ingress    -= floor(ingress);
      egress      = ingress + 1 - 2.*phase;
      temp.push_back(std::make_pair(ingress, egress));
    }
    return temp;
  }

  // Compute the radius of circle formed by LOSC in the plane of 
  // the lower outer rim of the disc
  const double RCONE_LO = std::max(0., TANI*(-HOUT-r.z()));

  // Circle encloses rim, so no intersection
  if(RCONE_LO >= RXY + rdisc2) return temp;

  // Compute the radius of circle formed by LOSC in the plane of 
  // the upper outer rim of the disc
  const double RCONE_HI = TANI*(HOUT-r.z());

  // Circle disjoint from rim, so no intersection
  if(RXY >= RCONE_HI + rdisc2) return temp;

  // For the moment we pretend that the disc has no hole at its centre, so
  // that we are simply interested in the phases over which eclipse occurs. 
  // At this point we are guaranteed that this will happen. All events are
  // symmetrically located around a phase defined by x and y only which will
  // be calculated at the end. We therefore just find the half range which
  // is called 'eclipse_phase' below.

  double eclipse_phase = 0;
  if(RXY + RCONE_LO <= rdisc2){

    // Cone swept out by line of sight always inside lower face so total eclipse
    eclipse_phase = 0.5;

  }else if(RXY <= rdisc2){

    // Points that project close to the z axis which are only 
    // partially obscured by the disc hovering above them.
    // this means they must be below -HOUT    
    eclipse_phase = cut_phase(RXY, RCONE_LO, rdisc2);

  }else{

    // Points further from the z axis than the outer rim of the disc that will be eclipsed.     
 
    if(Subs::sqr(RCONE_HI) + Subs::sqr(rdisc2) >= Subs::sqr(RXY) &&
       Subs::sqr(RCONE_LO) + Subs::sqr(rdisc2) <= Subs::sqr(RXY)){
     // In this case it is the curved outer disc rim that sets the limit
      eclipse_phase = asin(rdisc2/RXY)/Constants::TWOPI;

    }else if(Subs::sqr(RCONE_HI) + Subs::sqr(rdisc2) < Subs::sqr(RXY)){
      // In this case it is upper outer rim that sets the limit
      eclipse_phase = cut_phase(RXY, RCONE_HI, rdisc2);

    }else{
      // In this case it is lower outer rim that sets the limit
      eclipse_phase = cut_phase(RXY, RCONE_LO, rdisc2);
    }
  }

  // At this point we have covered all cases for the eclipse, whilst ignoring the
  // possibility of seeing the point through the hole in the middle of the disc.
  // Now let's calculate the 'appear_phase' if any.

  // First compute height of disc at inner boundary
  const double HIN = height*pow(rdisc1, beta);

  double appear_phase = -1;

  if(r.z() < -HOUT){
    // In this case the LOSC has to run through 4 circles which are the upper and
    // lower outer and inner rims.

    // First, the lower outer rim
    result = circle_eclipse(RXY, r.z(), -HOUT, rdisc2, TANI, phase);
    if(result == Roche::INSIDE){
      appear_phase = 0.5;
    }else if(result == Roche::CROSSING){
      appear_phase = phase;
    }

    // Second, the lower inner rim
    if(appear_phase > 0){
      result = circle_eclipse(RXY, r.z(), -HIN, rdisc1, TANI, phase);
      if(result == Roche::CROSSING){
	appear_phase = std::min(appear_phase, phase);
      }else if(result != Roche::INSIDE){
	appear_phase = -1.;
      }
    }

    // Third, the upper inner rim
    if(appear_phase > 0){
      result = circle_eclipse(RXY, r.z(), HIN, rdisc1, TANI, phase);
      if(result == Roche::CROSSING){
	appear_phase = std::min(appear_phase, phase);
      }else if(result != Roche::INSIDE){
	appear_phase = -1.;
      }
    }

    // Fourth, the upper outer rim
    if(appear_phase > 0){
      result = circle_eclipse(RXY, r.z(), HOUT, rdisc2, TANI, phase);
      if(result == Roche::CROSSING){
	appear_phase = std::min(appear_phase, phase);
      }else if(result != Roche::INSIDE){
	appear_phase = -1.;
      }
    }

  }else if(RXY < rdisc1){

    if(r.z() < -HIN){

      // Points hovering around underside of disc. Have to consider just three circles
      
      // First, the lower inner rim
      result = circle_eclipse(RXY, r.z(), -HIN, rdisc1, TANI, phase);
      if(result == Roche::INSIDE){
	appear_phase = 0.5;
      }else if(result == Roche::CROSSING){
	appear_phase = phase;
      }
      
      // Second, the upper inner rim
      if(appear_phase > 0){
	result = circle_eclipse(RXY, r.z(), HIN, rdisc1, TANI, phase);
	if(result == Roche::CROSSING){
	  appear_phase = std::min(appear_phase, phase);
	}else if(result != Roche::INSIDE){
	  appear_phase = -1.;
	}
      }
      
      // Third, the upper outer rim
      if(appear_phase > 0){
	result = circle_eclipse(RXY, r.z(), HOUT, rdisc2, TANI, phase);
	if(result == Roche::CROSSING){
	  appear_phase = std::min(appear_phase, phase);
	}else if(result != Roche::INSIDE){
	  appear_phase = -1.;
	}
      }

    }else if(r.z() < HIN){

      // Points inside hole in middle of disc. Have to consider just two circles
      
      // First, the upper inner rim
      result = circle_eclipse(RXY, r.z(), HIN, rdisc1, TANI, phase);
      if(result == Roche::INSIDE){
	appear_phase = 0.5;
      }else if(result == Roche::CROSSING){
	appear_phase = phase;
      }
      
      // Second, the upper outer rim
      if(appear_phase > 0){
	result = circle_eclipse(RXY, r.z(), HOUT, rdisc2, TANI, phase);
	if(result == Roche::CROSSING){
	  appear_phase = std::min(appear_phase, phase);
	}else if(result != Roche::INSIDE){
	  appear_phase = -1.;
	}
      }
    }
  }

  // Here is the central phase
  double phi0 = atan2(r.y(), -r.x())/Constants::TWOPI;

  if(appear_phase <= 0){
    ingress     = phi0 - eclipse_phase;
    ingress    -= floor(ingress);
    egress      = ingress + 2.*eclipse_phase;
    temp.push_back(std::make_pair(ingress, egress));
  }else if(appear_phase < eclipse_phase){
    ingress     = phi0 - eclipse_phase;
    ingress    -= floor(ingress);
    egress      = ingress + (eclipse_phase-appear_phase);
    temp.push_back(std::make_pair(ingress, egress));
    ingress     = phi0 + appear_phase;
    ingress    -= floor(ingress);
    egress      = ingress + (eclipse_phase-appear_phase);
    temp.push_back(std::make_pair(ingress, egress));
  }
  return temp;
}

/** Determines the nature of intersection of the LOSC cone to a point with a circle
 * centred on the z axis. 
 * \param rxy the distance of the point from the z axis
 * \param z   the z coordinate of the point
 * \param zcirc the z ordinate of the circle 
 * \param rcirc the radius of the circle
 * \param tani tangent of the inclination
 * \param result the phase at which any crossing occurs
 * \return the outcome of the intersection computation, 5 possibilities. See trm_roche.h
 */

Roche::CIRCLE circle_eclipse(double rxy, double z, double zcirc, double radius, double tani, double& phase){

  double cut_phase(double rxy, double rcone, double radius);

  // point above circle
  if(z >= zcirc) return Roche::ABOVE;

  const double RCONE = tani*(zcirc - z);
  
  // line-of-sight always outside the circle
  if(RCONE >= rxy + radius) return Roche::OUTSIDE;

  // line-of-sight circle separate from the circle
  if(rxy >= RCONE + radius) return Roche::SEPARATE;

  // line-of-sight always outside the circle
  if(rxy + RCONE <= radius) return Roche::INSIDE;

  // crossing case
  phase = cut_phase(rxy, RCONE, radius);

  return Roche::CROSSING;
  
}

// Helper function that is often needed. It calculates a phase corresponding to the
// point at which a circle of radius = 'radius' intersects one of radius = 'rcone'
// when their centres are separated by a distance = 'rxy'. 
// For it to work the following conditions must hold:
// 1) rxy + rcone > radius -- to guarantee that the rcone circle is not wholly inside the radius circle
// 2) rxy < radius + rcone -- to guarantee that the two circles are not completely separate
// 3) rcone < radius + rxy -- to guarantee that the rcone circle does not wholly enclose the radius circle
// to ensure that an intersection occurs. These will not be checked.

double cut_phase(double rxy, double rcone, double radius){

  // Temporary checks
  if(rxy + rcone <= radius)
    throw Roche::Roche_Error("rxy + rcone <= radius");
  if(rxy >= radius + rcone)
    throw Roche::Roche_Error("rxy >= radius + rcone");
  if(rcone >= radius + rxy)
    throw Roche::Roche_Error("rcone >= radius + rxy");

  return acos((Subs::sqr(rxy) + Subs::sqr(rcone) - Subs::sqr(radius))/(2.*rcone*rxy))/Constants::TWOPI;
}
