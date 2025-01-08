#include <cstdlib>
#include <cmath>
#include <iostream>
#include <string>
#include "trm/subs.h"
#include "trm/constants.h"
#include "trm/vec3.h"
#include "trm/roche.h"

using Subs::operator+;

/**
 * fourth_contact computes the fourth contact phase of a white dwarf being eclipsed by
 * a Roche-lobe filling secondary star. More correctly it computes last contact if there
 * are no second and third contacts.
 *
 * \param q       the mass ratio = M2/M1.
 * \param iangle  the orbital inclination in degrees
 * \param rw      the radius of the white dwarf scaled by a.
 * \param ntheta  the number of points to compute over the loer-right quadrant of the white dwarf to test for eclipse.
 * \param acc     step size parameter for 'blink'.
 * \param contact the contact phase 
 * \return true/false as to whether there was a contact phase or not
 */


bool wd_is_eclipsed_at_all(double q, double iangle, double rw, int ntheta, double acc, double phase);

bool Roche::fourth_contact(double q, double iangle, double rw, int ntheta, double acc, double& contact){

  double p1 = 0., p2 = 0.25;

  if(!wd_is_eclipsed_at_all(q,iangle,rw,ntheta,acc,p1)) return false;

  while(p2-p1 > 1.e-9){
    contact = (p1+p2)/2.;
    if(wd_is_eclipsed_at_all(q,iangle,rw,ntheta,acc,contact)){
      p1 = contact;
    }else{
      p2 = contact;
    }
  }
  return true;
}

bool wd_is_eclipsed_at_all(double q, double iangle, double rw, int ntheta, double acc, double phase){

  double cosi = cos(Constants::TWOPI*iangle/360.);
  double sini = sin(Constants::TWOPI*iangle/360.);

  Subs::Vec3 earth, x, y, r;
  double cosp = cos(Constants::TWOPI*phase);
  double sinp = sin(Constants::TWOPI*phase);
  
  earth.x() =  sini*cosp;
  earth.y() = -sini*sinp;
  earth.z() = cosi;
  
  // Compute two vectors spanning plane in sky that cuts through centre of white dwarf
  x.x() = sinp;
  x.y() = cosp;
  x.z() = 0.;
  
  y.x() = earth.y()*x.z() - earth.z()*x.y();
  y.y() = earth.z()*x.x() - earth.x()*x.z();
  y.z() = earth.x()*x.y() - earth.y()*x.x();
  
  double theta, cost, sint;
  for(int i=0; i<ntheta; i++){
    theta = Constants::TWOPI*(0.75 + 0.25*i/(ntheta-1));
    cost  = cos(theta);
    sint  = sin(theta);
    
    // Compute position vector of point on edge of white dwarf
    r.x()   = rw*(cost*x.x() + sint*y.x());
    r.y()   = rw*(cost*x.y() + sint*y.y());
    r.z()   = rw*(cost*x.z() + sint*y.z());
    
    // Test for eclipse
    if(Roche::blink(q,r,earth,acc)) return true;
  }
  return false;
}

bool wd_is_eclipsed_totally(double q, double iangle, double rw, int ntheta, double acc, double phase);

/**
 * third_contact computes the third contact phase of a white dwarf being eclipsed by
 * a Roche-lobe filling secondary star if there is one. 
 *
 * \param q       the mass ratio = M2/M1.
 * \param iangle  the orbital inclination in degrees
 * \param rw      the radius of the white dwarf scaled by a.
 * \param ntheta  the number of points to compute over the loer-right quadrant of the white dwarf to test for eclipse.
 * \param acc     step size parameter for 'blink'.
 * \param contact the contact phase 
 * \return true/false as to whether there was a contact phase or not
 */

bool Roche::third_contact(double q, double iangle, double rw, int ntheta, double acc, double& contact){

  double p1 = 0., p2 = 0.25;

  if(!wd_is_eclipsed_totally(q,iangle,rw,ntheta,acc,p1)) return false;

  while(p2-p1 > 1.e-9){
    contact = (p1+p2)/2.;
    if(wd_is_eclipsed_totally(q,iangle,rw,ntheta,acc,contact)){
      p1 = contact;
    }else{
      p2 = contact;
    }
  }
  return true;
}

bool wd_is_eclipsed_totally(double q, double iangle, double rw, int ntheta, double acc, double phase){

  double cosi = cos(Constants::TWOPI*iangle/360.);
  double sini = sin(Constants::TWOPI*iangle/360.);

  Subs::Vec3 earth, x, y, r;
  double cosp = cos(Constants::TWOPI*phase);
  double sinp = sin(Constants::TWOPI*phase);
  
  earth.x() =  sini*cosp;
  earth.y() = -sini*sinp;
  earth.z() = cosi;
  
  // Compute two vectors spanning plane in sky that cuts through centre of white dwarf
  x.x() = sinp;
  x.y() = cosp;
  x.z() = 0.;
  
  y.x() = earth.y()*x.z() - earth.z()*x.y();
  y.y() = earth.z()*x.x() - earth.x()*x.z();
  y.z() = earth.x()*x.y() - earth.y()*x.x();
  
  double theta, cost, sint;
  for(int i=0; i<ntheta; i++){
    theta = Constants::TWOPI*(0.25 + 0.25*i/(ntheta-1));
    cost  = cos(theta);
    sint  = sin(theta);
    
    // Compute position vector of point on edge of white dwarf
    r.x()   = rw*(cost*x.x() + sint*y.x());
    r.y()   = rw*(cost*x.y() + sint*y.y());
    r.z()   = rw*(cost*x.z() + sint*y.z());
    
    // Test for eclipse
    if(!Roche::blink(q,r,earth,acc)) return false;
  }
  return true;
}

