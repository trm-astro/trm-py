#include "trm/subs.h"
#include "trm/vec3.h"
#include "trm/roche.h"

/**
 * Works out whether the donor is irradiated by a given place on the
 * equator of the accretor.
 * \param q mass ratio = M2/M1
 * \param x x position of place on accretor
 * \param y y position of place on accretor
 * \return true if any of donor can see the point.
 */ 
  
bool Roche::irrad(double q, double x, double y){
  if(x >= 0.) return true;
  Subs::Vec3 r, e;
  r.x() = x;
  r.y() = y;
  r.z() = 0.;
  double d = r.length();
  e.x() =  y/d;
  e.y() = -x/d;
  e.z() = 0.;
  return blink(q,r,e,0.005);
}
