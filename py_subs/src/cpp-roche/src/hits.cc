#include "trm/subs.h"
#include "trm/constants.h"
#include "trm/binary.h"
#include "trm/roche.h"

/**
 * 'hits' works out whether the stream crosses a Roche equipotential.
 * If it does, it returns the position of the first point at which the crossing occurs.
 * If it does not, it returns the position of the point of minimum potential.
 * 
 * \param q the mass ratio = m2/m1
 * \param pot Roche potential. This must be more negative than the critical potential.
 * \param x returned x position of crossing point
 * \param y returned y position of crossing point
 * \return Returns true or false according to whether crossing occured occurred.
 */
  
bool Roche::hits(double q, double pot, double& x, double& y){

  const double EPS   = 1.e-9;   // Integration accuracy
  const double BCHOP = 1.e-8;   // Accuracy of location of critical points in time (1 orbit=2*Pi)

  // Initialise stream
  double ttry = 2.e-4, tdid, time = 0., tnext, told = 0;
  Subs::Vec3 r, v, rold, vold;
  
  strinit(q,r,v);

  // Set first two points
  rold = r;

  gsint(q, r, v, ttry, tdid, tnext, time, EPS);

  // integrate until distance from primary starts increasing
  while(rdot(r,v) < 0.){
    ttry = tnext;
    rold = r;
    vold = v;
    told = time;
    gsint(q, r, v, ttry, tdid, tnext, time, EPS);
  }

  // now binary chop to get closest point
  
  double tlo = 0, thi = tdid;
  while(thi-tlo > BCHOP){
    ttry = (tlo+thi)/2.;
    r    = rold;
    v    = vold;
    time = told;
    gsint(q, r, v, ttry, tdid, tnext, time, EPS);
    if(rdot(r,v) < 0.){
      tlo = ttry;
    }else{
      thi = ttry;
    }
  }
  if(rpot(q,r) >= pot){
    x = r.x();
    y = r.y();
    return false;
  }

  // Now work out where it hits
  double tmin = time;

  strinit(q,r,v);
  time  = 0.;
  tnext = 2.e-4;
  while(rpot(q,r) > pot){
    rold = r;
    vold = v;
    told = time;
    ttry = std::min(tnext,tmin-time);
    gsint(q, r, v, ttry, tdid, tnext, time, EPS);
  }

  // now binary chop again
  tlo = 0., thi = tdid;
  while(thi-tlo > BCHOP){
    ttry = (tlo+thi)/2.;
    r    = rold;
    v    = vold;
    time = told;
    gsint(q, r, v, ttry, tdid, tnext, time, EPS);
    if(rpot(q,r) > pot){
      tlo = ttry;
    }else{
      thi = ttry;
    }
  }
  x = r.x();
  y = r.y();
  return true;
}


/**
 * This version of 'hits' works out whether the stream crosses a Roche equipotential.
 * If it does, it returns the position and velocity of the first point at which the 
 * crossing occurs. If it does not it comes back with the position and velocity of the
 * minimum potential reached.
 * 
 * \param q the mass ratio = m2/m1
 * \param pot Roche potential. This must be more negative than the critical potential.
 * \param x  returned x position of crossing or minimum point (units of separation)
 * \param y  returned y position of crossing or minimum point (units of separation)
 * \param vx returned vx velocity of the crossing or minimum point (units of 2*Pi*separation/P, rotating frame) 
 * \param vy returned vy velocity of the crossing or minimum point (units of 2*Pi*separation/P, rotating frame)
 * \return Returns true or false according to whether crossing occured occurred.
 */
  
bool Roche::hits(double q, double pot, double& x, double& y, double& vx, double& vy){

  const double EPS   = 1.e-9;   // Integration accuracy
  const double BCHOP = 1.e-8;   // Accuracy of location of critical points in time (1 orbit=2*Pi)

  // Initialise stream
  double ttry = 2.e-4, tdid, time = 0., tnext, told = 0;
  Subs::Vec3 r, v, rold, vold;
  
  strinit(q,r,v);

  // Set first two points
  rold = r;

  gsint(q, r, v, ttry, tdid, tnext, time, EPS);

  // Integrate until distance from primary starts increasing
  while(rdot(r,v) < 0.){
    ttry = tnext;
    rold = r;
    vold = v;
    told = time;
    gsint(q, r, v, ttry, tdid, tnext, time, EPS);
  }

  // Now binary chop to find minimum point
  double tlo = 0, thi = tdid;
  while(thi-tlo > BCHOP){
    ttry = (tlo+thi)/2.;
    r    = rold;
    v    = vold;
    time = told;
    gsint(q, r, v, ttry, tdid, tnext, time, EPS);
    if(rdot(r,v) < 0.){
      tlo = ttry;
    }else{
      thi = ttry;
    }
  }

  // Fails to reach potential
  if(rpot(q,r) >= pot){
    x  = r.x();
    y  = r.y();
    vx = v.x();
    vy = v.y();
    return false;
  }

  // Does reach potential, work out where
  double tmin = time;

  strinit(q,r,v);
  time  = 0.;
  tnext = 2.e-4;
  while(rpot(q,r) > pot){
    rold = r;
    vold = v;
    told = time;
    ttry = std::min(tnext,tmin-time);
    gsint(q, r, v, ttry, tdid, tnext, time, EPS);
  }

  // now binary chop again  
  tlo = 0., thi = tdid;
  while(thi-tlo > BCHOP){
    ttry = (tlo+thi)/2.;
    r    = rold;
    v    = vold;
    time = told;
    gsint(q, r, v, ttry, tdid, tnext, time, EPS);
    if(rpot(q,r) > pot){
      tlo = ttry;
    }else{
      thi = ttry;
    }
  }
  x  = r.x();
  y  = r.y();
  vx = v.x();
  vy = v.y();
  return true;
}
