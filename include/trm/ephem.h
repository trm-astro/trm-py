#ifndef TRM_EPHEMERIS
#define TRM_EPHEMERIS

#include <cmath>
#include <iostream>
#include "trm/subs.h"

namespace Subs {

  //! Ephemeris class

  /** A class for describing ephemerides, computing phases
   * etc. It allows you to store a linear or quadratic ephemeris
   * and then handle standard computations such as deriving a time from
   * a phase, or with more difficulty, a phase from a time. It also stores
   * the type of timescale that the ephemeris refers to. This is for reference
   * only; it is up to the calling program to supply the correct times.
   *
   * The routines handle errors on the assumption that there are no covariances. 
   * It is up to the user to make sure that this is approximately the case. 
   */

  class Ephem {
  public:

    //! Type of ephemeris
    enum ETYPE {
      LINEAR,     /**< Linear ephemeris */
      QUADRATIC   /**< Quadratic ephemeris */
    };

    //! Type of timescale
    enum TSCALE {
      HJD,         /**< Heliocentric UTCs, stored as JD */
      HMJD,        /**< Heliocentric UTCs, stored as MJD */
      BJD,         /**< Barycentric TDB, stored as JD  */
      BMJD         /**< Barycentric TDB, stored as MJD  */
    };

    //! Default constructor
    Ephem() : etype(LINEAR), tscale(HJD), tzero(0.), per(1.), tzerr(0.), perr(0.) {};

    //! Linear ephemeris constructor
    Ephem(double T0, double period, TSCALE tscale);
 
    //! Linear ephemeris with errors
    Ephem(double T0, double period, double T0err, double pererr, TSCALE tscale);
						 
    //! Quadratic ephemeris
    Ephem(double T0, double period, double pdot, TSCALE tscale);

    //! Quadratic ephemeris with errors
    Ephem(double T0, double period, double pdot, double T0err, double pererr, double qderr, TSCALE tscale);
    
    //! Returns with phase for a given time
    double phase(double t) const;

    //! Returns with error on phase for a given time
    double pherr(double t) const;

    //! Returns with time for a given phase
    double time(double p) const;

    //! Returns with time error for a given phase
    double timerr(double p) const;

    //! Returns the type of times
    TSCALE get_tscale() const {return tscale;}
    
    //! Sets a linear ephemeris
    void set(double T0, double period, TSCALE tscale);

    //! Sets a quadratic ephemeris
    void set(double T0, double period, double pdot, TSCALE tscale);

    //! ASCII output of an Ephem
    friend std::ostream& operator<<(std::ostream& ost, const Ephem& ephem);

    //! ASCII input of an Ephem
    friend std::istream& operator>>(std::istream& ist, Ephem& ephem);
    
  private:

    ETYPE  etype;
    TSCALE tscale;
    double tzero, per,  quad;
    double tzerr, perr, qerr;

  };

};

#endif






