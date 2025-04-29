#ifndef TRM_BINARY
#define TRM_BINARY

#include <math.h>
#include "trm/star.h"
#include "trm/ephem.h"

namespace Subs {

  //! Class for representing binary stars

  /** Class for stars with names, positions and ephemerides.
   * It defines an object which contains information
   * on a binary star including its name, position and ephemeris.
   * It does this through inheritance of the Star and Ephem classes.
   */


  class Binary : public Star, public Ephem {

  public:
    
    //! Default constructor
    Binary() : Star(), Ephem() {};

    //! Constructor from a star and an ephemeris
    Binary(const Star& star, const Ephem& eph) : Star(star), Ephem(eph) {}

    //! Sets the values of the position and ephemeris
    void set(const Star& star, const Ephem& eph);
    
    //! Does this have an ephemeris?
    bool has_ephem() const {return true;}

    //! What is the timescale of the ephemeris
    Ephem::TSCALE get_tscale() const {return Ephem::get_tscale();}

    //! Returns with phase for a given time
    double phase(double t) const {return Ephem::phase(t);}

    //! Returns with error on phase for a given time
    double pherr(double t) const {return Ephem::pherr(t);}

    //! Returns with time for a given phase
    double time(double p) const {return Ephem::time(p);}

    //! Returns with error in time for a given phase
    double timerr(double p) const {return Ephem::timerr(p);}

    friend std::ostream& operator<<(std::ostream& ost, const Binary& binary);
    
  };

};

#endif










