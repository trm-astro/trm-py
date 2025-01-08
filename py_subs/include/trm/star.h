#ifndef TRM_STAR
#define TRM_STAR

#include <math.h>
#include "trm/subs.h"
#include "trm/position.h"
#include "trm/ephem.h"

namespace Subs {

  //! Represents a star
  /** Class to represent a star holding its name and position. It also has
   * methods that can be inherited if the star is a binary.
   */
  class Star : public Position {

  public:

    //! Default constructor
    Star() : Position(), nam() {};

    //! General constructor
    Star(const std::string& star, const Position& position)
      : Position(position), nam(star) {};
    
    //! Virtual destructor
    virtual ~Star() {};
    
    std::string name() const { return nam;};
    
    friend std::istream& operator>>(std::istream& ist, Star& pos);
    
    //! Does this have an ephemeris
    virtual bool has_ephem() const {return false;}

    //! What is the timescale of the ephemeris (if there is one)
    virtual Ephem::TSCALE get_tscale() const {throw Subs_Error("get_tscale(): star = " + nam + " has no ephemeris");}

    //! Returns an orbital phase given a time
    virtual double phase(double t) const {throw Subs_Error("phase(double): star = " + nam + " has no ephemeris");}

    //! Returns an error in the orbital phase given a time
    virtual double pherr(double t) const {throw Subs_Error("pherr(double): star = " + nam + " has no ephemeris");}

    //! Returns a time give a phase
    virtual double time(double p) const {throw Subs_Error("time(double): star = " + nam + " has no ephemeris");}

    //! Returns an error in the time give a phase
    virtual double timerr(double p) const {throw Subs_Error("timerr(double): star = " + nam + " has no ephemeris");}

  private:

    std::string nam;

  };

  std::ostream& operator<<(std::ostream& ost, const Star& star);
};

#endif










