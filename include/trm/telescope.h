#ifndef TRM_TELESCOPE
#define TRM_TELESCOPE

#include <string>
#include <math.h>
#include "trm/constants.h"
#include "trm/subs.h"

namespace Subs {

  //! Class for telescope positions etc.

  /** Telescope bundles together the position, height
   * and name of a telescope. 
   * 
   * A positive sign is used for northern latitudes and for *easterly*
   * longitudes (same as RA).
   */

  class Telescope{

  public:
    
    //! Default constructor
    Telescope() : telescope_name(""), site_name(""), lng(0.), 
      lat(0.), hgt(0.) {};
    
    //! General constructor
    Telescope(const std::string& tname, const std::string& sname, int longd, int longm, float longs, 
	      char eorw, int latd, int latm, float lats, char sorn, float height);

    //! Constructor from a string
    Telescope(const std::string& name);

    //! Set from a string
    void set(const std::string& name);

    //! Retrieve name of telescope
    const std::string& telescope() const {return telescope_name;}  

    //! Retrieve site name
    const std::string& site() const {return site_name;}  

    //! Longitude of telescope
    double longitude() const {return lng;}

    //! Latitude of telescope
    double latitude() const {return lat;}  

    //! Longitude of telescope in radians
    double longituder() const {return Constants::TWOPI*lng/360.;}

    //! Latitude of telescope in radians
    double latituder() const {return Constants::TWOPI*lat/360;}  

    //! Height of telescope
    float  height() const {return hgt;}

    //! Checks whether two telescopes have the same position and height.
    friend bool operator!=(const Telescope& tel1, const Telescope& tel2);

    //! Write telescope to a binary file
    void write(std::ofstream& s) const;

    //! Read telescope from a binary file
    void read(std::ifstream& s, bool swap_bytes);

    //! Skip Telescope in a binary file
    static void skip(std::ifstream& s, bool swap_bytes);    

    friend std::istream& operator>>(std::istream& ist, Telescope& pos);

    //! Error class for Telescope objects
    class Telescope_Error : public Subs_Error {
    public:
      Telescope_Error() : Subs_Error() {};
      Telescope_Error(const std::string& err) : Subs_Error(err) {};
    };
  
  private:
    std::string telescope_name, site_name;
    REAL8 lng, lat; // degrees
    REAL4  hgt; // metres
  };

  //! ASCII output of a Telescope object
  std::ostream& operator<<(std::ostream& ost, const Telescope& obs);

  //! ASCII input of a Telescope object
  std::istream& operator>>(std::istream& ost, Telescope& obs);

  //! Tests for inequality between two Telescopes
  bool operator!=(const Telescope& tel1, const Telescope& tel2);

};

#endif






