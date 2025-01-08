#ifndef TRM_POSITION
#define TRM_POSITION

#include <cmath>
#include <string>
#include <sofa.h>
#include "trm/subs.h"
#include "trm/constants.h"
#include "trm/date.h"
#include "trm/time.h"
#include "trm/telescope.h"
#include "trm/vec3.h"

namespace Subs {

    //! Structure for return of accurate target data
    struct Pinfo {
	double tcor_bar;    /**< Time in seconds that when added to the local gives equivalent at barycentre */
	double tcor_hel;    /**< Time in seconds that when added to the local gives equivalent at heliocentre */
	double vearth_bar;  /**< Apparent velocity of target owing to Earth's motion relative to barycentre, km/s */
	double vearth_hel;  /**< Apparent velocity of target owing to Earth's motion relative to heliocentre, km/s */
    };

    //! A class for positions of astronomical objects on the sky

    /**
     * Astronomical position class. This stores the positions of objects in ICRS
     * coordinates only but allows for any epoch to define the position given 
     * proper motions and radial velocity.
     */

    class Position {
    public:

	//! Default constructor
	Position() : r_a(0.), decl(0.), ep(2000.), pm_ra(0.), pm_dec(0.), prlx(0.), r_v(0.) {};

	//! Constructor for general RA, Dec & Epoch (but no proper motion)
	Position(int rah, int ram, float ras, char decsgn, int decd, int decm, float decs, double epoch);

	//! Constructor for general RA, Dec & Epoch (but no proper motion)
	Position(double ra, double dec, double epoch);

	//! Constructor from a string
	Position(const std::string& pos);

	//! Sets value from a string
	void set(const std::string& pos);
    
	//! Returns the RA in hours
	double ra() const {return r_a;}  

	//! Returns the declination in degrees
	double dec() const {return decl;}  

	//! Returns the epoch of the coordinates  
	double epoch() const {return ep;}

	//! Returns the proper motion in arcseconds of RA/yr
	double pmra() const {return pm_ra;}

	//! Returns the proper motion in arcseconds of dec/yr
	double pmdec() const {return pm_dec;}    

	//! Returns the parallax in arcseconds
	double parallax() const {return prlx;}    

	//! Returns the radial velocity in km/s (+ve receding)
	double rv() const {return r_v;}    

	//! Returns the RA in radians
	double rar() const {return Constants::TWOPI*r_a/24.;}  

	//! Returns the declination in radians
	double decr() const {return Constants::TWOPI*decl/360.;}

	//! Returns the proper motion in radians of RA/yr (i.e. large near poles)
	double pmrar() const {return Constants::TWOPI*pm_ra/std::cos(decr())/360./3600.;}

	//! Returns the proper motion in radians of dec/yr
	double pmdecr() const {return Constants::TWOPI*pm_dec/360./3600;}    

	//! Sets the RA in hours
	double& ra() {return r_a;}  

	//! Sets the proper motion in arcseconds of RA/yr
	float& pmra() {return pm_ra;}  

	//! Sets the declination in degrees
	double& dec() {return decl;}  

	//! Sets the proper motion in arcseconds of dec/yr
	float& pmdec() {return pm_dec;}  

	//! Sets the epoch
	double& epoch() {return ep;}  

	//! Sets the parallax in arcseconds
	float& parallax() {return prlx;}    

	//! Sets the radial velocity in km/s (+ve receding)
	float& rv() {return r_v;}    

	//! Set RA, Dec from their value in radians
	void set_radec(double rar, double decr){
	    r_a  = 24.*rar/Constants::TWOPI;
	    decl = 360.*decr/Constants::TWOPI;
	}

	//! Returns a string of the position
	std::string ra_dec() const;
    
	//! Update coordinates to a new epoch, correcting for space motion
	void   update(double epch);

	//! Update coordinates to match the epoch of another Position, correcting for space motion
	void   update(const Position& pos);

	//! Returns vector towards Position.
	Vec3 vect() const;

	//! Heliocentric time correction
	double tcorr_hel(const Time& time, const Telescope& tel) const;

	//! Barycentric time correction
	double tcorr_bar(const Time& time, const Telescope& tel) const;

	//! Computes velocity and time correction info
	Pinfo pinfo(const Time& time, const Telescope& tel) const;

	//! Sets the position to that of the Sun
	void set_to_sun(const Time& time, const Telescope& tel);

	//! Computes observational information
	Altaz  altaz(const Time& time, const Telescope& tel) const;

	//! Write position to a binary file
	void write(std::ofstream& s) const;

	//! Read position from a binary file
	void read(std::ifstream& s, bool swap_bytes);

	//! Skip Position in a binary file
	static void skip(std::ifstream& s);    

	friend std::istream& operator>>(std::istream& ist, Position& pos);
    
	//! Exception class for Positions
	class Position_Error : public Subs::Subs_Error {
	public:
	    //! Default constructor
	    Position_Error() : Subs::Subs_Error() {};
	    //! Constructor from a string
	    Position_Error(const std::string& str) : Subs::Subs_Error(str) {}
	};
    
    private:

	// r_a  = right ascension, hours of RA, ICRS
	// decl = declination, degrees, ICRS
	// ep   = epoch of coordinates for applying proper motion.
	REAL8 r_a, decl, ep;

	// Proper motions in arcsec per year of RA and Dec.
	REAL4 pm_ra, pm_dec;

	// Parallax (arcsec) and radial velocity (km/s)
	REAL4 prlx, r_v;

    };
  
    //! ASCII output
    std::ostream& operator<<(std::ostream& ost, const Subs::Position& pos);

    //! ASCII input
    std::istream& operator>>(std::istream& ist, Position& pos);

    //! Take dot product of two positions
    double dot(const Subs::Position& pos1, const Subs::Position& pos2);
  
};

#endif










