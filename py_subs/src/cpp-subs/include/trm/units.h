#ifndef TRM_UNITS
#define TRM_UNITS

#include <cmath>
#include <string>
#include "trm/subs.h"
#include "trm/constants.h"
#include "trm/fraction.h"

namespace Subs {

    //! A class for the specification of physical units

    /**
     * The Subs::Units class is used to specify physical units. Units can be described in terms
     * of the fundamental quantities mass, length, time, electric charge, solid angle and their
     * SI equivalents metres, seconds, kilograms, Coulombs, steradians. 
     *
     * Any given unit can be specified by a series of exponents for each of these and a 
     * scaling factor to get to it from the SI version. e.g. an Angstrom could be represented
     * by (0,1,0,0,0) and 1.e-10, while a Watt would be (1,2,-3,0,0) and 1. This is how the 
     * Units class operates. It also includes a facility to recognise the simple form of some
     * units. Thus it would report a unit as a Watt rather than kg m**2/s**2. The exponents
     * are internally represented by the ratios of two integers.
     */

    class Units {
    public:

	//! Enumeration of recognised units for constructor
	// If anything is added to this, an equivalent line must be inserted
	// in the check_string routine in Units.cc and at other obvious
	// points in that file
	enum UNITS {
	    ANGSTROM,          /**< length, Angstroms */
	    NANNOMETRE,        /**< length, nannometres */     
	    MICRON,            /**< length, micrometres */
	    CENTIMETRE,        /**< length, centimetres */
	    METRE,             /**< length, metres */
	    KILOMETRE,         /**< length, kilometres */
	    SOLAR_RADIUS,      /**< length, solar radii */
	    AU,                /**< length, Astronomical units */
	    PARSEC,            /**< length, parsecs */
	    KILOPARSEC,        /**< length, kiloparsecs */
	    MEGAPARSEC,        /**< length, megaparsecs */
	    GRAM,              /**< mass, grams */
	    KILOGRAM,          /**< mass, kilograms */
	    SOLAR_MASS,        /**< mass, solar masses */
	    SECOND,            /**< time, seconds */
	    MILLISECOND,       /**< time, milliseconds */
	    MICROSECOND,       /**< time, microseconds */
	    MINUTE,            /**< time, minutes */
	    HOUR,              /**< time, hours */
	    DAY,               /**< time, days */
	    YEAR,              /**< time, years */
	    FLAMBDA_A,         /**< flambda, ergs/s/cm**2/A */
	    FLAMBDA_NM,        /**< flambda, ergs/s/cm**2/nm */
	    FNU_MJY,           /**< fnu, milliJanskys = 10**-26 W/m**2/Hz */
	    AB_MAG,            /**< 16.4 - 2.5*log10(fnu, mJy) */
	    KMS,               /**< speed, km/s */
	    MPS,               /**< speed, m/s */
	    HERTZ,             /**< frequency, cycles/second */
	    CYCPERDAY,         /**< frequency, cycles/day */
	    PERSEC,            /**< frequency, per second, e.g. for count rates */
	    PIXEL,             /**< pixel number, dimensionless */
	    PHASE,             /**< cycles of orbital phase etc, dimensionless */
	    NONE,              /**< dimensionless quantity */
	    USER,              /**< user-defined unit */
	};

	//! Default constructor
	Units() : 
	    mass(0), length(0), time(0), charge(0), solid_angle(0), scale(1.), 
	    pgplot_name(), pgplot_name_set(true) {}

	//! Constructor from enumerated types
	Units(UNITS units, const std::string& name="");

	//! Constructor from a string
	Units(const std::string& sunit);

	//! Gets the PGPLOT name
	std::string get_pgplot_name() const {
	    if(pgplot_name_set)
		return pgplot_name;
	    else
		return "";
	}

	//! Explicitly sets the PGPLOT name
	void set_pgplot_name(const std::string& pgplot_name){
	    this->pgplot_name = pgplot_name;
	    pgplot_name_set = true;
	} 

	//! Sets the PGPLOT name automatically.
	void set_pgplot_name();

	//! Unset the PGPLOT name 
	/** This function leaves the Units name in an undefined state so that
	 * the next attempt to retrieve will cause it to be created automatically.
	 */
	void unset_pgplot_name(){
	    pgplot_name_set = false;
	}

	//! Returns long string completely defining units
	std::string to_string() const;

	//! Get the scale factor
	double get_scale() const { return scale;}

	friend bool operator==(const Subs::Units& unit1, const Subs::Units& unit2);

	friend bool operator!=(const Subs::Units& unit1, const Subs::Units& unit2);

	friend Units operator/(const Subs::Units& unit1, const Subs::Units& unit2);

	friend Units operator/(double factor, const Subs::Units& unit);

	friend Units operator*(const Subs::Units& unit1, const Subs::Units& unit2);

	friend std::ostream& operator<<(std::ostream& ostr, const Subs::Units& unit);

	//! Get scaling factor to turn current units into new units
	double get_scale_factor(const Subs::Units& units) const;

	friend double get_scale_factor(const Subs::Units& unit1, const Subs::Units& unit2);

	//! Rescale 
	void rescale(double factor);

	//! Binary write
	void write(std::ofstream& ostr) const;

	//! Binary read
	void read(std::ifstream& istr, bool swap_bytes);

	//! Binary skip
	static void skip(std::ifstream& istr, bool swap_bytes);
    
	//! Exception class for Units
	class Units_Error : public Subs::Subs_Error {
	public:
	    //! Default constructor
	    Units_Error() : Subs::Subs_Error() {};
	    //! Constructor from a string
	    Units_Error(const std::string& str) : Subs::Subs_Error(str) {}
	};
    
	//! Checks whether a string is one of the recognised units.
	static bool check_string(const std::string& check, UNITS& units);

	//! Returns true if the object has dimensions.
	bool has_dimensions() const;

    private:

	Fraction mass;
	Fraction length;
	Fraction time;
	Fraction charge;
	Fraction solid_angle;
	double scale;
	std::string pgplot_name;
	bool pgplot_name_set;

	struct Info {
	    Info(const std::string& name, const std::string& definition, UNITS units) :
		name(name), definition(definition), units(units) {}
	    std::string name;
	    std::string definition;
	    UNITS units;
	};

    };

    //! Are two units equivalent?
    bool operator==(const Units& unit1, const Units& unit2);

    //! Are two units not equivalent?
    bool operator!=(const Units& unit1, const Units& unit2);

    //! Divide two units
    Units operator/(const Units& unit1, const Units& unit2);

    //! Inverts units while changing scale
    Units operator/(double factor, const Units& unit);

    //! Multiply two units
    Units operator*(const Units& unit1, const Units& unit2);

    //! ASCII output
    std::ostream& operator<<(std::ostream& ostr, const Units& unit);

    //! Get scaling factor to turn units1 numbers into units2 numbers
    double get_scale_factor(const Units& unit1, const Units& unit2);
  
};

#endif










