#ifndef TRM_TIME
#define TRM_TIME

#include <string> 
#include <ctime> 
#include <sofa.h>
#include "trm/date.h"
#include "trm/subs.h"

namespace Subs {

    class Telescope; 
    class Vec3;

    /** Class to represent times with high precision. The times stored are taken to represent
     * UTC,  that is coordinated UT, which tracks TAI (atomic time) to an integral number of seconds).
     * Various functions allow this to be converted to other forms of time.
     */

    class Time : public Date {

    public:
    
	//! represents time as hours, minutes, seconds, fractions of secs
    
	struct HMS{
	    //! Number of hours of the day (0 to 23)
	    int hour;
	    //! Number of minutes of the hour (0 to 59)
	    int min;
	    //! Number of seconds of the minute (0 to 59)
	    int sec;
	    //! Fraction of a second (0. to 0.9999999999)
	    float fsec;
	};
    
	//! Default constructor, sets to 0 hours
	Time() : Date(), hour_(0.) {};

	//! constructor from a date and decimal hour
	Time(int day_, int month_, int year_, double hour=0.);

	//! constructor from a date and decimal hour
	Time(int day_, Date::Month month_, int year_, double hour=0.);

	//! constructor from a date and hours, minutes, secs
	Time(int day_, int month_, int year_, int hour, int min, double sec);

	//! constructor from a date and hours, minutes, secs
	Time(int day_, Date::Month month_, int year_, int hour, int min, double sec);

	//! constructor from a Date and decimal hour
	Time(const Date& date, double hour=0.f);

	//! constructor from a string
	Time(const std::string& time);

	//! constructor from a decimal MJD
	Time(double mjd);
    
	//! returns the decimal hour
	double hour() const {return hour_;}

	//! returns the hours, minutes, seconds structure
	HMS hms() const;

	//! returns the decimal modified Julian day
	double mjd() const;  

	void hack_report() const;

	//! returns the Julian epoch
	double jepoch() const;  

	//! returns the difference between terrestrial time and UTC (TT-UTC) in seconds
	double dtt() const;  

	//! returns terrestrial time (MJD TT)
	double tt() const;  

	//! returns Barycentric Coordinate Time minus TT (TDB-TT)
	double dtdb(const Telescope& tel) const;

	//! returns Barycentric Coordinate Time
	double tdb(const Telescope& tel) const;

	//! computes earth position relative to barycentre
	Vec3 earth_pos_bar(const Telescope& tel) const;

	//! computes earth position relative to heliocentre
	Vec3 earth_pos_hel(const Telescope& tel) const;

	//! computes earth position and velocity relative to heliocentre and barycentre
	void earth(const Telescope& tel, Vec3& ph, Vec3& vh, Vec3& pb, Vec3& vb) const;

	//! computes useful info (altitude etc) for the Sun
	Altaz sun(const Telescope& tel) const;

	//! returns a string for output
	std::string str() const;
    
	//! Sets the decimal hour
	void   set_hour(double hour);

	//! Sets the time with day, month etc
	void   set(int day_, Date::Month month_, int year_, int hour, int min, double sec);

	//! Sets the time with day, month etc
	void   set(int day_, int month_, int year_, int hour, int min, double sec);

	//! Sets the time from a string
	void   set(const std::string& time);

	//! Sets the time from a decimal modified Julian day
	void   set(double mjd);

	//! Sets to the current time
	void   set();

	//! Sets the time corresponding to the start of 'date'
	void   set(const Date& date, double hour=0.);
    
	//! Adds decimal hours to the time
	void   add_hour(double hour);

	//! Adds decimal seconds to the time
	void   add_second(double second);

	//! Returns GMST
	double GMST() const;

	//! Converts GMT to GST
	void   GMTtoGST();
    
	friend bool operator>(const Time& t1, const Time& t2);

	friend bool operator<(const Time& t1, const Time& t2);

	friend bool operator>=(const Time& t1, const Time& t2);

	friend bool operator<=(const Time& t1, const Time& t2);

	friend bool operator==(const Time& t1, const Time& t2);

	friend double operator-(const Time& t1, const Time& t2);
    
	//! Binary input
	void read(std::ifstream& s, bool swap_bytes);

	//! Binary output
	void write(std::ofstream& s) const;

	//! Binary skip
	static void skip(std::ifstream& s);
    
	//! Exception class for Time objects
	class Time_Error : public Subs_Error {
	public:
	    //! Default constructor
	    Time_Error() : Subs_Error() {};
	    //! Constructor from a string
	    Time_Error(const std::string& str) : Subs_Error(str) {}
	};
    
	//! Tests a time for validity
	static void valid_time(int hour, int minute, double second);
  
    private:
    
	double hour_;
    
    };
  
    //! ASCII output
    std::ostream& operator<<(std::ostream& ost, const Time& time);

    //! ASCII input
    std::istream& operator>>(std::istream& ist, Time& time);

    //! Is t1 > t2 ?
    bool operator>(const Time& t1, const Time& t2);
  
    //! Is t1 < t2 ?
    bool operator<(const Time& t1, const Time& t2);
  
    //! Is t1 >= t2 ?
    bool operator>=(const Time& t1, const Time& t2);
  
    //! Is t1 <= t2 ?
    bool operator<=(const Time& t1, const Time& t2);
  
    //! Is t1 == t2 ?
    bool operator==(const Time& t1, const Time& t2);
  
    //! Returns  t1 - t2 in seconds
    double operator-(const Time& t1, const Time& t2);

};

#endif



















