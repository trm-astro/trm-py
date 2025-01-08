#ifndef TRM_DATE
#define TRM_DATE

#include <iostream>
#include <string>
#include <sofa.h>
#include "trm/subs.h"

// Define MJD0 as the 0 for the MJD system, i.e. 1858 Nov 17.5
#define MJD0 2400000.5

namespace Subs {

  //! Class for representing dates.

  /** Simple class for storage and manipulation of dates, based upon
   * the Gregorian calendar. This allows one to add days, compare dates
   * etc.
   */

  class Date {

  public:

    //! Controls the ASCII output of a date
    static int print_method;
    
    //! Defines the months
    enum Month {Jan=1,Feb,Mar,Apr,May,Jun,Jul,Aug,Sep,Oct,Nov,Dec};
    
    //! Constructs default date 
    Date() : mjd_(50000) {};

    //! Constructs a date
    Date(int day_, int month_, int year_);

    //! Constructs a date
    Date(int day_, Month month_, int year_);

    //! Constructs a date from a string
    Date(const std::string& date);
    
    //! Adds a number of days to a date
    void   add_day(int nday);
    
    //! Returns the day of month
    int day() const;

    //! Returns the month
    int month() const;

    //! Returns the year
    int year() const;

    //! Returns the day, month & year
    void date(int& day_, int& month_, int& year_) const;

    //! Returns the modified julian day number
    double mjd() const;

    //! Returns a date string
    std::string str() const;

    //! Returns the day of week as a string
    std::string day_of_week() const;
    
    //! Returns the day of the week as an integer
    int int_day_of_week() const;

    //! Sets a date
    void  set(int day_, int month_, int year_);

    //! Sets a date from a string
    void   set(const std::string& date);

    //! Sets a date from a modified Julian day number
    void   set(int mjd);
    
    friend bool operator>(const Date& d1, const Date& d2);

    friend bool operator<(const Date& d1, const Date& d2);

    friend bool operator>=(const Date& d1, const Date& d2);

    friend bool operator<=(const Date& d1, const Date& d2);

    friend bool operator==(const Date& d1, const Date& d2);

    friend int operator-(const Date& d1, const Date& d2);
    
    //! Binary input of a Date
    void read(std::ifstream& s, bool swap_bytes);

    //! Binary output of a Date
    void write(std::ofstream& s) const;

    //! Binary skip of a Date
    static void skip(std::ifstream& s);
    
    //! Error class for Dates
    class Date_Error : public Subs::Subs_Error {
    public:
      //! Default constructor
      Date_Error() : Subs::Subs_Error() {};
      //! Constructor from a string
      Date_Error(const std::string& str) : Subs::Subs_Error(str) {}
    };

    //! Number of days to add to account for leap years
    static int leapyear(int year);
  
    //! Tests the validity of the supplied date
    static void valid_date(int day, int month, int year);
    
  private:
    
    INT4 mjd_;

  };

  //! ASCII output of a Date
  std::ostream& operator<<(std::ostream& ost, const Date& date);
  
  //! ASCII output of a month
  std::ostream& operator<<(std::ostream& ost, const Date::Month& m);
  
  //! ASCII input of a Date
  std::istream& operator>>(std::istream& ist, Date& date);

  //! Is d1 > d2 ?
  bool operator>(const Date& d1, const Date& d2);

  //! Is d1 < d2 ?
  bool operator<(const Date& d1, const Date& d2);

  //! Is d1 >= d2?
  bool operator>=(const Date& d1, const Date& d2);
  
  //! Is d1 <= d2?
  bool operator<=(const Date& d1, const Date& d2);
  
  //! Is d1 == d2?
  bool operator==(const Date& d1, const Date& d2);
  
  //! Returns d1 - d2 in days
  int operator-(const Date& d1, const Date& d2);

};

#endif



















