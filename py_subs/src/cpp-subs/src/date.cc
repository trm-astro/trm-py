#include "trm/date.h"

#include <iomanip>
#include <sstream>
#include "trm/subs.h"


int Subs::Date::print_method;

/** Constructs a Date from a day, month and year. The values will be range-checked.
 * \param day_   the day of the month (1-31, but checked for each month)
 * \param month_ the month as an enum for clarity.
 * \param year_  the year e.g. 1999, 2002
 */
Subs::Date::Date(int day_, Month month_, int year_){
    set(day_,int(month_),year_);
}

/** Constructs a Date from a day, month and year. The values will be range-checked.
 * \param day_   the day of the month (1-31, but checked for each month)
 * \param month_ the month (1-12)
 * \param year_  the year e.g. 1999, 2002
 */
Subs::Date::Date(int day_, int month_, int year_) {
    set(day_,month_,year_);
}

/** Constructs a Date from a string of the form "01 Mar 2010"
 * \param date   the string containing the date
 */
Subs::Date::Date(const std::string& date){
    set(date);
}

/** Adds a number of days to a Date
 * \param nday the number of days to b added
 */
void Subs::Date::add_day(int nday){
    mjd_ += nday;
}
 
/** Returns the day of the month
 */
int Subs::Date::day() const{
    int d,m,y,status;
    double fd;
	status = iauJd2cal(MJD0, double(mjd_)+0.1,&y, &m, &d, &fd);
    return d;
}
 
/** Returns the month of the year
 */
int Subs::Date::month() const{
    int d,m,y,status;
    double fd;
    status = iauJd2cal(MJD0, double(mjd_)+0.1, &y, &m, &d, &fd);
    return m;
}

/** Returns the year
 */
int Subs::Date::year() const{
    int d,m,y,status;
    double fd;
    status = iauJd2cal(MJD0, double(mjd_)+0.1, &y, &m, &d, &fd);
    return y;
}

/** Returns the day, month and year. Because of the way the date is stored
 * it is more efficient to use this function than separate calls to day()
 * month() and year().
 * \param day_   the day of the month
 * \param month_ the month of the year
 * \param year_  the year
 */
void Subs::Date::date(int& day_, int& month_, int& year_) const {
    int status;
    double fd;
	status = iauJd2cal(MJD0, double(mjd_)+0.1,&year_, &month_, &day_, &fd);
}

/** Returns modified Julian day number,
 * i.e. JD - 2400000.5 at the start of the day 
 * in question.
 */
double Subs::Date::mjd() const {
return mjd_;
}

/** Set the date
 * \param day_   the day of the month (1-31)
 * \param month_ the month of the year (1-12)
 * \param year_ the day of the month (e.g. 2003)
 */

void Subs::Date::set(int day_, int month_, int year_){
    int status;
    double mj;
	double mjd;
    valid_date(day_,month_,year_);
	// mjd here is 2400000.5 and mj is the fractional part of the day
	status = iauCal2jd(year_, month_, day_, &mjd, &mj);
	
	//mjd_ = int(mj+0.1);
	// Set mjd_ internally to be the modified Julian day number (float)
	mjd_ = mj;
}

/** Set the date from a string
 * \param date   a date string. e.g. 1 Jan 2002 or 17 Nov 1961 or 17/11/1961
 */
void Subs::Date::set(const std::string& date){
    std::istringstream ist(date);
    int day_, month_, year_;
    std::string mname;
    ist >> day_;
    if(!ist) throw Date_Error("Subs::Date::set(const std::string&): failed to read date = " + date);
    char c;
    ist.get(c);
    if(!ist) throw Date_Error("Subs::Date::set(const std::string&): failed to read date = " + date);
    if(c == '/'){
	ist >> month_;
	if(!ist) throw Date_Error("Subs::Date::set(const std::string&): failed to read date = " + date);
	ist.ignore();
	ist >> year_;
	if(!ist) throw Date_Error("Subs::Date::set(const std::string&): failed to read date = " + date);
    }else{
	ist >> mname >> year_;  

	mname = Subs::toupper(mname);
	if(mname == "JAN"){
	    month_ = Jan;
	}else if(mname == "FEB"){
	    month_ = Feb;
	}else if(mname == "MAR"){
	    month_ = Mar;
	}else if(mname == "APR"){
	    month_ = Apr;
	}else if(mname == "MAY"){
	    month_ = May;
	}else if(mname == "JUN"){
	    month_ = Jun;
	}else if(mname == "JUL"){
	    month_ = Jul;
	}else if(mname == "AUG"){
	    month_ = Aug;
	}else if(mname == "SEP"){
	    month_ = Sep;
	}else if(mname == "OCT"){
	    month_ = Oct;
	}else if(mname == "NOV"){
	    month_ = Nov;
	}else if(mname == "DEC"){
	    month_ = Dec;
	}else{
	    throw Date_Error(std::string("Subs::Date::set(const std::string&): unrecognised month = ") + mname);
	}
    }
    set(day_,month_,year_);
}

/** Set the date from an MJD (JD-2400000.5)
 * \param mjd MJD at start of date in question
 */
void Subs::Date::set(int mjd){
    mjd_ = mjd;
}

/** Format a date as a string
 */
std::string Subs::Date::str() const {

    int d, m, y;
    date(d,m,y);

    std::ostringstream ost;

    ost << std::setw(2) << std::setfill('0') << d;
  
    if(print_method == 2){
	ost << "/" << std::setw(2) << std::setfill('0') << m << "/";
    }else{
	switch(m){
	    case Jan:
		ost << " Jan ";
		break;
	    case Feb:
		ost << " Feb ";
		break;
	    case Mar:
		ost << " Mar ";
		break;
	    case Apr:
		ost << " Apr ";
		break;
	    case May:
		ost << " May ";
		break;
	    case Jun:
		ost << " Jun ";
		break;
	    case Jul:
		ost << " Jul ";
		break;
	    case Aug:
		ost << " Aug ";
		break;
	    case Sep:
		ost << " Sep ";
		break;
	    case Oct:
		ost << " Oct ";
		break;
	    case Nov:
		ost << " Nov ";
		break;
	    case Dec:
		ost << " Dec ";
		break;
	    default:
		throw Date_Error("Subs::Date::str(): unrecognised month");
	}
    }
    ost << std::setw(4) << y << std::setfill(' ');
    return ost.str();
}

/** Computes the day of the week corresponding to the Date
 * \return A string representing the day of the week, as in Monday
 */
std::string Subs::Date::day_of_week() const {

    std::string str;

    switch((int(mjd()+0.1) + 2) % 7){
	case 0:
	    str = "Monday";
	    break;
	case 1:
	    str =  "Tuesday";
	    break;
	case 2:
	    str =  "Wednesday";
	    break;
	case 3:
	    str =  "Thursday";
	    break;
	case 4:
	    str =  "Friday";
	    break;
	case 5:
	    str =  "Saturday";
	    break;
	case 6:
	    str =  "Sunday";
	    break;
	default:
	    str =  "day";
    }
    return str;
}

/** Computes an integer representing the day of the week, starting at 0 for
 * Sunday, 1 for Monday, etc. The switch on the Saturday/Sunday boundary matches
 * times from the GPS.
 * \return Integer from 0 to 6, 0 = Sunday.
 */
int Subs::Date::int_day_of_week() const {
    return (int(mjd()+0.1) + 3) % 7;
}

/** Reads a date stored in binary format
 */
void Subs::Date::read(std::ifstream& s, bool swap_bytes){
    INT4 mj;
    s.read((char*)&mj,sizeof(INT4));
    if(!s) throw Date_Error("void Subs::Date::read(std::ifstream&): error reading date");
    if(swap_bytes) mj = Subs::byte_swap(mj);
    if(mj < 0 || mj > 1000000)
	throw Date_Error("void Subs::Date::read(std::ifstream&): error date out of range");

    mjd_ = mj;
}

void Subs::Date::write(std::ofstream& s) const {
    s.write((char*)&mjd_,sizeof(INT4));
    if(!s) throw Date_Error("void Subs::Date::write(std::ofstream&) const: error writing date");
}

void Subs::Date::skip(std::ifstream& s){
    s.ignore(sizeof(INT4));
    if(!s) throw Date_Error("void Subs::Date::skip(std::ifstream&): error skipping date");
}

/** Returns 0 or 1 to add to the number of days in February according to whether
 * or not 'year' is a leap year. Leap years are if 'year' is divisble by 400
 * or divible by 4, but not the latter when they are divisible by 100 (with the 400
 * rule taking precedence).
 * \param year the year in question
 * \return 0 or 1, the number of days to add to February
 */

int Subs::Date::leapyear(int year){
    if(year % 400 == 0) return 1;
    if(year % 100 == 0) return 0;
    if(year % 4 == 0) return 1;
    return 0;
}

/** Defines valid dates.
 * \param day the day of month
 * \param month the month of the year
 * \param year the year
 * \exception Throws a Subs::Date::Date_Error if the date given is invalid.
 */

void Subs::Date::valid_date(int day, int month, int year){

    if(year < -4699)
	throw Date_Error("void Subs::Date::valid_date(int, int, int): invalid date. Year = " + Subs::str(year) + " is less than -4699");

    if(month < 1 || month > 12)
	throw Date_Error("void Subs::Date::valid_date(int, int, int): invalid date. Month = " + Subs::str(month) + " is out of range 1 to 12");

    int mx = 0;
    switch(month){
	case Subs::Date::Feb:
	    mx = 28 + leapyear(year);
	    break;
	case Subs::Date::Sep: case Subs::Date::Apr: case Subs::Date::Jun: case Subs::Date::Nov:
	    mx = 30;
	    break;
	case Subs::Date::Jan: case Subs::Date::Mar: case Subs::Date::May: 
	case Subs::Date::Jul: case Subs::Date::Aug: case Subs::Date::Oct: case Subs::Date::Dec:
	    mx = 31;
    }
    if(day < 1 || day > mx)
	throw Date_Error("void Subs::Date::valid_date(int, int, int): invalid date. Day = " + Subs::str(day) + " is out of range 1 to " + Subs::str(mx) );
}

/** ASCII input. Reads dates of the form 17/11/1961 and (2) '17 Nov 1961' 
 * (3 letter month names only).
 */

std::istream& Subs::operator>>(std::istream& ist, Subs::Date& date){

    if(!ist) return ist;

    int day_, month_, year_;
    std::string mname;
    ist >> day_;
    if(!ist) return ist;
    char c;
    ist.get(c);
    if(!ist) return ist;
    if(c == '/'){
	ist >> month_;
	if(!ist) return ist;
	ist.ignore();
	ist >> year_;
	if(!ist) return ist;
    }else{
	ist >> mname >> year_;  

	mname = Subs::toupper(mname);
	if(mname == "JAN"){
	    month_ = Date::Jan;
	}else if(mname == "FEB"){
	    month_ = Date::Feb;
	}else if(mname == "MAR"){
	    month_ = Date::Mar;
	}else if(mname == "APR"){
	    month_ = Date::Apr;
	}else if(mname == "MAY"){
	    month_ = Date::May;
	}else if(mname == "JUN"){
	    month_ = Date::Jun;
	}else if(mname == "JUL"){
	    month_ = Date::Jul;
	}else if(mname == "AUG"){
	    month_ = Date::Aug;
	}else if(mname == "SEP"){
	    month_ = Date::Sep;
	}else if(mname == "OCT"){
	    month_ = Date::Oct;
	}else if(mname == "NOV"){
	    month_ = Date::Nov;
	}else if(mname == "DEC"){
	    month_ = Date::Dec;
	}else{
	    ist.setstate(std::ios::failbit);
	    return ist;
	}
    }
    try{
	date.set(day_,month_,year_);
    }
    catch(...){
	ist.setstate(std::ios::failbit);
    }
    return ist;
}

std::ostream& Subs::operator<<(std::ostream& ost, const Subs::Date& date){

    if(!ost) return ost;

    int day, month, year;
    date.date(day,month,year);

    // First output the day of the month

    ost.setf(std::ios::right);
    ost << std::setw(2) << std::setfill('0') << day;

    // Now the month of the year

    if(Subs::Date::print_method == 2){
	ost << "/" << std::setw(2) << std::setfill('0') << month << "/";
    }else{
	switch(month){
	    case Subs::Date::Jan:
		ost << " Jan ";
		break;
	    case Subs::Date::Feb:
		ost << " Feb ";
		break;
	    case Subs::Date::Mar:
		ost << " Mar ";
		break;
	    case Subs::Date::Apr:
		ost << " Apr ";
		break;
	    case Subs::Date::May:
		ost << " May ";
		break;
	    case Subs::Date::Jun:
		ost << " Jun ";
		break;
	    case Subs::Date::Jul:
		ost << " Jul ";
		break;
	    case Subs::Date::Aug:
		ost << " Aug ";
		break;
	    case Subs::Date::Sep:
		ost << " Sep ";
		break;
	    case Subs::Date::Oct:
		ost << " Oct ";
		break;
	    case Subs::Date::Nov:
		ost << " Nov ";
		break;
	    case Subs::Date::Dec:
		ost << " Dec ";
		break;
	    default:
		ost << " UNKNOWN MONTH ";
	}
    }
    ost << std::setw(4) << year << std::setfill(' ');
    return ost;
}

//! Defines 'greater than' for two Dates
bool Subs::operator>(const Subs::Date& d1, const Subs::Date& d2){
    return d1.mjd_ > d2.mjd_;
}

//! Defines 'less than' for two Dates
bool Subs::operator<(const Subs::Date& d1, const Subs::Date& d2){
    return d1.mjd_ < d2.mjd_;
}

//! Defines 'greater than or equal to' for two Dates
bool Subs::operator>=(const Subs::Date& d1, const Subs::Date& d2){
    return d1.mjd_ >= d2.mjd_;
}

//! Defines 'less than or equal to' for two Dates
bool Subs::operator<=(const Subs::Date& d1, const Subs::Date& d2){
    return d1.mjd_ <= d2.mjd_;
}

//! Defines 'equals' for two Dates
bool Subs::operator==(const Subs::Date& d1, const Subs::Date& d2){
    return d1.mjd_ == d2.mjd_;
}

//! Subtracts two dates, defining the operation '-'
int Subs::operator-(const Subs::Date& d1, const Subs::Date& d2){
    return d1.mjd_ - d2.mjd_;
}












