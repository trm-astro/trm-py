#include "trm/units.h"
#include <cstdio>
#include <iomanip>

/** Static member function to check a string for matching one of the 
 * recognised unit types. The check is case insensitive and unique
 * shortenings are allowed. This is meant for interactive use and will
 * list all recognised units if a unique match is not found.
 * \param check the string to check
 * \param units the enumeration representing the units when matched
 * \return true if a unique match is found.
 */
bool Subs::Units::check_string(const std::string& check, UNITS& units){
  static bool first = true;
  static std::vector<Info> inform;
  if(first){
    first = false;
    inform.push_back(Info("ANGSTROM",   "length = 1e-10 metres", ANGSTROM));
    inform.push_back(Info("NM",         "length = 1e-9 metres", NANNOMETRE));
    inform.push_back(Info("MICRON",     "length = 1e-6 metres", MICRON));
    inform.push_back(Info("CM",         "length = 1e-2 metres", CENTIMETRE));
    inform.push_back(Info("METRE",      "length = 1 metre", METRE));
    inform.push_back(Info("KM",         "length = 1e3 metres", KILOMETRE));
    inform.push_back(Info("RSUN",       "length = 1 solar radius", SOLAR_RADIUS));
    inform.push_back(Info("AU",         "length = 1 astronomical unit", AU));
    inform.push_back(Info("PC",         "length = 1 parsec", PARSEC));
    inform.push_back(Info("KPC",        "length = 1000 parsecs", KILOPARSEC));
    inform.push_back(Info("MPC",        "length = 1e6 parsecs", MEGAPARSEC));
    inform.push_back(Info("GM",         "mass = 1.e-3 kilograms", GRAM));
    inform.push_back(Info("KG",         "mass = 1 kilogram", KILOGRAM));
    inform.push_back(Info("MSUN",       "mass = 1 solar mass", SOLAR_MASS));
    inform.push_back(Info("SECOND",     "time = 1 second", SECOND));
    inform.push_back(Info("MILLISECOND","time = 0.001 second", MILLISECOND));
    inform.push_back(Info("MICROSECOND","time = 1e.-6 second", MICROSECOND));
    inform.push_back(Info("MINUTE",     "time = 60 seconds", MINUTE));
    inform.push_back(Info("HOUR",       "time = 3600 seconds", HOUR));
    inform.push_back(Info("DAY",        "time = 86400 seconds", DAY));
    inform.push_back(Info("YEAR",       "time = 1 year", YEAR));
    inform.push_back(Info("FLAMBDA_A",  "flux density, ergs/s/cm**2/A", FLAMBDA_A));
    inform.push_back(Info("FLAMBDA_NM", "flux density, ergs/s/cm**2/nm", FLAMBDA_NM));
    inform.push_back(Info("FNU_MJY",    "flux density, milliJanskys", FNU_MJY));
    inform.push_back(Info("AB_MAG",     "AB magnitude", AB_MAG));
    inform.push_back(Info("MPS",        "speed = 1 m/s",  MPS));
    inform.push_back(Info("KMS",        "speed = 1 km/s", KMS));
    inform.push_back(Info("HERTZ",      "frequency, cycles/second", HERTZ));
    inform.push_back(Info("CYCPERDAY",  "frequency, cycles/day", CYCPERDAY));
    inform.push_back(Info("PIXEL",      "pixel number, pixel", PIXEL));
    inform.push_back(Info("PHASE",      "phase, cycles", PHASE));
    inform.push_back(Info("NONE",       "dimensionless quantity", NONE));
    inform.push_back(Info("USER",       "a user-defined unit", USER));
  }

  // Now search for matching string, which must be unique.
  std::string trial = Subs::toupper(check);
  int nmatch = 0;

  std::vector<Info>::const_iterator cit, cit_save;
  for(cit = inform.begin(); cit != inform.end(); cit++){
    if(cit->name == trial){
      cit_save = cit;
      nmatch++;
      if(nmatch > 1) break;
    }
  }

  if(nmatch == 1){
    units = cit_save->units;
    return true;
  }
  if(nmatch > 1)
    std::cerr << "Multiple matches were found to the std::string = " << check << std::endl;
  else
    std::cerr << "No matches were found to the std::string = " << check << std::endl;

  std::cerr << "\nList of options follows:\n\n";
  const int NSPRINTF=1024;
  char sprint_out[NSPRINTF];
  for(size_t i=0; i<inform.size(); i++){
      sprintf(sprint_out, "%2i) %-20s %s", int(i+1), inform[i].name.c_str(), inform[i].definition.c_str());
      std::cerr << sprint_out << std::endl;
  }
  return false;
}

/* Constructor from a string which can either be a recognised type as in "ANGSTROM" or
 * built from scratch as in "user 1 0 -1 0 0 \"kg s\\d-1\\u\"". Here the string 'user'
 * identifies a roll-your-own definition, while the exponents refer to mass, length, time,
 * electric charge and solid angle.
 * \param sunit a string to define the units
 * \exception This will throw a Units_Error exception if it cannot interpret the string.
 */
Subs::Units::Units(const std::string& sunit){
  std::istringstream istr(sunit);
  std::vector<std::string> vtrial = Subs::read_line(istr);
  if(vtrial.size() == 0) 
    throw Units_Error("Subs::Units::Units(const std::string&): exception: std::string = " + sunit + " is invalid.");

  UNITS units;
  if(Subs::Units::check_string(vtrial[0], units)){
    if(units == USER){

      // interpret lines of the form:
      // user 1 -1/3 0 0 1 1.e-10 "Funny units"
      if(vtrial.size() != 8)
	throw Units_Error("Subs::Units::Units(const std::string&): user defined option = " + sunit + " has incorrect number of entries.");

      try{
	mass        = Fraction(vtrial[1]);
	length      = Fraction(vtrial[2]);
	time        = Fraction(vtrial[3]);
	charge      = Fraction(vtrial[4]);
	solid_angle = Fraction(vtrial[5]);
      }
      catch(const std::string& err){
	throw Units_Error("Subs::Units::Units(const std::string&): failed to read fractional exponent");
      }
      std::istringstream sstr(vtrial[6]);
      sstr >> scale;
      if(!sstr)
	throw Units_Error("Subs::Units::Units(const std::string&): failed to read scale factor = " + vtrial[6]);
      
      if(Subs::toupper(vtrial[7]) == "UNDEFINED"){
	pgplot_name = "";
	pgplot_name_set = false;
      }else{
	pgplot_name = vtrial[7];
	pgplot_name_set = true;
      }
      
    }else{
      *this = Units(units);
    }
  }else{
    throw Units_Error("Subs::Units::Units(const std::string&): std::string = " + sunit + " not interpretable.");
  }
}

/** This constructor provides a quick means to set units from an enumerated type
 * \param units the enum value
 * \param name blank for the program to choose it.
 */
Subs::Units::Units(UNITS units, const std::string& name){

  mass        = 0;
  length      = 0;
  time        = 0;
  charge      = 0;
  solid_angle = 0;

  if(name != ""){
    pgplot_name = name;
    pgplot_name_set = true;
  }else{
    pgplot_name_set = false;
  }

  switch(units){
  case ANGSTROM:
    length = 1;
    scale  = 1.e-10;
    if(!pgplot_name_set) pgplot_name = "\\A";
    break;
  case NANNOMETRE:
    length = 1;
    scale  = 1.e-9;
    if(!pgplot_name_set) pgplot_name = "nm";
    break;
  case MICRON:
    length = 1;
    scale  = 1.e-6;
    if(!pgplot_name_set) pgplot_name = "\\gmm";
    break;
  case CENTIMETRE:
    length = 1;
    scale  = 1.e-2;
    if(!pgplot_name_set) pgplot_name = "cm";
    break;
  case METRE:
    length = 1;
    scale  = 1.;
    if(!pgplot_name_set) pgplot_name = "m";
    break;
  case KILOMETRE:
    length = 1;
    scale  = 1.e3;
    if(!pgplot_name_set) pgplot_name = "km";
    break;
  case SOLAR_RADIUS:
    length = 1;
    scale  = Constants::RSUN;
    if(!pgplot_name_set) pgplot_name = "R\\d\(2281)\\u";
    break;
  case AU:
    length = 1;
    scale  = Constants::AU;
    if(!pgplot_name_set) pgplot_name = "AU";
    break;
  case PARSEC:
    length = 1;
    scale  = Constants::PC;
    if(!pgplot_name_set) pgplot_name = "pc";
    break;
  case KILOPARSEC:
    length = 1;
    scale  = 1000.*Constants::PC;
    if(!pgplot_name_set) pgplot_name = "kpc";
    break;
  case MEGAPARSEC:
    length = 1;
    scale  = 1000000.*Constants::PC;
    if(!pgplot_name_set) pgplot_name = "Mpc";
    break;
  case GRAM:
    mass   = 1;
    scale  = 1.e-3;
    if(!pgplot_name_set) pgplot_name = "g";
    break;
  case KILOGRAM:
    mass   = 1;
    scale  = 1.;
    if(!pgplot_name_set) pgplot_name = "kg";
    break;
  case SOLAR_MASS:
    mass   = 1;
    scale  = Constants::MSUN;
    if(!pgplot_name_set) pgplot_name = "M\\d\\(2281)\\u";
    break;
  case SECOND:
    time   = 1;
    scale  = 1.;
    if(!pgplot_name_set) pgplot_name = "s";
    break;
  case MILLISECOND:
    time   = 1;
    scale  = 1.e-3;
    if(!pgplot_name_set) pgplot_name = "10\\u-3\\d s";
    break;
  case MICROSECOND:
    time   = 1;
    scale  = 1.e-6;
    if(!pgplot_name_set) pgplot_name = "10\\u-6\\d s";
    break;
  case MINUTE:
    time   = 1;
    scale  = 60.;
    if(!pgplot_name_set) pgplot_name = "min";
    break;
  case HOUR:
    time   = 1;
    scale  = 3600.;
    if(!pgplot_name_set) pgplot_name = "hr";
    break;
  case DAY:
    time   = 1;
    scale  = Constants::DAY;
    if(!pgplot_name_set) pgplot_name = "day";
    break;
  case YEAR:
    time   = 1;
    scale  = Constants::YEAR;
    if(!pgplot_name_set) pgplot_name = "yr";
    break;
  case FLAMBDA_A:
    mass   =  1;
    length = -1;
    time   = -3;
    scale  = 1e7;
    if(!pgplot_name_set) pgplot_name = "ergs s\\u-1\\d cm\\u-2\\d A\\u-1\\d";
    break;
  case FLAMBDA_NM:
    mass   =  1;
    length = -1;
    time   = -3;
    scale  = 1e6;
    if(!pgplot_name_set) pgplot_name = "ergs s\\u-1\\d cm\\u-2\\d nm\\u-1\\d";
    break;
  case FNU_MJY:
    mass   =  1;
    time   = -2;
    scale  = 1e-29;
    if(!pgplot_name_set) pgplot_name = "mJy";
    break;
  case AB_MAG:
    scale  = 1;
    if(!pgplot_name_set) pgplot_name = "AB";
    break;
  case MPS:
    length =  1;
    time   = -1;
    scale  =  1;
    if(!pgplot_name_set) pgplot_name = "m s\\u-1\\d";
    break;
  case KMS:
    length =  1;
    time   = -1;
    scale  = 1e3;
    if(!pgplot_name_set) pgplot_name = "km s\\u-1\\d";
    break;
  case HERTZ:
    length =  0;
    time   = -1;
    scale  = 1.;
    if(!pgplot_name_set) pgplot_name = "Hz";
    break;
  case CYCPERDAY:
    length =  0;
    time   = -1;
    scale  = 1./Constants::DAY;
    if(!pgplot_name_set) pgplot_name = "cycles/day";
    break;
  case PERSEC:
    length =  0;
    time   = -1;
    scale  = 1.;
    if(!pgplot_name_set) pgplot_name = "sec\\u-1\\d";
    break;
  case PIXEL:
    scale  = 1.;
    if(!pgplot_name_set) pgplot_name = "pixels";
    break;
  case PHASE:
    scale  = 1.;
    if(!pgplot_name_set) pgplot_name = "cycles";
    break;
  case NONE:
    scale  = 1.;
    break;
  default:
    throw Units_Error("Subs::Units: unrecognised units in enumerated constructor");
  }
  pgplot_name_set = true;
}

/** This function comes back with true or false according to whether two
 * units are equivalent. They are equivalent if all the exponents of the
 * fundamental quantities match up. Scaling factors need not be the same.
 * i.e. m/s and km/s are regarded as the same.
 * \param unit1 the first unit
 * \param unit2 the second unit
 * \return true if the two units are equivalent
 */
bool Subs::operator==(const Subs::Units& unit1, const Subs::Units& unit2){

  return (
	  unit1.mass        == unit2.mass        &&
	  unit1.length      == unit2.length      &&
	  unit1.time        == unit2.time        &&
	  unit1.charge      == unit2.charge      &&
	  unit1.solid_angle == unit2.solid_angle
	  );
}

/** ASCII output of Units
 * \param ostr the output stream
 * \param unit the Units
 */
std::ostream& Subs::operator<<(std::ostream& ostr, const Subs::Units& unit){
  ostr << "mltcs = " << unit.mass 
       << ", " << unit.length
       << ", " << unit.time
       << ", " << unit.charge
       << ", " << unit.solid_angle;
  if(unit.pgplot_name_set){
    ostr << ", pgplot name = " << unit.pgplot_name;
  }else{
    ostr << ", pgplot name not set";
  }
  ostr << ", scale = " << unit.scale;
  return ostr;
}

/** This function comes back with true or false according to whether two
 * units are not equivalent. They are not equivalent if any of the exponents of the
 * fundamental quantities fail to match. 
 * \param unit1 the first unit
 * \param unit2 the second unit
 * \return true if the two units are not equivalent
 */
bool Subs::operator!=(const Subs::Units& unit1, const Subs::Units& unit2){
  return !operator==(unit1, unit2);
}

/** Divides two units. This is always possible.
 * \param unit1 the first unit
 * \param unit2 the second unit
 * \return the new unit
 */
Subs::Units Subs::operator/(const Subs::Units& unit1, const Subs::Units& unit2){

  Units new_unit;

  new_unit.length      = unit1.length      - unit2.length;
  new_unit.mass        = unit1.mass        - unit2.mass;
  new_unit.time        = unit1.time        - unit2.time;
  new_unit.charge      = unit1.charge      - unit2.charge ;
  new_unit.solid_angle = unit1.solid_angle - unit2.solid_angle;
  if(unit2.scale == 0)
    throw Units::Units_Error("Subs::operator(Subs::Units&, Subs::Units&): zero scaling factor in units");
  new_unit.scale       = unit1.scale / unit2.scale;

  new_unit.pgplot_name_set = false;
  return new_unit;
}

/** Inverts units while changing scale. This allows a 
 * statement of the form new = 3600 / old. If 'old' were in seconds
 * then the new would be in cycles per hour in this case.
 * \param factor the factor to change the scale by
 * \param units the units to invert
 */
Subs::Units Subs::operator/(double factor, const Subs::Units& unit){
  Units new_unit;

  new_unit.length      = - unit.length;
  new_unit.mass        = - unit.mass;
  new_unit.time        = - unit.time;
  new_unit.charge      = - unit.charge ;
  new_unit.solid_angle = - unit.solid_angle;
  if(unit.scale == 0)
    throw Units::Units_Error("Subs::operator(double, Subs::Units&): zero scaling factor in units");
  new_unit.scale       = 1 / unit.scale / factor;

  new_unit.pgplot_name_set = false;
  return new_unit;
}

/** Multiplies two units. This is always possible.
 * \param unit1 the first unit
 * \param unit2 the second unit
 * \return the new unit
 */
Subs::Units Subs::operator*(const Subs::Units& unit1, const Subs::Units& unit2){

  Units new_unit;

  new_unit.length      = unit1.length      + unit2.length;
  new_unit.mass        = unit1.mass        + unit2.mass;
  new_unit.time        = unit1.time        + unit2.time;
  new_unit.charge      = unit1.charge      + unit2.charge ;
  new_unit.solid_angle = unit1.solid_angle + unit2.solid_angle;

  new_unit.scale       = unit1.scale * unit2.scale;
  new_unit.pgplot_name_set = false;
  return new_unit;
}

/** When combining units you may want to translate one set of numbers to have the same 
 * units as another. This may imply a rescaling. This routine returns the scale factor
 * needed to do this.
 * \param units the units to convert to 
 * \return the scale factor which translates quantities of the current type
 * to the new units. The units must be equivalent of course.
 */
double Subs::Units::get_scale_factor(const Subs::Units& units) const {
  if(units.scale == 0)
    throw Units::Units_Error("double Subs::Units::get_scale_factor(Subs::Units&): new units have zero scaling factor");
  if(*this != units)
    throw Units::Units_Error("double Subs::Units::get_scale_factor(Subs::Units&): new units not compatible with old");

  return scale / units.scale;
}

/** When combining units you may want to translate one set of numbers to have the same 
 * units as another. This may imply a rescaling. This routine returns the scale factor
 * needed to convert one set of units to another
 * \param units1 the units to convert from 
 * \param units2 the units to convert to 
 * \return the scale factor which translates quantities of the current type
 * to the new units. The units must be equivalent of course.
 */
double Subs::get_scale_factor(const Subs::Units& units1, const Subs::Units& units2) {
  if(units2.scale == 0)
    throw Units::Units_Error("double Subs::get_scale_factor(Subs::Units&, Subs::Units&): second units have zero scaling factor");
  if(units1 != units2)
    throw Units::Units_Error("double Subs::get_scale_factor(Subs::Units&, Subs::Units&): units are not compatible with each other");

  return units1.scale / units2.scale;
}


/** This function sets the name to a form suited to plotting with
 * the PGPLOT package. It attempts to recognise simple forms of units,
 * but failing that it will come back with a full string. As a short-cut,
 * nothing is done if an internal flag indicates that the name is already set.
 */
void Subs::Units::set_pgplot_name(){

  if(pgplot_name_set)
    return;
  else
    pgplot_name_set = true;

  if(mass == 1 && length == 0 && time == 0 && charge  == 0 && solid_angle == 0){

    // Mass
    if(fabs(scale-1.) < 1.e-10){
      pgplot_name = "kg";
    }else if(fabs(scale-1.e-3) < 1.e-13){
      pgplot_name = "g";
    }else if(fabs(scale-Constants::MSUN) < 1.e-10*Constants::MSUN){
      pgplot_name = "M\\d\\(2281)\\u";
    }else{
      pgplot_name_set = false;
    }

  }else if(mass == 0 && length == 1 && time == 0 && charge  == 0 && solid_angle == 0){

    // Length
    if(fabs(scale-1.) < 1.e-10){
      pgplot_name = "m";
    }else if(fabs(scale-1.e-2) < 1.e-12){
      pgplot_name = "cm";
    }else if(fabs(scale-1.e3) < 1.e-7){
      pgplot_name = "km";
    }else if(fabs(scale-Constants::AU) < 1.e-10*Constants::AU){
      pgplot_name = "AU";
    }else if(fabs(scale-Constants::PC) < 1.e-10*Constants::PC){
      pgplot_name = "pc";
    }else if(fabs(scale-Constants::RSUN) < 1.e-10*Constants::RSUN){
      pgplot_name = "R\\d\\(2281)\\u";
    }else{
      pgplot_name_set = false;
    }    

  }else if(mass == 0 && length == 0 && time == 1 && charge  == 0 && solid_angle == 0){

    // Time
    if(fabs(scale-1.) < 1.e-10){
      pgplot_name = "s";
    }else if(fabs(scale-Constants::DAY) < 1.e-10*Constants::DAY){
      pgplot_name = "day";
    }else if(fabs(scale-Constants::MINUTE) < 1.e-10*Constants::MINUTE){
      pgplot_name = "min";
    }else if(fabs(scale-Constants::YEAR) < 1.e-10*Constants::YEAR){
      pgplot_name = "yr";
    }else if(fabs(scale-Constants::HOUR) < 1.e-10*Constants::HOUR){
      pgplot_name = "hr";
    }else{
      pgplot_name_set = false;
    }    

  }else if(mass == 0 && length == 0 && time == -1 && charge  == 0 && solid_angle == 0){

    // Frequency
    if(fabs(scale-1.) < 1.e-10){
      pgplot_name = "s\\u-1\\d";
    }else if(fabs(scale-1./Constants::DAY) < 1.e-10/Constants::DAY){
      pgplot_name = "day\\u-1\\d";
    }else if(fabs(scale-1./Constants::MINUTE) < 1.e-10/Constants::MINUTE){
      pgplot_name = "min\\u-1\\d";
    }else if(fabs(scale-1./Constants::YEAR) < 1.e-10/Constants::YEAR){
      pgplot_name = "yr\\u-1\\d";
    }else if(fabs(scale-1./Constants::HOUR) < 1.e-10/Constants::HOUR){
      pgplot_name = "hr\\u-1\\d";
    }else{
      pgplot_name_set = false;
    }    

  }else if(mass == 1 && length == 0 && time == -3 && charge  == 0 && solid_angle == 0){

    // Flux
    if(fabs(scale-1.) < 1.e-10){
      pgplot_name = "W m\\u-2\\d";
    }else if(fabs(scale-1.e-3) < 1.e-13){
      pgplot_name = "ergs s\\u-1\\d cm\\u-2\\d";
    }else{
      pgplot_name_set = false;
    }    

  }else if(mass == 1 && length == -1 && time == -3 && charge  == 0 && solid_angle == 0){

    // Flux density (flambda)
    if(fabs(scale-1.) < 1.e-10){
      pgplot_name = "W m\\u-2\\d m\\u-1\\d";
    }else if(fabs(scale-1.e9) < 1.e-1){
      pgplot_name = "W m\\u-2\\d nm\\u-1\\d";
    }else if(fabs(scale-1.e6) < 1.e-1){
      pgplot_name = "W m\\u-2\\d \\gmm\\u-1\\d";
    }else if(fabs(scale-1.e10) < 1){
      pgplot_name = "W m\\u-2\\d A\\u-1\\d";
    }else if(fabs(scale-1.e7) < 1.e-3){
      pgplot_name = "ergs s\\u-1\\d cm\\u-2\\d A\\u-1\\d";
    }else{
      pgplot_name_set = false;
    }    

  }else if(mass == 1 && length == 0 && time == -2 && charge  == 0 && solid_angle == 0){

    // Flux density (fnu)
    if(fabs(scale-1.) < 1.e-10){
      pgplot_name = "W m\\u-2\\d Hz\\u-1\\d";
    }else if(fabs(scale-1.e-26) < 1.e-1){
      pgplot_name = "mJy";
    }else{
      pgplot_name_set = false;
    }    

  }

  if(!pgplot_name_set){

    // If we are here, we have failed to recognize the name which must 
    // therefore be built up the hard way
    if(fabs(scale-1.) > 1.e-10) 
      pgplot_name = Subs::str(scale) + " ";
    else
      pgplot_name = "";
    
    bool space = false;
    if(mass > 0){
      space = true;
      pgplot_name += "kg\\u" + mass.get_string()        + "\\d";
    }
    if(length > 0){
      if(space) pgplot_name += " ";
      space = true;
      pgplot_name += "m\\u"  + length.get_string()      + "\\d";
    }
    if(time > 0){
      if(space) pgplot_name += " ";
      space = true;
      pgplot_name += "s\\u"  + time.get_string()        + "\\d";
    }
    if(charge > 0){
      if(space) pgplot_name += " ";
      space = true;
      pgplot_name += "C\\u"  + charge.get_string()      + "\\d";
    }
    if(solid_angle > 0){
      if(space) pgplot_name += " ";
      space = true;
      pgplot_name += "sr\\u" + solid_angle.get_string() + "\\d";
    }
    
    if(mass < 0){
      if(space) pgplot_name += " ";
      space = true;
      pgplot_name += "kg\\u" + mass.get_string()        + "\\d";
    }
    if(length < 0){
      if(space) pgplot_name += " ";
      space = true;
      pgplot_name += "m\\u"  + length.get_string()      + "\\d";
    }
    if(time < 0){
      if(space) pgplot_name += " ";
      space = true;
      pgplot_name += "s\\u"  + time.get_string()        + "\\d";
    }
    if(charge < 0){
      if(space) pgplot_name += " ";
      space = true;
      pgplot_name += "C\\u"  + charge.get_string()      + "\\d";
    }
    if(solid_angle < 0){
      if(space) pgplot_name += " ";
      space = true;
      pgplot_name += "sr\\u" + solid_angle.get_string() + "\\d";
    }
    pgplot_name_set = true;
  }
}


/** Writes Units in binary format
 * \param ostr the output stream
 */
void Subs::Units::write(std::ofstream& ostr) const {
  mass.write(ostr);
  length.write(ostr);
  time.write(ostr);
  charge.write(ostr);
  solid_angle.write(ostr);

  ostr.write((char*)&scale,   sizeof(double));  
  write_string(ostr, pgplot_name);
  ostr.write((char*)&pgplot_name_set, sizeof(bool));  

  if(!ostr) throw Units_Error("Subs::Units::write(std::ofstream&): failure during write");
}

/** Reads a Units in binary format
 * \param istr the input stream
 */
void Subs::Units::read(std::ifstream& istr, bool swap_bytes) {
  mass.read(istr, swap_bytes);
  length.read(istr, swap_bytes);
  time.read(istr, swap_bytes);
  charge.read(istr, swap_bytes);
  solid_angle.read(istr, swap_bytes);
  istr.read((char*)&scale,   sizeof(double));
  if(swap_bytes) scale = Subs::byte_swap(scale);
  read_string(istr, pgplot_name, swap_bytes);
  istr.read((char*)&pgplot_name_set, sizeof(bool));  
  if(!istr) throw Units_Error("Subs::Units::read(std::ifstream&): failure during read");
}

/** Skips a Units in binary format
 * \param istr the input stream
 */
void Subs::Units::skip(std::ifstream& istr, bool swap_bytes) {
  Fraction::skip(istr);
  Fraction::skip(istr);
  Fraction::skip(istr);
  Fraction::skip(istr);
  Fraction::skip(istr);
  istr.ignore(sizeof(double));  
  skip_string(istr, swap_bytes);
  istr.ignore(sizeof(bool));  
  if(!istr) throw Units_Error("Subs::Units::skip(std::ifstream&): failure during skip");
}

/** Checks whether any of the fundamental dimensions appear
 * in the object
 */
bool Subs::Units::has_dimensions() const {
  return (mass != 0 || length != 0 || time != 0 || charge != 0 || solid_angle != 0);
}
  
std::string Subs::Units::to_string() const {
  std::ostringstream ostr;
  ostr << "user " << mass << " " << length << " " << time << " " << charge << " " << solid_angle;
  ostr << " " << scale << " \"" << pgplot_name << "\"";
  return ostr.str();
}
