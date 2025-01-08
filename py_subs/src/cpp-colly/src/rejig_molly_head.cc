#include <string>
#include "trm/subs.h"
#include "trm/header.h"
#include "trm/colly.h"

/**
 * rejig_molly_head rejigs a standard molly header read by 
 * read_molly_head, basically making it look nicer, by, for 
 * example, combining positional information into one 
 * directory. You should not use rejig_molly_head if you are 
 * going to use write_molly_head to write out headers.
 *
 * \param head the header to rejig.
 * \param warn true to print warnings of header items not found
 * although expected; false to suppress these warnings
 */

void Colly::rejig_molly_head(Subs::Header& head, bool warn){

  try{
    head.set("Time", new Subs::Hdirectory("Timing information"));
    int day    = head["Day"]->get_int();
    int month  = head["Month"]->get_int();
    int year   = head["Year"]->get_int();
    double utc = head["UTC"]->get_double();
    head.erase("Day");
    head.erase("Month");
    head.erase("Year");
    head.erase("UTC");
    bool add;
    if((add = utc > 24.)) utc -= 24.;
    Subs::Time time(day,Subs::Date::Month(month),year,utc);
    time.add_day(1);
    head.set("Time.UTC", new Subs::Htime(time,"UTC of mid-exposure"));
    head.rename("Dwell","Time.Dwell");
    head["Time.Dwell"]->set_comment("Exposure time (seconds)");
  }
  catch(const std::string& s){
    if(warn) std::cerr << s << std::endl;
  }

#ifdef DEBUG
    std::cerr << "Rejigged times" << std::endl;
#endif

  try{
    head.rename("HJD","Time.HJD");
    head["Time.HJD"]->set_comment("Heliocentric Julian day at "
				  "centre of exposure");

    head.rename("RJD","Time.JD");
    head["Time.JD"]->set_comment("Julian day at centre of exposure");

    head.rename("Delay","Time.JD-HJD");
    head["Time.JD-HJD"]->set_comment("JD minus heliocentric JD");
  }
  catch(const std::string& s){
    if(warn) std::cerr << s << std::endl;
  }

#ifdef DEBUG
  std::cerr << "Rejigged JDs" << std::endl;
#endif

  Subs::Header::Hnode *cp;

  try{
    head.set("Position", new Subs::Hdirectory("Positional information"));
    double ra, dec, eq;

    // Allow for a bit of variation in the names here
    
    if((cp = head.find("RA"))->has_data()){
	ra   = cp->value->get_double();
	head.erase(cp);
    }else{
	throw Colly_Error("Colly::rejig_molly_head: No RA present");
    }

    if((cp = head.find("Dec"))->has_data() || 
       (cp = head.find("DEC"))->has_data()){
	dec  = cp->value->get_double();
	head.erase(cp);
    }else{
	throw Colly_Error("Colly::rejig_molly_head: No Dec present");
    }
    
    // I have some spectra where the Equinox is mistakenly
    // called Epoch

    if((cp = head.find("Equinox"))->has_data() || 
       (cp = head.find("EQUINOX"))->has_data() ||
       (cp = head.find("Epoch"))->has_data()){
	eq   = cp->value->get_double();
	head.erase(cp);
    }else{
	throw Colly_Error("Colly::rejig_molly_head: No Equinox present");
    }

    head.set("Position.Celestial", new Subs::Hposition(Subs::Position(ra,dec,eq),"Celestial coordinates"));
    head.set("Position.Galactic",  new Subs::Hdirectory("Galactic coordinates"));
    head.set("Position.Horizon",   new Subs::Hdirectory("Horizon coordinates"));

    head.rename("Gal_latitude","Position.Galactic.Latitude");
    head["Position.Galactic.Latitude"]->set_comment("Galactic latitude (degrees)");

    head.rename("Gal_longitude","Position.Galactic.Longitude");
    head["Position.Galactic.Latitude"]->set_comment("Galactic longitude (degrees)");

    head.rename("Sidereal_time","Position.Horizon.Sidereal");
    head["Position.Horizon.Sidereal"]->set_comment("Local sidereal time (hours)");

    head.rename("Hour_angle","Position.Horizon.HA");
    head["Position.Horizon.HA"]->set_comment("Hour angle");

    head.rename("Airmass","Position.Horizon.Airmass");
    head["Position.Horizon.Airmass"]->set_comment("Airmass");
  }
  catch(const std::string& s){
    if(warn) std::cerr << s << std::endl;
  }

#ifdef DEBUG
    std::cerr << "Rejigged position" << std::endl;
#endif

  try{
    head.set("Observatory", new Subs::Hdirectory("Observatory & telescope information"));

    head.rename("Site","Observatory.Site");
    head.rename("Telescope","Observatory.Telescope");

    head.rename("Longitude","Observatory.Longitude");
    head["Observatory.Longitude"]->set_comment("Degrees north");

    head.rename("Latitude","Observatory.Latitude");
    head["Observatory.Longitude"]->set_comment("Degrees west");

    head.rename("Height","Observatory.Height");
    head["Observatory.Height"]->set_comment("Height in metres");
  }
  catch(const std::string& s){
    if(warn) std::cerr << s << std::endl;
  }

#ifdef DEBUG
    std::cerr << "Rejigged observatory" << std::endl;
#endif

  try{
    head.set("Miscellaneous", new Subs::Hdirectory("Miscellaneous information"));

    head.rename("Record","Miscellaneous.Record");
    head["Miscellaneous.Record"]->set_comment("Record number");

    head.rename("Vearth","Miscellaneous.Vearth");
    head["Miscellaneous.Vearth"]->set_comment("Apparent velocity of target"
                                              " due to the Earth's motion");

    head.rename("Night","Miscellaneous.Night");
    head["Miscellaneous.Night"]->set_comment("Grey scattering factor on top of"
	       			             " Rayleigh scattering");

    head.rename("Extract_position","Miscellaneous.Extract");
    head["Miscellaneous.Extract"]->set_comment("Spatial position of extraction");

    head.rename("Arc(s)_used","Miscellaneous.Arcs");
    head["Miscellaneous.Arcs"]->set_comment("Record number(s) of arc(s) used"
					    " for wavelength calibration");

    head.rename("Rayleigh","Miscellaneous.Rayleigh");
    head["Miscellaneous.Rayleigh"]->set_comment("Scale factor for Rayleigh"
						" scattering used");

    head.rename("Grey","Miscellaneous.Grey");
    head["Miscellaneous.Grey"]->set_comment("Grey scattering factor on top of"
					    " Rayleigh scattering");

    head.rename("Flux_star","Miscellaneous.Flux");
    head["Miscellaneous.Flux"]->set_comment("Record number of flux star used"
					    " for flux calibration");
  }
  catch(const std::string& s){
    if(warn) std::cerr << s << std::endl;
  }

#ifdef DEBUG
    std::cerr << "Rejigged miscellaneous" << std::endl;
#endif

}


  

  






