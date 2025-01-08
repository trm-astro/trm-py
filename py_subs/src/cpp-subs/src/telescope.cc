#include <iostream>
#include "trm/subs.h"
#include "trm/telescope.h"

Subs::Telescope::Telescope(const std::string& tname, const std::string& sname,
			   int longd, int longm, float longs, char eorw, 
			   int latd, int latm, float lats, char sorn, 
			   float height) : 
  telescope_name(tname), site_name(sname) {
  
  lng = fabs(double(longd))+fabs(double(longm))/60.+fabs(double(longs))/3600.;
  if(eorw == 'w' || eorw == 'W'){
    lng = - lng;
  }else if(eorw != 'e' && eorw != 'E'){
    throw Telescope_Error("Subs::Telescope::Telescope(const std::string&, const std::string&,"
			  "int, int, float, char, int, int, float, char, float): " 
			  " longitude sign must be one of e, E, w, or W");
  }
  if(lng < -180. || lng > +180.)
    throw Telescope_Error("Subs::Telescope::Telescope(const std::string&, const std::string&,"
			  "int, int, float, char, int, int, float, char, float): " 
			  " longitude must lie in range -180 to +180");  
  
  lat  = fabs(double(latd))+fabs(double(latm))/60.+fabs(double(lats))/3600.;
  if(sorn == 's' || sorn == 'S' || sorn == '-'){
    lat = - lat;
  }else if(sorn != 'n' && sorn != 'N' && sorn != '+'){
    throw Telescope_Error("Subs::Telescope::Telescope(const std::string&, const std::string&,"
			  "int, int, float, char, int, int, float, char, float): " 
			  " latitude sign must be one of s, S, -, n, N or +");
  }
  if(lat < -90. || lat > +90.)
    throw Telescope_Error("Subs::Telescope::Telescope(const std::string&, const std::string&,"
			  "int, int, float, char, int, int, float, char, float): " 
			  " latitude must lie in range -90 to +90");

  hgt = height;
}

Subs::Telescope::Telescope(const std::string& name){
  this->set(name);
}

void Subs::Telescope::set(const std::string& name){

  if(Subs::toupper(name) == "WHT"){
    *this = Telescope("WHT","La Palma",17,52,53.9,'W',28,45,38.3,'N',2332.);

  }else if(Subs::toupper(name) == "INT"){
    *this = Telescope("INT","La Palma",17,52,40.,'W',28,45,43.,'N',2336.);

  }else if(Subs::toupper(name) == "SOUTHAMPTON"){
    *this = Telescope("Southampton", "Southampton",01,00,00,'W',51,00,00,'N', 20.);

  }else if(Subs::toupper(name) == "SAAO"){
    *this = Telescope("SAAO", "Sutherland",20,48,38.5,'E',32,22,46.,'S', 1798.);

  }else if(Subs::toupper(name) == "NTT"){
    *this = Telescope("NTT", "La Scilla",70,44,00.,'W',29,15,00.,'S', 2400.);

  }else if(Subs::toupper(name) == "VLT"){
    *this = Telescope("VLT", "Paranal",70,24,9.9,'W',24,37,30.3,'S',2635.);

  }else if(Subs::toupper(name) == "MAUI"){
    *this = Telescope("Maui", "Faulkes, Maui",156,15,29.,'W',20,42,30.,'N', 3058.);

  }else if(Subs::toupper(name) == "UKIRT"){
    *this = Telescope("UKIRT", "Mauna Kea", 155,28,13.,'W',19,49,21.,'N', 4200.);

  }else if(Subs::toupper(name) == "CALAR ALTO"){
    *this = Telescope("Calar Alto", "Spain",02,32,46.5,'W',37,13,25.,'N',2168.);

  }else if(Subs::toupper(name) == "WROCLAW"){
    *this = Telescope("Wroclaw", "Poland",17,05,00.,'E',51,07,00.,'N',1000.);

  }else if(Subs::toupper(name) == "JODRELL"){
    *this = Telescope("Jodrell", "Jodrell Bank",02,18,25.7,'W',53,14,10.5,'N',77.);

  }else if(Subs::toupper(name) == "TNT"){
      *this = Telescope("TNT", "Doi Inthanon",98,28,00.0,'E',18,34,00.0,'N',2457.);

  }else{
    throw Telescope_Error("Name: \"" + name + 
			  "\" unrecognised in Telescope::set(const std::string&).\n\n"
			  "Recognised telescope names (case insensitive) are:\n\n"
			  " 1) WHT\n"
			  " 2) INT\n"
			  " 3) Southampton\n"
			  " 4) SAAO\n"
			  " 5) NTT\n"
			  " 6) VLT\n"
			  " 7) Maui\n"
			  " 8) UKIRT\n"
			  " 9) Calar Alto\n"
			  "10) Wroclaw\n"
			  "11) Jodrell\n"
			  "12) TNT\n"
			  );
  }
}  

std::ostream& Subs::operator<<(std::ostream& ost, const Telescope& obs){
  double lng = obs.longitude();
  std::string eorw = "East";
  if(lng < 0){
    eorw = "West";
    lng = -lng;
  }
  double lat = obs.latitude();
  std::string sorn = "North";
  if(lat < 0){
    sorn = "South";
    lat = -lat;
  }
  ost << "Telescope: " << obs.telescope() << ", site: " << obs.site() 
      << ", " << lng << " deg. " << eorw << ", " << lat << " deg. " 
      << sorn << ", height = " << obs.height() << "m";
  return ost;
}

bool Subs::operator!=(const Subs::Telescope& tel1, const Subs::Telescope& tel2) {
  return (tel1.lng != tel2.lng || tel1.lat != tel2.lat || tel1.hgt != tel2.hgt);
} 

void Subs::Telescope::read(std::ifstream& s, bool swap_bytes){

  std::string tname, sname;
  read_string(s, tname, swap_bytes);
  read_string(s, sname, swap_bytes);
  if(!s) throw Telescope_Error("Subs::Telescope::read(std::ifstream&): error reading names");

  REAL8 longit, latit;
  s.read((char*)&longit,sizeof(REAL8));
  if(swap_bytes) longit = Subs::byte_swap(longit);
  if(s && (longit < -360. || longit >= 360.)) 
    throw Telescope_Error("Subs::Telescope::read(std::ifstream&): longitude out of range -360 to +360");

  s.read((char*)&latit,sizeof(REAL8));
  if(swap_bytes) latit = Subs::byte_swap(latit);
  if(s && (latit < -90. || latit > 90.)) 
    throw Telescope_Error("Subs::Telescope::read(std::ifstream&): latitude out of range -90 to +90");

  REAL4 height;
  s.read((char*)&height,sizeof(REAL4));
  if(swap_bytes) height = Subs::byte_swap(height);

  if(!s) throw Telescope_Error("Read error in Subs::Telescope::read(std::ifstream&)");

  telescope_name = tname;
  site_name = sname;
  lng    = longit;
  lat    = latit;
  hgt    = height;
}

void Subs::Telescope::write(std::ofstream& s) const {

  write_string(s, telescope_name);
  write_string(s, site_name);
  if(!s) throw Telescope_Error("void Subs::Telescope::write(std::ofstream&) const: error writing names");

  s.write((char*)&lng, sizeof(REAL8));
  s.write((char*)&lat, sizeof(REAL8));
  s.write((char*)&hgt, sizeof(REAL4));
  if(!s) throw Telescope_Error("void Subs::Telescope::write(std::ofstream&) const: error writing position");

}

void Subs::Telescope::skip(std::ifstream& s, bool swap_bytes) {
  skip_string(s,swap_bytes); 
  skip_string(s,swap_bytes); 
  s.ignore(sizeof(REAL8));
  s.ignore(sizeof(REAL8));
  s.ignore(sizeof(REAL4));
}


/** ASCII input of a telescope in the form 
 * Telescope name
 * Site name
 * 120.67672 -34.787682 2345
 * \param ist input stream
 * \param pos the position to load into.
*/

std::istream& Subs::operator>>(std::istream& ist, Telescope& tel){

  std::string tname, sname;
  double longit, latit;
  float height;

  // Read each item
  if(!ist) return ist;
  ist >> tname;
  if(!ist) return ist;
  ist >> sname;
  if(!ist) return ist;
  ist >> longit >> latit >> height; 
  if(!ist) return ist;
  if(longit < -360. || longit >= 360.) 
    throw Telescope::Telescope_Error("Subs::operator>>(std::istream&, Telescope&): longitude out of range -360 to +360");
  if(latit < -90. || latit > 90.) 
    throw Telescope::Telescope_Error("Subs::operator>>(std::istream&, Telescope&): latitude out of range -90 to +90");
  
  tel.telescope_name = tname;
  tel.site_name = sname;
  tel.lng  = longit;
  tel.lat  = latit;
  tel.hgt  = height;

  return ist;
}

