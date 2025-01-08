#include <cmath>
#include <iomanip>
#include <sstream>
#include "trm/subs.h"
#include "trm/telescope.h"
#include "trm/position.h"

/** Constructs a Position of arbitrary RA and Dec (proper motion, parallax and radial velocity
 * all set = 0)
 * \param rah hours of RA
 * \param ram minutes of RA
 * \param ras seconds of RA
 * \param decsgn sign of declination, either '-' or '+'. Must be entered.
 * \param decd degrees of declination
 * \param decm minutes of declination
 * \param decs seconds of declination
 * \param epoch epoch of coordinates
 */

Subs::Position::Position(int rah, int ram, float ras, char decsgn, int decd, int decm, float decs, double epoch){

    r_a  = fabs(double(rah))+fabs(double(ram))/60.+fabs(double(ras))/3600.;
    if(r_a > 24.)
	throw Position_Error("Subs::Position::Position(int, int, float, char, int, int, float, double): RA > 24");
    decl = fabs(double(decd))+fabs(double(decm))/60.+fabs(double(decs))/3600.;
    if(decl > 90.)
	throw Position_Error("Subs::Position::Position(int, int, float, char, int, int, float, double): "
			     "declination out of range -90 to 90");

    if(decsgn == '-'){
	decl = -decl;
    }else if(decsgn != '+'){
	throw Position_Error("Subs::Position::Position(int, int, float, char, int, int, float, double): "
			     "declination sign must equal either '-' or '+'");
    }
    ep = epoch;  
    pm_ra = pm_dec = prlx = r_v = 0.f;
}

/** Constructs a Position of arbitrary RA and Dec (proper motion, parallax and radial velocity
 * all set = 0)
 * \param rah  Right Ascension, decimal hours
 * \param decd Declination, degrees
 * \param epoch epoch of coordinates
 */

Subs::Position::Position(double rah, double decd, double epoch){
    if(rah > 24.) throw Position_Error("Subs::Position::Position(double, double, double): RA > 24");
    if(decd< -90 || decd > 90.) 
	throw Position_Error("Subs::Position::Position(double, double, double): declination out of range -90 to 90");

    r_a   = rah;
    decl  = decd;
    ep    = epoch;
    pm_ra = pm_dec = prlx = r_v = 0.f;
}

Subs::Position::Position(const std::string& pos){
    std::istringstream ist(pos);
    ist >> *this;  
    if(!ist) 
	throw Position_Error("Subs::Position::Position(const std::string&): failed to construct from the std::string = " + pos);
}

void Subs::Position::set(const std::string& pos){
    std::istringstream ist(pos);
    ist >> *this;  
    if(!ist) 
	throw Position_Error("void Subs::Position::set(const std::string&): failed to set from the std::string = " + pos);
}

/** Returns a string of the current position.
 */

std::string Subs::Position::ra_dec() const {
    int   rah  = int(floor(ra()));
    int   ram  = int(floor(60.*(ra()-rah)));
    int   ras  = int(floor(3600.*(ra()-(rah+ram/60.))));
    int   rafs = int(floor(36000.*(ra()-(rah+ram/60.+ras/3600.))+0.5));

    if(rafs == 10){
	ras++;
	rafs = 0;
    }
    if(ras == 60){
	ram++;
	ras = 0;
    }
    if(ram == 60){
	rah++;
	ram = 0;
    }
    if(rah == 24) rah = 0;

    double d = dec();
    char sign = '+';
    if(d < 0.){
	sign = '-';
	d = -d;
    }
    int  decd  = int(floor(d));
    int  decm  = int(floor(60.*(d-double(decd))));
    int  decs  = int(floor(3600.*(d-(decd+decm/60.))));
    int  decfs = int(floor(36000.*(d-(decd+decm/60.+decs/3600.))+0.5));

    if(decfs == 10){
	decs++;
	decfs = 0;
    }
    if(decs == 60){
	decm++;
	decs = 0;
    }
    if(decm == 60){
	decd++;
	decm = 0;
    }

    std::ostringstream ost;
    ost.setf(std::ios::right);
    ost << std::setw(2) << std::setfill('0') << rah << ":" 
	<< std::setw(2) << std::setfill('0') << ram << ":" 
	<< std::setw(2) << std::setfill('0') << ras << "."
	<< std::setw(1) << std::setfill('0') << rafs 
	<< sign << std::setw(2) << std::setfill('0') << decd << ":" 
	<< std::setw(2) << std::setfill('0') << decm << ":" 
	<< std::setw(2) << std::setfill('0') << decs << "."
	<< std::setw(1) << std::setfill('0') << decfs << std::setfill(' ');  
    return ost.str();
}

/** Apply proper motion of the star to update position to the epoch
 *  of another position. It is assumed the coordinates are ICRS. 
 * \param pos the Position defining the new epoch
 */
void Subs::Position::update(const Position& pos){
    update(pos.epoch());
}

/** Apply proper motion / space velocity coorections to update the position from its current 
 *  epoch to a new one. It is assumed the coordinates are ICRS (which makes things much simpler). 
 * \param epch the epoch to update the coordinates to.
 */
void Subs::Position::update(double epch){

    if(epch == epoch()) return;

    // go from julian epoch to MJD
    double _MJD0, MJD_ep1, MJD_ep2;
    // _MJD0 is the MJD zero point at 2400000.5
    iauEpj2jd(epoch(), &_MJD0, &MJD_ep1);
    iauEpj2jd(epch, &_MJD0, &MJD_ep2);

    // Apply proper motion; space velocity to go from epoch() to epch
    double ra, dec, _pmr, _pmd, _px, _rv_;
    int code = 0;
    code = iauStarpm(
        rar(), decr(), pmrar(), pmdecr(), parallax(), rv(),
        MJD0, MJD_ep1, MJD0, MJD_ep2,
        &ra, &dec, &_pmr, &_pmd, &_px, &_rv_);
    
    // codes are non recoverable
    if(code != 0) throw Position_Error("Subs::Position::update(double): error code = " + Subs::str(code));
    // update the position
    set_radec(ra, dec);
    // update the epoch
    epoch() = epch;
}

/** Returns a vector corresponding to Position. This is calculated
 * directly from the ra and dec; if you need to allow for space motion you should
 * first update the position. The frame  is the BCRS.
 */

Subs::Vec3 Subs::Position::vect() const {
    double v[3];
    iauS2c(rar(), decr(), v);
    Subs::Vec3 vec(v);
    return Subs::Vec3(v);
}

/** Returns barycentric light travel time correction. This calculates the position of the observatory
 * defined by 'tel' and the target in the BCRS frame and then takes the dot product divided by the speed of light
 * to get correction that needs to be added to the time observed at Earth to get the time observed at the
 * solar system barycentre.
 * \param time the time in question.
 * \param tel  the relevant telescope
 * \return The routine returns a time in seconds that when added to the local time
 * gives the equivalent time at the barycentre of the Solar system.
 */
double Subs::Position::tcorr_bar(const Time& time, const Telescope& tel) const {

    // Compute observatory vector in BCRS
    Vec3 pobs = time.earth_pos_bar(tel);

    // Update the position to the current epoch
    Position pos = *this;
    pos.update(time.jepoch());

    return dot(pos.vect(), pobs)/Constants::C;
}

/** Returns heliocentric light travel time correction. This calculates the position of the observatory
 * defined by 'tel' and the target in the BCRS frame and then takes the dot product divided by the speed of light
 * to get correction that needs to be added to the time observed at Earth to get the time observed at the
 * heliocentre.
 * \param time the time in question.
 * \return The routine returns a time in seconds that when added to the local time
 * gives the equivalent time at the centre of the Sun
 */
double Subs::Position::tcorr_hel(const Time& time, const Telescope& tel) const {

    // Compute observatory vector in BCRS
    Vec3 pobs = time.earth_pos_hel(tel);

    // Update the position to the current epoch
    Position pos = *this;
    pos.update(time.jepoch());

    return dot(pos.vect(), pobs)/Constants::C;

}

/** Returns helio- and bary-centric time and velocity corrections. This is faster than calling both
 * tcorr_hel and tcorr_bar one by one.
 * \param time the time in question.
 * \param tel  the relevant telescope
 * \return The routine returns a small structure containing time and velocity corrections.
 */
Subs::Pinfo Subs::Position::pinfo(const Time& time, const Telescope& tel) const {

    // Update the position to the current epoch, convert to a vector
    Position pos = *this;
    pos.update(time.jepoch());
    Vec3 ptarg = pos.vect();

    // Compute full panoply of observatory positions and velocities
    Vec3 ph, vh, pb, vb;
    time.earth(tel, ph, vh, pb, vb);

    Pinfo temp;
    temp.tcor_hel    =  dot(ptarg,ph)/Constants::C;
    temp.tcor_bar    =  dot(ptarg,pb)/Constants::C;
    temp.vearth_hel  = -dot(ptarg,vh)/1000;
    temp.vearth_bar  = -dot(ptarg,vb)/1000;

    return temp;
}

/** Sets the position to that of the Sun at the supplied time.
 * \param time the time
 * \param tel  the telescope
 */
void Subs::Position::set_to_sun(const Time& time, const Telescope& tel) {

    Subs::Vec3 velhel, poshel = time.earth_pos_hel(tel);

    // Make into position of the Sun rather than the Earth.
    poshel = -poshel;

    double d[3];
    poshel.get(d);
    double ras, decs;

    // convert vector to spherical
    iauC2s(d,&ras,&decs);

    // Set the internal variables. 
    epoch() = time.jepoch();
    ra()    = 24.*ras/Constants::TWOPI;
    dec()   = 360.*decs/Constants::TWOPI;
    pmra()  = pmdec() = parallax() = rv() = 0.f;
}

/** This computes altitude, azimuth etc. Note that the polar motion X, Y and 
 * UT1-UTC are all set to 0 inside this routine, while observing conditions are 
 * set to default values so it is not to be used for accurate work.
 * \param time the time of interest
 * \param tel  the telescope
 * \return Structure containing the bits and pieces of information.
 */

Subs::Altaz Subs::Position::altaz(const Time& time, const Telescope& tel) const {

    // Apply space motions
    Position pos = *this;
    pos.update(time.jepoch());
    int corr = 0;

    // compute data, using defaults for many values.
    double utc = time.mjd();
    double aob, zob, hob, dob, rob;
    const double T    = 290.;    // ambient temp, K
    const double P    = 1013.25; // ambient pressure, mbar
    const double RH   = 0.2;     // relative humidity (0-1)
    const double WAVE = 0.55;    // observing wavelength microns
    const double TLR  = 0.0065;  // lapse rate, K/metre

    iauAtio13(pos.rar(), pos.decr(), 
              MJD0, utc, 0., 
              tel.longituder(), tel.latituder(), tel.height(), 
              0., 0., P, T-273.15, RH, WAVE, 
              &aob, &zob, &hob, &dob, &rob);
    // put in the expected ra and dec ranges
    // aob -pi to pi
    if (aob < -M_PI) aob += Constants::TWOPI;
    if (aob > M_PI) aob -= Constants::TWOPI;
    // zob 0 to 2pi
    if (zob < 0.) zob += Constants::TWOPI;
    if (zob > Constants::TWOPI) zob -= Constants::TWOPI;
    // hob -pi to pi
    if (hob < -M_PI) hob += Constants::TWOPI;
    if (hob > M_PI) hob -= Constants::TWOPI;
    // dob -pi/2 to pi/2
    double correction = 0.;
    if (dob < -Constants::PI/2.) {
        corr = dob + Constants::PI/2;
        dob = Constants::PI/2 + corr;
    } else
    if (dob > Constants::PI/2) {
        corr = dob - Constants::PI/2;
        dob = Constants::PI/2 - corr;
    }
    // rob 0 to 2pi
    if (rob < 0.) rob += Constants::TWOPI;
    if (rob > Constants::TWOPI) rob -= Constants::TWOPI;

    // compute refraction
    double ref;
    // SOFA refactor
    double refa, refb;
    iauRefco(P, T-273.15, RH, WAVE, &refa, &refb);
    ref = refa * tan(zob) + refb * pow(tan(zob), 3.);
    double zvac = zob + ref;

    Altaz temp;
    temp.ha       = 24.*hob/Constants::TWOPI;
    temp.alt_true = 90.-360.*zvac/Constants::TWOPI;
    temp.alt_obs  = 90.-360.*zob/Constants::TWOPI;
    temp.az       = 360.*aob/Constants::TWOPI;
    // calculate airmass
    double SECZM1 = 1.0 / (cos(std::min(1.52, std::abs(zob)))) - 1.0;
    temp.airmass  = 1.0 + SECZM1 * (0.9981833 - SECZM1 * (0.002875 + 0.0008083 * SECZM1));

    // Compute pa
    temp.pa = 360.*iauHd2pa(hob,pos.decr(),tel.latituder())/Constants::TWOPI;
    temp.pa = temp.pa > 0. ? temp.pa : 360.+temp.pa;
    return temp;
}

/** ASCII output of a position.
 */

std::ostream& Subs::operator<<(std::ostream& ost, const Position& pos){
    int   rah  = int(floor(pos.ra()));
    int   ram  = int(floor(60.*(pos.ra()-rah)));
    int   ras  = int(floor(3600.*(pos.ra()-(rah+ram/60.))));
    int   rafs = int(floor(36000.*(pos.ra()-(rah+ram/60.+ras/3600.))+0.5));

    if(rafs == 10){
	ras++;
	rafs = 0;
    }
    if(ras == 60){
	ram++;
	ras = 0;
    }
    if(ram == 60){
	rah++;
	ram = 0;
    }
    if(rah == 24) rah = 0;

    double d = pos.dec();
    char sign = '+';
    if(d < 0.){
	sign = '-';
	d = -d;
    }
    int  decd  = int(floor(d));
    int  decm  = int(floor(60.*(d-double(decd))));
    int  decs  = int(floor(3600.*(d-(decd+decm/60.))));
    int  decfs = int(floor(36000.*(d-(decd+decm/60.+decs/3600.))+0.5));

    if(decfs == 10){
	decs++;
	decfs = 0;
    }
    if(decs == 60){
	decm++;
	decs = 0;
    }
    if(decm == 60){
	decd++;
	decm = 0;
    }

    ost.setf(std::ios::right);
    ost << "RA: " 
	<< std::setw(2) << std::setfill('0') << rah << ":" 
	<< std::setw(2) << std::setfill('0') << ram << ":" 
	<< std::setw(2) << std::setfill('0') << ras << "."
	<< std::setw(1) << std::setfill('0') << rafs 
	<< ", Dec: " << sign 
	<< std::setw(2) << std::setfill('0') << decd << ":" 
	<< std::setw(2) << std::setfill('0') << decm << ":" 
	<< std::setw(2) << std::setfill('0') << decs << "."
	<< std::setw(1) << std::setfill('0') << decfs << std::setfill(' ')  
	<< ", Ep: " << pos.epoch();

    return ost;
}

/** ASCII input of a position of the form "01:10:12.2 +20:00:34.34", or indeed
 * "01 10 12.2 +20 00 34.34".
 * You can optionally include proper motion etc data if you add a string of the form
 * " P 0.66 -0.45 0.034 -21.3 2000.0"  where the P flags that we are referring to proper motion
 * and there follows the proper motion in arcseconds per year in RA and Dec, the parallax in
 * arcseconds and the radial velocity in km/s and the epoch of the coordinates.
 * \param ist input stream
 * \param pos the position to load into.
 */

std::istream& Subs::operator>>(std::istream& ist, Position& pos){

    if(!ist) return ist;
  
    int rah, ram, decd, decm;
    char decsgn, c;
    double ras, decs;
    double epoch;

    // Read each item, skipping over one character between some of them to allow
    // for colons or blanks
    ist >> rah;
    if(!ist) return ist;
    ist.ignore();
    if(!ist) return ist;
    ist >> ram;
    if(!ist) return ist;
    ist.ignore();
    if(!ist) return ist;
    ist >> ras >> decsgn >> decd; 
    if(!ist) return ist;
    ist.ignore();
    if(!ist) return ist;
    ist >> decm;
    if(!ist) return ist;
    ist.ignore();
    ist >> decs;
    if(!ist) return ist;
    c = ' ';
    while(ist && (c == ' ' || c == '\t')){
	ist.get(c);
    }
    if(ist && c == 'P'){
	ist >> pos.pm_ra >> pos.pm_dec >> pos.prlx >> pos.r_v >> epoch;
	if(!ist) return ist;

    }else{
	if(ist){
	    ist.unget();
	}else if(ist.eof()){
	    ist.clear();
	}
	pos.pm_ra = pos.pm_dec = pos.prlx = pos.r_v = 0;
	epoch = 2000.;
    }

    pos.r_a  = fabs(double(rah))+fabs(double(ram))/60.+fabs(double(ras))/3600.;
    pos.decl = fabs(double(decd))+fabs(double(decm))/60.+fabs(double(decs))/3600.;
    if(decsgn == '-'){
	pos.decl = -pos.decl;
    }else if(decsgn != '+'){
	throw Position::Position_Error("Declination sign error on input");
    }
    pos.ep = epoch;
    return ist;
}

double Subs::dot(const Position& pos1, const Position& pos2){
    if(pos1.epoch() != pos2.epoch()){
	std::cerr << "WARNING: taking dot product of positions with different epoch!\n";
	std::cerr << "         ... no correction is made.\n";
    }
    double rd1 = Constants::TWOPI*pos1.dec()/360.;
    double rd2 = Constants::TWOPI*pos2.dec()/360.;
    return sin(rd1)*sin(rd2)+cos(rd1)*cos(rd2)*
	cos(Constants::TWOPI*(pos1.ra()-pos2.ra())/24.);
}

void Subs::Position::read(std::ifstream& s, bool swap_bytes){
    REAL8 ra, dec, epoch;
    s.read((char*)&ra,sizeof(REAL8));
    if(swap_bytes) ra = Subs::byte_swap(ra);
    if(ra < 0. || ra >= 24.) throw Position_Error("RA out of range in Subs::Position::read(std::ifstream&, bool)");

    s.read((char*)&dec,sizeof(REAL8));
    if(swap_bytes) dec = Subs::byte_swap(dec);
    if(dec < -90. || dec > 90.) throw Position_Error("Dec out of range in Subs::Position::read(std::ifstream&, bool)");

    s.read((char*)&epoch,sizeof(REAL8));
    if(swap_bytes) epoch = Subs::byte_swap(epoch);
    if(epoch < 1900. || epoch > 2100.) throw Position_Error("Epoch out of range in Subs::Position::read(std::ifstream&, bool)");

    REAL4 rapm, decpm, parallax, radvel;
    s.read((char*)&rapm,    sizeof(REAL4));
    s.read((char*)&decpm,   sizeof(REAL4));
    s.read((char*)&parallax,sizeof(REAL4));
    s.read((char*)&radvel,  sizeof(REAL4));
    if(swap_bytes){
	rapm     = Subs::byte_swap(rapm);
	decpm    = Subs::byte_swap(decpm);
	parallax = Subs::byte_swap(parallax);
	radvel   = Subs::byte_swap(radvel);
    }
    if(!s) throw Position_Error("Read error in Subs::Position::read(std::ifstream&, bool)");
    r_a    = ra;
    decl   = dec;
    ep     = epoch;
    pm_ra  = rapm;
    pm_dec = decpm;
    prlx   = parallax;
    r_v    = radvel;
}

void Subs::Position::write(std::ofstream& s) const {
    s.write((char*)&r_a,  sizeof(REAL8));
    s.write((char*)&decl, sizeof(REAL8));
    s.write((char*)&ep,   sizeof(REAL8));
    if(!s) throw Position_Error("void Subs::Position::write(std::ofstream&) const: error writing position");

    s.write((char*)&pm_ra,  sizeof(REAL4));
    s.write((char*)&pm_dec, sizeof(REAL4));
    if(!s) throw Position_Error("void Subs::Position::write(std::ofstream&) const: error writing proper motions");

    s.write((char*)&prlx,  sizeof(REAL4));
    s.write((char*)&r_v,   sizeof(REAL4));
    if(!s) throw Position_Error("void Subs::Position::write(std::ofstream&) const: error writing parallax and radial velocity");

}

void Subs::Position::skip(std::ifstream& s) {
    s.ignore(sizeof(REAL8));
    s.ignore(sizeof(REAL8));
    s.ignore(sizeof(REAL8));
    s.ignore(sizeof(REAL4));
    s.ignore(sizeof(REAL4));
    s.ignore(sizeof(REAL4));
    s.ignore(sizeof(REAL4));
}
