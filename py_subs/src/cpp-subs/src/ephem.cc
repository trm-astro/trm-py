#include <iomanip>
#include "trm/subs.h"
#include "trm/ephem.h"

Subs::Ephem::Ephem(double T0, double period, TSCALE tscale) 
  : etype(LINEAR), tscale(tscale), tzero(T0), per(period), tzerr(0.), perr(0.) {}

Subs::Ephem::Ephem(double T0, double period, double T0err, double pererr, TSCALE tscale) : 
  etype(LINEAR), tscale(tscale), tzero(T0), per(period), tzerr(T0err), perr(pererr) {}

Subs::Ephem::Ephem(double T0, double period, double pdot, TSCALE tscale) : 
  etype(QUADRATIC), tscale(tscale), tzero(T0), per(period), quad(pdot), tzerr(0.), perr(0.), qerr(0.) {}

Subs::Ephem::Ephem(double T0, double period, double pdot, double T0err, double pererr, double qderr, TSCALE tscale) : 
  etype(QUADRATIC), tscale(tscale), tzero(T0), per(period), quad(pdot), tzerr(T0err), perr(pererr), qerr(qderr) {}

double Subs::Ephem::phase(double t) const {
  double ph = (t-tzero)/per;
  if(etype == QUADRATIC){
    double pold = ph+1.;
    while(fabs(ph-pold) > 1.e-6){
      pold = ph;
      ph   = (t-tzero-quad*ph*ph)/per;
    }
  }
  return ph;
}

double Subs::Ephem::pherr(double t) const {
  double ph   = phase(t);
  double terr = sqr(tzerr) + sqr(perr*ph);
  if(etype == QUADRATIC) terr += sqr(qerr*ph*ph);
  return sqrt(terr)/per;
}

double Subs::Ephem::time(double p) const {
  if(etype == QUADRATIC){
    return tzero+p*(per+quad*p);
  }else{
    return tzero+p*per;
  }
}

double Subs::Ephem::timerr(double p) const {
  if(etype == QUADRATIC){
    return sqrt(Subs::sqr(tzerr) + Subs::sqr(p*perr) + Subs::sqr(qerr*p*p));
  }else{
    return sqrt(Subs::sqr(tzerr) + Subs::sqr(p*perr));
  }
}

void Subs::Ephem::set(double T0, double period, TSCALE tscale){
  etype  = LINEAR;
  tzero = T0;
  per   = period;
  this->tscale = tscale;
}

void Subs::Ephem::set(double T0, double period, double pdot, TSCALE tscale){
  etype  = QUADRATIC;
  tzero = T0;
  per   = period;
  quad  = pdot;
  this->tscale = tscale;
}

std::ostream& Subs::operator<<(std::ostream& ost, const Ephem& ephem){

  switch(ephem.tscale) {
  case Ephem::HJD:
    ost << "HJD";
    break;
  case Ephem::HMJD:
    ost << "HMJD";
    break;
  case Ephem::BJD:
    ost << "BJD";
    break;
  case Ephem::BMJD:
    ost << "BMJD";
    break;
  default:
    throw Subs_Error("Subs::operator<<(std::ostream&, const Ephem&): unrecognised value for timescale");
  }

  if(ephem.etype == Ephem::LINEAR)
    ost << " linear ";
  else
    ost << " quadratic ";
  
  ost << std::setprecision(15) << ephem.tzero << " " << ephem.tzerr << " " << ephem.per   << " " << ephem.perr;

  if(ephem.etype == Ephem::QUADRATIC)
    ost << " " << std::setprecision(15) << ephem.quad << " " << ephem.qerr;

  return ost;
}

std::istream& Subs::operator>>(std::istream& ist, Ephem& ephem){

  if(!ist) return ist;

  std::string stscale;
  ist >> stscale;
  if(!ist) return ist;

  stscale = Subs::toupper(stscale);
  if(stscale == "HJD"){
    ephem.tscale = Ephem::HJD;
  }else if(stscale == "HMJD"){
    ephem.tscale = Ephem::HMJD;
  }else if(stscale == "BJD"){
    ephem.tscale = Ephem::BJD;
  }else if(stscale == "BMJD"){
    ephem.tscale = Ephem::BMJD;
  }else{
    ist.setstate(std::ios_base::failbit);
    return ist;
  }

  std::string stype;
  ist >> stype;
  if(!ist) return ist;
  stype = Subs::toupper(stype);
  if(stype == "LINEAR"){
    ephem.etype = Ephem::LINEAR;
  }else{
    ephem.etype = Ephem::QUADRATIC;
  }
  
  ist >> ephem.tzero >> ephem.tzerr >> ephem.per >> ephem.perr;
  if(!ist) return ist;

  if(ephem.etype == Ephem::QUADRATIC) ist >> ephem.quad >> ephem.qerr;

  return ist;

}
