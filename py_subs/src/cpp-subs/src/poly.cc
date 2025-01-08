#include "trm/poly.h"

/** Constructor of a poly representing a linear scale running from
 * 1 to npix
 */
Subs::Poly::Poly(int npix) : Array1D<double>(2), norm(true), middle((npix-1)/2.), hrange((npix-1)/2.) {
  buff[0] = (npix+1)/2.;
  buff[1] = (npix-1)/2.;
}

/** Constructor of a linear polynomial running from a value of xs at pixel -0.5 (left edge of first one)
 * to xe at pixel (npix-0.5). The xs and xe values are also used to define the scaling.
 */
Subs::Poly::Poly(double xs, double xe, int npix) : Array1D<double>(2), norm(true), middle((xs+xe)/2.), hrange(abs(xe-xs)/2.) {
  buff[1] = hrange*(xe-xs)/npix;
  buff[0] = xs - buff[1]*(-0.5-middle)/hrange;
}

//! Another general constructor
Subs::Poly::Poly(bool normal, double xs, double xe, const Array1D<double>& coeff) : 
  Array1D<double>(coeff), norm(normal), middle((xs+xe)/2.), hrange(abs(xe-xs)/2.) {}

//! Another general constructor
Subs::Poly::Poly(bool normal, double xs, double xe, const std::vector<double>& coeff) : 
  Array1D<double>(coeff), norm(normal), middle((xs+xe)/2.), hrange(abs(xe-xs)/2.) {}


Subs::Poly::Poly(const Poly& obj) : Array1D<double>(obj), norm(obj.norm), middle(obj.middle), hrange(obj.hrange) {}

Subs::Poly::Poly(int npoly, bool normal, double xs, double xe) : 
  Array1D<double>(npoly), norm(normal), middle((xs+xe)/2.), hrange(abs(xe-xs)/2.) {}

double Subs::Poly::get_value(double x) const {
  if(this->size() == 0)
    throw Subs_Error("Subs::Poly::Get_value(double): undefined operation on null poly");
  double sum   = buff[0];
  if(this->size() > 1){
    // Nasty little bit designed to save a multiplication.
    double fac   = (x-middle)/hrange;
    double power = fac;
    for(int i=1; i<this->size()-1; i++){
      sum   += buff[i]*power;
      power *= fac;
    }
    sum   += buff[this->size()-1]*power;
  }
  if(norm)
    return sum;
  else
    return exp(sum);
}

double Subs::Poly::get_deriv(double x) const {
  double fac = (x-middle)/hrange;
  double power = fac;
  double sum   = buff[1];
  int nminus = size()-1;
  for(int i=2; i<nminus; i++){
    sum   += i*buff[i]*power;
    power *= fac;
  }
  sum += nminus*buff[nminus]*power;
  if(norm)
    return sum/hrange;
  else
    return sum*get_value(x)/hrange;
}

/** This routine tries to find the value of X equivalent to a given value
 * of the polynomial. The polynomial is assumed to be monotonic. It carries our
 * a Newton-Raphson starting from xinit, so the poly should be reasonably
 * well-behaved and xguess should be good or else there will be trouble.
 * \param value the value to find the equivalent X for
 * \param xguess initial value of X
 * \param acc the precision to find X to.
 */
double Subs::Poly::get_x(double value, double xguess, double acc) const {
  
  const int ITMAX = 500;
  
  double xold = xguess + 2.*acc;
  int nit = 0;
  while(abs(xold-xguess) > acc && nit < ITMAX){
    xold = xguess;
    nit++;
    xguess -= (get_value(xguess)-value)/get_deriv(xguess);
  }
  if(nit == ITMAX)
    std::cerr << "WARNING: hit maximum iterations in Subs::Poly::get_x(double, double, double, double, double, double)" << std::endl;
  
  return xguess;
}


void Subs::Poly::write(std::ofstream& s) const {
  s.write((char*)&norm,   sizeof(bool));
  s.write((char*)&middle, sizeof(double));
  s.write((char*)&hrange, sizeof(double));
  Array1D<double>::write(s);
}

void Subs::Poly::skip(std::ifstream& s, bool swap_bytes) {
  s.ignore(sizeof(bool));
  s.ignore(sizeof(double));
  s.ignore(sizeof(double));
  Array1D<double>::skip(s, swap_bytes);
}

void Subs::Poly::read(std::ifstream& s, bool swap_bytes){
  s.read((char*)&norm,   sizeof(bool));
  s.read((char*)&middle, sizeof(double));
  s.read((char*)&hrange, sizeof(double));
  Array1D<double>::read(s, swap_bytes);
}

std::ostream& Subs::operator<<(std::ostream& s, const Subs::Poly& poly){
  if(!s) return s;
  s << poly.norm << " " << poly.middle << " " << poly.hrange << " ";
  poly.ascii_output(s);
  return s;
}

std::istream& Subs::operator>>(std::istream& s, Subs::Poly& poly){
  if(!s) return s;
  s >> poly.norm >> poly.middle >>  poly.hrange;
  poly.ascii_input(s);
  return s;
}
