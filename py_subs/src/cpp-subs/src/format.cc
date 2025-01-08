#include "trm/format.h"


Subs::Format::Format(int precision) : 
  precision_(precision), width_(0), format_(std::ios::fmtflags(0)), upper(false), 
  showpnt(true), fill_char(' '), fadjust(std::ios::left) {}


void Subs::Format::precision(int p) {
  precision_ = p;
}

void Subs::Format::scientific() {
  format_ = std::ios::scientific;
}

void Subs::Format::fixed(){
  format_ = std::ios::fixed;
}

void Subs::Format::general(){
  format_ = std::ios::fmtflags(0);
}

void Subs::Format::width(int w){
  width_  = w;
}

void Subs::Format::showpoint(){
  showpnt = true;
}

void Subs::Format::noshowpoint(){
  showpnt = false;
}

void Subs::Format::uppercase() {
  upper = true;
}

void Subs::Format::lowercase(){
  upper = false;
}

void Subs::Format::fill(char c){
  fill_char = c;
}

void Subs::Format::left(){
  fadjust = std::ios::left;
}


void Subs::Format::right(){
  fadjust = std::ios::right;
}

void Subs::Format::internal(){
  fadjust = std::ios::internal;
}

//! Extractor for doubles
std::ostream& Subs::operator<<(std::ostream& ostr, const Subs::Bound_form_d& bf){
  std::ostringstream  s;
  s.precision(bf.form.precision_);
  s.width(bf.form.width_);
  s.fill(bf.form.fill_char);
  s.setf(bf.form.format_, std::ios::floatfield);
  s.setf(bf.form.fadjust, std::ios::adjustfield);
  if(bf.form.upper) s.setf(std::ios::uppercase);
  s << bf.val;
  return ostr << s.str();
}

//! Extractor for floats
std::ostream& Subs::operator<<(std::ostream& ostr, const Subs::Bound_form_f& bf){
  std::ostringstream  s;
  s.precision(bf.form.precision_);
  s.width(bf.form.width_);
  s.fill(bf.form.fill_char);
  s.setf(bf.form.format_, std::ios::floatfield);
  s.setf(bf.form.fadjust, std::ios::adjustfield);
  if(bf.form.upper) s.setf(std::ios::uppercase);
  if(bf.form.showpnt) s.setf(std::ios::showpoint);
  s << bf.val;
  return ostr << s.str();
}

//! Extractor for strings
std::ostream& Subs::operator<<(std::ostream& ostr, const Subs::Bound_form_s& bf){
  std::ostringstream s;
  s.width(bf.form.width_);
  s.fill(bf.form.fill_char);
  s.setf(bf.form.fadjust, std::ios::adjustfield);
  if(bf.form.upper) s.setf(std::ios::uppercase);
  s << bf.val;
  return ostr << s.str();
}

//! Returns combination of format and double object
Subs::Bound_form_d Subs::Format::operator()(double d) const {
  return Bound_form_d(*this, d);
}

//! Returns combination of format and float object
Subs::Bound_form_f Subs::Format::operator()(float f) const {
  return Bound_form_f(*this, f);
}

//! Returns combination of format and string object
Subs::Bound_form_s Subs::Format::operator()(const std::string& s) const {
    return Bound_form_s(*this, s);
}

//! Returns combination of format and string object
Subs::Bound_form_s Subs::Format::operator()(const std::string& s, int width) {
    this->width(width);
    return Bound_form_s(*this, s);
}



















