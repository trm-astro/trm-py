//
// header item class; see .h file for more details.
//

#include <iomanip>
#include "trm/hitem.h"

int    Subs::Hitem::pmode            = 3;
size_t Subs::Hitem::value_width      = 30; 
size_t Subs::Hitem::type_width       = 12; 
size_t Subs::Hitem::double_precision = 17; 
size_t Subs::Hitem::float_precision  = 9; 

void Subs::Hitem::set_default_output(){
    pmode            = 1;
    value_width      = 20; 
    double_precision = 17; 
    float_precision  = 9;
}

Subs::Hitem::HITEM_TYPE Subs::Hitem::read_type(std::ifstream& s, bool swap_bytes) {
    int itype;
    s.read((char*)&itype,sizeof(int));
    if(!s) throw Hitem_Error("Subs::Hitem::read_type(std::ifstream&): failed to read type");
    if(swap_bytes) itype = Subs::byte_swap(itype);
    return Hitem::HITEM_TYPE(itype);
}


void Subs::Hitem::write_type(std::ofstream& s) const {
    int type = this->type_code();
    s.write((char*)&type,sizeof(int));
    if(!s) throw Hitem_Error("Subs::Hitem::write_type(std::ofstream&): failed to write type");
}

// Default versions of virtual functions

void Subs::Hitem::get_value(char &cval) const {
    throw Hitem_Error("Cannot retrieve a char from " + type() + "-type header items");
}

void Subs::Hitem::get_value(unsigned char &cval) const {
    throw Hitem_Error("Cannot retrieve an unsigned char from " + type() + "-type header items");
}

void Subs::Hitem::get_value(double &dval) const {
    throw Hitem_Error("Cannot retrieve a double from " + type() + "-type header items");
}

void Subs::Hitem::get_value(int &ival) const {
    throw Hitem_Error("Cannot retrieve an int from " + type() + "-type header items");
}

void Subs::Hitem::get_value(unsigned int &ival) const {
    throw Hitem_Error("Cannot retrieve an unsigned int from " + type() + "-type header items");
}

void Subs::Hitem::get_value(long int &ival) const {
    throw Hitem_Error("Cannot retrieve a long int from " + type() + "-type header items");
}

void Subs::Hitem::get_value(unsigned long int &ival) const {
    throw Hitem_Error("Cannot retrieve an unsigned long int from " + type() + "-type header items");
}

void Subs::Hitem::get_value(float &fval) const {
    throw Hitem_Error("Cannot retrieve a float from " + type() + "-type header items");
}

void Subs::Hitem::get_value(bool &bval) const {
    throw Hitem_Error("Cannot retrieve a bool from " + type() + "-type header items");
}

void Subs::Hitem::get_value(Date &dval) const {
    throw Hitem_Error("Cannot retrieve a Date from " + type() + " header items");
}

void Subs::Hitem::get_value(Time &tval) const {
    throw Hitem_Error("Cannot retrieve a Time from " + type() + "-type header items");
}

void Subs::Hitem::get_value(Position &pval) const {
    throw Hitem_Error("Cannot retrieve a Position from " + type() + "-type header items");
}

void Subs::Hitem::get_value(Telescope &tval) const {
    throw Hitem_Error("Cannot retrieve a Telescope from " + type() + "-type header items");
}

void Subs::Hitem::get_value(std::vector<double>& val) const {
    throw Hitem_Error("Cannot retrieve a std::vector<double> from " + type() + "-type header items");
}

char Subs::Hitem::get_char() const {
    throw Hitem_Error("Cannot retrieve a char from " + type() + "-type header items");
}

unsigned char Subs::Hitem::get_uchar() const {
    throw Hitem_Error("Cannot retrieve an unsigned char from " + type() + "-type header items");
}

double Subs::Hitem::get_double() const {
    throw Hitem_Error("Cannot retrieve a double from " + type() + "-type header items");
}

int Subs::Hitem::get_int() const {
    throw Hitem_Error("Cannot retrieve an int from " + type() + "-type header items");
}

unsigned int Subs::Hitem::get_uint() const {
    throw Hitem_Error("Cannot retrieve an unsigned int from " + type() + "-type header items");
}

long int Subs::Hitem::get_lint() const {
    throw Hitem_Error("Cannot retrieve a long int from " + type() + "-type header items");
}

unsigned long int Subs::Hitem::get_ulint() const {
    throw Hitem_Error("Cannot retrieve an unsigned long int from " + type() + "-type header items");
}

float Subs::Hitem::get_float() const {
    throw Hitem_Error("Cannot retrieve a float from " + type() + "-type header items");
}

bool Subs::Hitem::get_bool() const {
    throw Hitem_Error("Cannot retrieve a bool from " + type() + "-type header items");
}

const Subs::Date& Subs::Hitem::get_date() const {
    throw Hitem_Error("Cannot retrieve a Date from " + type() + "-type header items");
}

const Subs::Time& Subs::Hitem::get_time() const {
    throw Hitem_Error("Cannot retrieve a Time from " + type() + "-type header items");
}

const Subs::Position& Subs::Hitem::get_position() const {
    throw Hitem_Error("Cannot retrieve a Position from " + type() + "-type header items");
}

const Subs::Telescope& Subs::Hitem::get_telescope() const {
    throw Hitem_Error("Cannot retrieve a Telescope from " + type() + "-type header items");
}

const std::vector<double>& Subs::Hitem::get_dvector() const {
    throw Hitem_Error("Cannot retrieve a std::vector<double> from " + type() + "-type header items");
}

const std::vector<int>& Subs::Hitem::get_ivector() const {
    throw Hitem_Error("Cannot retrieve a std::vector<int> from " + type() + "-type header items");
}

const std::vector<float>& Subs::Hitem::get_fvector() const {
    throw Hitem_Error("Cannot retrieve a std::vector<float> from " + type() + "-type header items");
}

void Subs::Hitem::set_value(const char &cval){
    throw Hitem_Error("Cannot store a char in " + type() +  "-type header items.");
}

void Subs::Hitem::set_value(const unsigned char &cval){
    throw Hitem_Error("Cannot store an unsigned char in " + type() + "-type header items.");
}

void Subs::Hitem::set_value(const double &dval){
    throw Hitem_Error("Cannot store a double in " + type() + "-type header items.");
}

void Subs::Hitem::set_value(const int &ival){
    throw Hitem_Error("Cannot store an int in " + type() + "-type header items.");
}

void Subs::Hitem::set_value(const unsigned int &ival){
    throw Hitem_Error("Cannot store an unsigned int in " + type() + "-type header items.");
}

void Subs::Hitem::set_value(const long int &ival){
    throw Hitem_Error("Cannot store a long int in " + type() + "-type header items.");
}

void Subs::Hitem::set_value(const unsigned long int &ival){
    throw Hitem_Error("Cannot store an unsigned long int in " + type() + "-type header items.");
}

void Subs::Hitem::set_value(const float &fval){
    throw Hitem_Error("Cannot store a float in " + type() + "-type header items.");
}

void Subs::Hitem::set_value(const bool &bval){
    throw Hitem_Error("Cannot store a bool in " + type() + "-type header items.");
}

void Subs::Hitem::set_value(const Date &dval){
    throw Hitem_Error("Cannot store a Date in " + type() + "-type header items.");
}

void Subs::Hitem::set_value(const Time &tval){
    throw Hitem_Error("Cannot store a Time in " + type() + "-type header items.");
}

void Subs::Hitem::set_value(const Position& pval){
    throw Hitem_Error("Cannot store a Position in " + type() + "-type header items.");
}

void Subs::Hitem::set_value(const Telescope& tval){
    throw Hitem_Error("Cannot store a Telescope in " + type() + "-type header items.");
}

void Subs::Hitem::set_value(const std::vector<double>& val){
    throw Hitem_Error("Cannot store a std::vector<double> in " + type() + "-type header items.");
}

// Output

std::ostream& Subs::operator<<(std::ostream& s, const Subs::Hitem* obj){
    s.setf(std::ios::left);
    if(Subs::Hitem::pmode > 0) obj->print(s);
    if(Subs::Hitem::pmode > 1){
	s << " /" << obj->type();
	if(obj->type().length() < Hitem::type_width-2){
	    s.setf(std::ios::left, std::ios::adjustfield);
	    s.width(Hitem::type_width-obj->type().length() - 1);
	}
	s << "/";
    }
    if(Subs::Hitem::pmode > 2) s << " " << obj->get_comment();

    return s;
}

// Hchar
 
void Subs::Hchar::print(std::ostream& s) const {
    s.fill(' ');
    s.width(value_width);
    s.setf(std::ios::left, std::ios::adjustfield);
    s << value;
}

void Subs::Hchar::write(std::ofstream& s) const {
    write_string(s,comment);
    s.write((char*)&value,sizeof(char));
    if(!s) throw Hitem_Error("Subs::Hchar::write(std::ofstream&): failed to write item");
}

void Subs::Hchar::write_ascii(std::ofstream& s) const {
    s << "C " << value << " " << comment << std::endl;
    if(!s) throw Hitem_Error("Subs::Hchar::write_ascii(std::ofstream&): failed to write item");
}

void Subs::Hchar::read(std::ifstream& s, bool swap_bytes) {
    read_string(s,comment,swap_bytes);
    s.read((char*)&value,sizeof(char));
    if(!s) throw Hitem_Error("Subs::Hchar::read(std::ifstream&): failed to read item");
}

void Subs::Hchar::skip(std::ifstream& s, bool swap_bytes) {
    skip_string(s,swap_bytes);
    s.ignore(sizeof(char));
    if(!s) throw Hitem_Error("Subs::Hchar::skip(std::ifstream&): failed to skip item");
}

void Subs::Hchar::get_value(std::string &sval) const {
    std::ostringstream ostr;
    ostr << value;
    if(!ostr)
	throw Hitem_Error("void Subs::Hchar::get_value(std::string&) const: error translating a char into a string.");
    sval = ostr.str();
}

std::string Subs::Hchar::get_string() const {
    std::ostringstream ostr;
    ostr << value;
    if(!ostr)
	throw Hitem_Error("void Subs::Hchar::get_string() const: error translating a char into a string.");
    return ostr.str();
}

void Subs::Hchar::set_value(const std::string &sval){
    std::istringstream istr(sval);
    istr >> value;
    if(!istr)
	throw Hitem_Error("void Subs::Hitem::set_value(const std::string&): failed to translate \"" + sval + "\" into a char");
}

// Hdouble
 
void Subs::Hdouble::print(std::ostream& s) const {
    s.fill(' ');
    s.width(value_width);
    s.setf(std::ios::left, std::ios::adjustfield);
    s.precision(double_precision < value_width ? double_precision : value_width);
    s << value;
}

void Subs::Hdouble::write(std::ofstream& s) const {
    write_string(s,comment);
    s.write((char*)&value,sizeof(double));
    if(!s) throw Hitem_Error("Subs::Hdouble::write(std::ofstream&): failed to write item");
}

void Subs::Hdouble::write_ascii(std::ofstream& s) const {
    s << "D " << value << " " << comment << std::endl;
    if(!s) throw Hitem_Error("Subs::Hdouble::write_ascii(std::ofstream&): failed to write item");
}

void Subs::Hdouble::read(std::ifstream& s, bool swap_bytes) {
    read_string(s,comment,swap_bytes);
    s.read((char*)&value,sizeof(double));
    if(!s) throw Hitem_Error("Subs::Hdouble::read(std::ifstream&): failed to read item");
    if(swap_bytes) value = Subs::byte_swap(value);
}

void Subs::Hdouble::skip(std::ifstream& s, bool swap_bytes) {
    skip_string(s,swap_bytes);
    s.ignore(sizeof(double));
    if(!s) throw Hitem_Error("Subs::Hdouble::skip(std::ifstream&): failed to skip item");
}

void Subs::Hdouble::get_value(std::string &sval) const {
    std::ostringstream ostr;
    ostr << std::setprecision(17) << value;
    if(!ostr)
	throw Hitem_Error("void Subs::Hdouble::get_value(std::string&) const: error translating a double into a string.");
    sval = ostr.str();
}

std::string Subs::Hdouble::get_string() const {
  
    std::ostringstream ostr;
    ostr << std::setprecision(17) << value;
    if(!ostr)
	throw Hitem_Error("void Subs::Hdouble::get_string() const: error translating a double into a string.");
    return ostr.str();
}

void Subs::Hdouble::set_value(const std::string &sval){
    std::istringstream istr(sval);
    istr >> value;
    if(!istr)
	throw Hitem_Error("void Subs::Hitem::set_value(const std::string&): failed to translate \"" + sval + "\" into a double");
}

// int

void Subs::Hint::print(std::ostream& s) const {
    s.fill(' ');
    s.width(value_width);
    s.setf(std::ios::left, std::ios::adjustfield);
    s << value;
}

void Subs::Hint::write(std::ofstream& s) const {
    write_string(s,comment);
    s.write((char*)&value,sizeof(INT4));
    if(!s) throw Hitem_Error("Subs::Hint::write(std::ofstream&): failed to write item");
}

void Subs::Hint::write_ascii(std::ofstream& s) const {
    s << "I " << value << " " << comment << std::endl;
    if(!s) throw Hitem_Error("Subs::Hint::write_ascii(std::ofstream&): failed to write item");
}

void Subs::Hint::read(std::ifstream& s, bool swap_bytes) {
    read_string(s,comment,swap_bytes);
    s.read((char*)&value,sizeof(INT4));
    if(!s) throw Hitem_Error("Subs::Hint::read(std::ifstream&): failed to read item");
    if(swap_bytes) value = Subs::byte_swap(value);
}

void Subs::Hint::skip(std::ifstream& s, bool swap_bytes) {
    skip_string(s,swap_bytes);
    s.ignore(sizeof(INT4));
    if(!s) throw Hitem_Error("Subs::Hint::skip(std::ifstream&): failed to skip item");
}

void Subs::Hint::get_value(std::string &sval) const {
    std::ostringstream ostr;
    ostr << value;
    if(!ostr)
	throw Hitem_Error("void Subs::Hint::get_value(std::string&) const: error translating an int into a string.");
    sval = ostr.str();
}

std::string Subs::Hint::get_string() const {
  
    std::ostringstream ostr;
    ostr << value;
    if(!ostr)
	throw Hitem_Error("void Subs::Hint::get_string() const: error translating an int into a string.");
    return ostr.str();
}

void Subs::Hint::set_value(const std::string &sval){
    std::istringstream istr(sval);
    istr >> value;
    if(!istr)
	throw Hitem_Error("void Subs::Hitem::set_value(const std::string&): failed to translate \"" + sval + "\" into an int");
}

// unsigned int

void Subs::Huint::print(std::ostream& s) const {
    s.fill(' ');
    s.width(value_width);
    s.setf(std::ios::left, std::ios::adjustfield);
    s << value;
}

void Subs::Huint::write(std::ofstream& s) const {
    write_string(s,comment);
    s.write((char*)&value,sizeof(UINT4));
    if(!s) throw Hitem_Error("Subs::Huint::write(std::ofstream&): failed to write item");
}

void Subs::Huint::write_ascii(std::ofstream& s) const {
    s << "UI " << value << " " << comment << std::endl;
    if(!s) throw Hitem_Error("Subs::Huint::write_ascii(std::ofstream&): failed to write item");
}

void Subs::Huint::read(std::ifstream& s, bool swap_bytes) {
    read_string(s,comment,swap_bytes);
    s.read((char*)&value,sizeof(UINT4));
    if(!s) throw Hitem_Error("Subs::Huint::read(std::ifstream&): failed to read item");
    if(swap_bytes) value = Subs::byte_swap(value);
}

void Subs::Huint::skip(std::ifstream& s, bool swap_bytes) {
    skip_string(s,swap_bytes);
    s.ignore(sizeof(UINT4));
    if(!s) throw Hitem_Error("Subs::Huint::skip(std::ifstream&): failed to skip item");
}

void Subs::Huint::get_value(std::string &sval) const {
    std::ostringstream ostr;
    ostr << value;
    if(!ostr)
	throw Hitem_Error("void Subs::Huint::get_value(std::string&) const: error translating an unsigned int into a string.");
    sval = ostr.str();
}

std::string Subs::Huint::get_string() const {
  
    std::ostringstream ostr;
    ostr << value;
    if(!ostr)
	throw Hitem_Error("void Subs::Huint::get_string() const: error translating an unsigned int into a string.");
    return ostr.str();
}

void Subs::Huint::set_value(const std::string &sval){
    std::istringstream istr(sval);
    istr >> value;
    if(!istr)
	throw Hitem_Error("void Subs::Hitem::set_value(const std::string&): failed to translate \"" + sval + "\" into an unsigned int");
}

// unsigned short int

void Subs::Husint::print(std::ostream& s) const {
    s.fill(' ');
    s.width(value_width);
    s.setf(std::ios::left, std::ios::adjustfield);
    s << value;
}

void Subs::Husint::write(std::ofstream& s) const {
    write_string(s,comment);
    s.write((char*)&value,sizeof(UINT2));
    if(!s) throw Hitem_Error("Subs::Husint::write(std::ofstream&): failed to write item");
}

void Subs::Husint::write_ascii(std::ofstream& s) const {
    s << "UI " << value << " " << comment << std::endl;
    if(!s) throw Hitem_Error("Subs::Husint::write_ascii(std::ofstream&): failed to write item");
}

void Subs::Husint::read(std::ifstream& s, bool swap_bytes) {
    read_string(s,comment,swap_bytes);
    s.read((char*)&value,sizeof(UINT2));
    if(!s) throw Hitem_Error("Subs::Husint::read(std::ifstream&): failed to read item");
    if(swap_bytes) value = Subs::byte_swap(value);
}

void Subs::Husint::skip(std::ifstream& s, bool swap_bytes) {
    skip_string(s,swap_bytes);
    s.ignore(sizeof(UINT2));
    if(!s) throw Hitem_Error("Subs::Husint::skip(std::ifstream&): failed to skip item");
}

void Subs::Husint::get_value(std::string &sval) const {
    std::ostringstream ostr;
    ostr << value;
    if(!ostr)
	throw Hitem_Error("void Subs::Husint::get_value(std::string&) const: error translating an unsigned short int into a string.");
    sval = ostr.str();
}

std::string Subs::Husint::get_string() const {
  
    std::ostringstream ostr;
    ostr << value;
    if(!ostr)
	throw Hitem_Error("void Subs::Husint::get_string() const: error translating an unsigned short int into a string.");
    return ostr.str();
}

void Subs::Husint::set_value(const std::string &sval){
    std::istringstream istr(sval);
    istr >> value;
    if(!istr)
	throw Hitem_Error("void Subs::Hitem::set_value(const std::string&): failed to translate \"" + sval + "\" into an unsigned short int");
}

// long int

void Subs::Hlint::print(std::ostream& s) const {
    s.fill(' ');
    s.width(value_width);
    s.setf(std::ios::left, std::ios::adjustfield);
    s << value;
}

void Subs::Hlint::write(std::ofstream& s) const {
    write_string(s,comment);
    s.write((char*)&value,sizeof(INT4));
    if(!s) throw Hitem_Error("Subs::Hlint::write(std::ofstream&): failed to write item");
}

void Subs::Hlint::write_ascii(std::ofstream& s) const {
    s << "LI " << value << " " << comment << std::endl;
    if(!s) throw Hitem_Error("Subs::Hlint::write_ascii(std::ofstream&): failed to write item");
}

void Subs::Hlint::read(std::ifstream& s, bool swap_bytes) {
    read_string(s,comment,swap_bytes);
    s.read((char*)&value,sizeof(INT4));
    if(!s) throw Hitem_Error("Subs::Hlint::read(std::ifstream&, bool): failed to read item");
    if(swap_bytes) value = Subs::byte_swap(value);
}

void Subs::Hlint::skip(std::ifstream& s, bool swap_bytes) {
    skip_string(s,swap_bytes);
    s.ignore(sizeof(INT4));
    if(!s) throw Hitem_Error("Subs::Hlint::skip(std::ifstream&, bool): failed to skip item");
}

void Subs::Hlint::get_value(std::string &sval) const {
    std::ostringstream ostr;
    ostr << value;
    if(!ostr)
	throw Hitem_Error("void Subs::Hlint::get_value(std::string&) const: error translating a long int into a string.");
    sval = ostr.str();
}

std::string Subs::Hlint::get_string() const {
  
    std::ostringstream ostr;
    ostr << value;
    if(!ostr)
	throw Hitem_Error("void Subs::Hlint::get_string() const: error translating a long int into a string.");
    return ostr.str();
}

void Subs::Hlint::set_value(const std::string &sval){
    std::istringstream istr(sval);
    istr >> value;
    if(!istr)
	throw Hitem_Error(std::string("void Subs::Hitem::set_value(const std::string&): failed to translate \"") + sval + 
			  std::string("\" into a long int"));
}

// unsigned long int

void Subs::Hulint::print(std::ostream& s) const {
    s.fill(' ');
    s.width(value_width);
    s.setf(std::ios::left, std::ios::adjustfield);
    s << value;
}

void Subs::Hulint::write(std::ofstream& s) const {
    write_string(s,comment);
    s.write((char*)&value,sizeof(UINT4));
    if(!s) throw Hitem_Error("Subs::Hulint::write(std::ofstream&): failed to write item");
}

void Subs::Hulint::write_ascii(std::ofstream& s) const {
    s << "ULI " << value << " " << comment << std::endl;
    if(!s) throw Hitem_Error("Subs::Hulint::write_ascii(std::ofstream&): failed to write item");
}

void Subs::Hulint::read(std::ifstream& s, bool swap_bytes) {
    read_string(s,comment,swap_bytes);
    s.read((char*)&value,sizeof(UINT4));
    if(!s) throw Hitem_Error("Subs::Hulint::read(std::ifstream&, bool): failed to read item");
    if(swap_bytes) value = Subs::byte_swap(value);
}

void Subs::Hulint::skip(std::ifstream& s, bool swap_bytes) {
    skip_string(s,swap_bytes);
    s.ignore(sizeof(UINT4));
    if(!s) throw Hitem_Error("Subs::Hulint::skip(std::ifstream&, bool): failed to skip item");
}

void Subs::Hulint::get_value(std::string &sval) const {
    std::ostringstream ostr;
    ostr << value;
    if(!ostr)
	throw Hitem_Error("void Subs::Hulint::get_value(std::string&) const: error translating an unsigned long int into a string.");
    sval = ostr.str();
}

std::string Subs::Hulint::get_string() const {
    std::ostringstream ostr;
    ostr << value;
    if(!ostr)
	throw Hitem_Error("void Subs::Hulint::get_string() const: error translating an unsigned long int into a string.");
    return ostr.str();
}

void Subs::Hulint::set_value(const std::string &sval){
    std::istringstream istr(sval);
    istr >> value;
    if(!istr)
	throw Hitem_Error(std::string("void Subs::Hitem::set_value(const std::string&): failed to translate \"") + sval + 
			  std::string("\" into an unsigned long int"));
}

// float
void Subs::Hfloat::print(std::ostream& s) const {
    s.fill(' ');
    s.width(value_width);
    s.setf(std::ios::left, std::ios::adjustfield);
    s.precision(float_precision < value_width ? float_precision : value_width);
    s << value;
}

void Subs::Hfloat::write(std::ofstream& s) const {
    write_string(s,comment);
    s.write((char*)&value,sizeof(REAL4));
    if(!s) throw Hitem_Error("Subs::Hfloat::write(std::ofstream&): failed to write item");
}

void Subs::Hfloat::write_ascii(std::ofstream& s) const {
    s << "F " << value << " " << comment << std::endl;
    if(!s) throw Hitem_Error("Subs::Hfloat::write_ascii(std::ofstream&): failed to write item");
}

void Subs::Hfloat::read(std::ifstream& s, bool swap_bytes) {
    read_string(s,comment,swap_bytes);
    s.read((char*)&value,sizeof(REAL4));
    if(!s) throw Hitem_Error("Subs::Hfloat::read(std::ifstream&): failed to read item");
    if(swap_bytes) value = Subs::byte_swap(value);
}

void Subs::Hfloat::skip(std::ifstream& s, bool swap_bytes) {
    skip_string(s, swap_bytes);
    s.ignore(sizeof(REAL4));
    if(!s) throw Hitem_Error("Subs::Hfloat::skip(std::ifstream&): failed to skip item");
}

void Subs::Hfloat::get_value(std::string &sval) const {
    std::ostringstream ostr;
    ostr << std::setprecision(9) << value;
    if(!ostr)
	throw Hitem_Error("void Subs::Hfloat::get_value(std::string&) const: error translating a float into a string.");
    sval = ostr.str();
}

std::string Subs::Hfloat::get_string() const {
    std::ostringstream ostr;
    ostr << std::setprecision(9) << value;
    if(!ostr)
	throw Hitem_Error("void Subs::Hfloat::get_string() const: error translating a float into a string.");
    return ostr.str();
}

void Subs::Hfloat::set_value(const std::string &sval){
    std::istringstream istr(sval);
    istr >> value;
    if(!istr)
	throw Hitem_Error(std::string("void Subs::Hitem::set_value(const std::string&): failed to translate \"") + sval + 
			  std::string("\" into a float"));
}

// boolean
void Subs::Hbool::print(std::ostream& s) const {
    s.fill(' ');
    s.width(value_width);
    s.setf(std::ios::left, std::ios::adjustfield);
    s << value;
}

// 'bool' variables not well defined in terms of bytes so write/read as UCHAR
void Subs::Hbool::write(std::ofstream& s) const {
    write_string(s,comment);
    UCHAR c = UCHAR(value);
    s.write((char*)&c,sizeof(UCHAR));
    if(!s) throw Hitem_Error("Subs::Hbool::write(std::ofstream&): failed to write item");
}

void Subs::Hbool::write_ascii(std::ofstream& s) const {
    s << "B ";
    if(value)
	s << "TRUE" << " " << comment << std::endl;
    else
	s << "TRUE" << " " << comment << std::endl;
    if(!s) throw Hitem_Error("Subs::Hbool::write_ascii(std::ofstream&): failed to write item");
}

void Subs::Hbool::read(std::ifstream& s, bool swap_bytes) {
    read_string(s,comment,swap_bytes);
    UCHAR c;
    s.read((char*)&c,sizeof(UCHAR));
    if(!s) throw Hitem_Error("Subs::Hbool::read(std::ifstream&): failed to read item");
    value = bool(c);
}

void Subs::Hbool::skip(std::ifstream& s, bool swap_bytes) {
    skip_string(s, swap_bytes);
    s.ignore(sizeof(UCHAR));
    if(!s) throw Hitem_Error("Subs::Hbool::skip(std::ifstream&): failed to skip item");
}

void Subs::Hbool::get_value(std::string &sval) const {
    std::ostringstream ostr;
    ostr << value;
    if(!ostr)
	throw Hitem_Error("void Subs::Hbool::get_value(std::string&) const: error translating a bool into a string.");
    sval = ostr.str();
}

std::string Subs::Hbool::get_string() const {
    std::ostringstream ostr;
    ostr << value;
    if(!ostr)
	throw Hitem_Error("void Subs::Hbool::get_string() const: error translating a bool into a string.");
    return ostr.str();
}

void Subs::Hbool::set_value(const std::string &sval){
    std::istringstream istr(sval);
    istr >> value;
    if(!istr)
	throw Hitem_Error(std::string("void Subs::Hitem::set_value(const std::string&): failed to translate \"") + sval + 
			  std::string("\" into a bool"));
}

// string
void Subs::Hstring::print(std::ostream& s) const {
    s.fill(' ');
    s.setf(std::ios::left, std::ios::adjustfield);
    std::string::size_type smax = std::min(value.length()-1, value.find_last_not_of(" \t"))+1;
    s << value.substr(0,smax);
    if(smax < value_width){
	s.width(value_width-smax);
	s << " ";
    }
}

void Subs::Hstring::write(std::ofstream& s) const {
    write_string(s,comment);
    write_string(s,value);
}

void Subs::Hstring::write_ascii(std::ofstream& s) const {
    s << "S " << value << " # " << comment << std::endl;
    if(!s) throw Hitem_Error("Subs::Hstring::write_ascii(std::ofstream&): failed to write item");
}

void Subs::Hstring::read(std::ifstream& s, bool swap_bytes) {
    read_string(s,comment,swap_bytes);
    read_string(s,value,swap_bytes);
}

void Subs::Hstring::skip(std::ifstream& s, bool swap_bytes) {
    skip_string(s, swap_bytes);
    skip_string(s, swap_bytes);
}

// directory
void  Subs::Hdirectory::get_value(std::string &sval) const {
    throw Hitem_Error("void  Subs::Hdirectory::get_value(std::string& ) const: error, this should not have been called");
}

std::string Subs::Hdirectory::get_string() const {
    throw Hitem_Error("void  Subs::Hdirectory::get_string() const: error, this should not have been called");
}

void  Subs::Hdirectory::set_value(const std::string &sval){
    throw Hitem_Error("void  Subs::Hdirectory::set_value(const std::string&): error, this should not have been called");
}

void Subs::Hdirectory::print(std::ostream& s) const {}

void Subs::Hdirectory::write(std::ofstream& s) const {
    write_string(s,comment);
}

void Subs::Hdirectory::write_ascii(std::ofstream& s) const {
    s << "DIR " << comment << std::endl;
    if(!s) throw Hitem_Error("Subs::Hdirectory::write_ascii(std::ofstream&): failed to write item");
}

void Subs::Hdirectory::read(std::ifstream& s, bool swap_bytes) {
    read_string(s,comment,swap_bytes);
}

void Subs::Hdirectory::skip(std::ifstream& s, bool swap_bytes) {
    skip_string(s,swap_bytes);
}

// Dates
void Subs::Hdate::print(std::ostream& s) const {
    s << std::setfill(' ') << value;
}

void Subs::Hdate::write(std::ofstream& s) const {
    write_string(s,comment);
    value.write(s);
}

void Subs::Hdate::write_ascii(std::ofstream& s) const {
    s << "DATE " << value << " " << comment << std::endl;
    if(!s) throw Hitem_Error("Subs::Hdate::write_ascii(std::ofstream&): failed to write item");
}

void Subs::Hdate::read(std::ifstream& s, bool swap_bytes) {
    read_string(s,comment,swap_bytes);
    value.read(s,swap_bytes);
}

void Subs::Hdate::skip(std::ifstream& s, bool swap_bytes) {
    skip_string(s,swap_bytes);
    Date::skip(s);
}

void Subs::Hdate::get_value(std::string &sval) const {
    std::ostringstream ostr;
    ostr << value;
    if(!ostr)
	throw Hitem_Error("void Subs::Hdate::get_value(std::string&) const: error translating a Date into a string.");
    sval = ostr.str();
}

std::string Subs::Hdate::get_string() const {
  
    std::ostringstream ostr;
    ostr << value;
    if(!ostr)
	throw Hitem_Error("void Subs::Hdate::get_string() const: error translating a Date into a string.");
    return ostr.str();
}

void Subs::Hdate::set_value(const std::string &sval){
    std::istringstream istr(sval);
    istr >> value;
    if(!istr)
	throw Hitem_Error(std::string("void Subs::Hitem::set_value(const std::string&): failed to translate \"") + sval + 
			  std::string("\" into a Date"));
}

// Times
void Subs::Htime::print(std::ostream& s) const {
    s << std::setfill(' ') << value << "   ";
}

void Subs::Htime::write(std::ofstream& s) const {
    write_string(s,comment);
    value.write(s);
}

void Subs::Htime::write_ascii(std::ofstream& s) const {
    s << "T " << value << " " << comment << std::endl;
    if(!s) throw Hitem_Error("Subs::Htime::write_ascii(std::ofstream&): failed to write item");
}

void Subs::Htime::read(std::ifstream& s, bool swap_bytes) {
    read_string(s,comment,swap_bytes);
    value.read(s,swap_bytes);
}

void Subs::Htime::skip(std::ifstream& s, bool swap_bytes) {
    skip_string(s,swap_bytes);
    Time::skip(s);
}

void Subs::Htime::get_value(std::string &sval) const {  
    std::ostringstream ostr;
    ostr << value;
    if(!ostr)
	throw Hitem_Error("void Subs::Htime::get_value(std::string&) const: error translating a Time into a string.");
    sval = ostr.str();
}

std::string Subs::Htime::get_string() const {
    std::ostringstream ostr;
    ostr << value;
    if(!ostr)
	throw Hitem_Error("void Subs::Htime::get_string() const: error translating a Time into a string.");
    return ostr.str();
}

void Subs::Htime::set_value(const std::string &sval){
    std::istringstream istr(sval);
    istr >> value;
    if(!istr)
	throw Hitem_Error(std::string("void Subs::Hitem::set_value(const std::string&): failed to translate \"") + sval + 
			  std::string("\" into a Time"));
}

// Positions
void Subs::Hposition::print(std::ostream& s) const {
    s << std::setfill(' ') << value;
}

void Subs::Hposition::write(std::ofstream& s) const {
    write_string(s,comment);
    value.write(s);
}

void Subs::Hposition::write_ascii(std::ofstream& s) const {
    s << "POS " << value << " " << comment << std::endl;
    if(!s) throw Hitem_Error("Subs::Hposition::write_ascii(std::ofstream&): failed to write item");
}

void Subs::Hposition::read(std::ifstream& s, bool swap_bytes) {
    read_string(s,comment,swap_bytes);
    value.read(s,swap_bytes);
}

void Subs::Hposition::skip(std::ifstream& s, bool swap_bytes) {
    skip_string(s,swap_bytes);
    Position::skip(s);
}

void Subs::Hposition::get_value(std::string &sval) const {
    std::ostringstream ostr;
    ostr << value;
    if(!ostr)
	throw Hitem_Error("void Subs::Hposition::get_value(std::string&) const: error translating a Position into a string.");
    sval = ostr.str();
}

std::string Subs::Hposition::get_string() const {
    std::ostringstream ostr;
    ostr << value;
    if(!ostr)
	throw Hitem_Error("void Subs::Hposition::get_string() const: error translating a Position into a string.");
    return ostr.str();
}

void Subs::Hposition::set_value(const std::string &sval){
    std::istringstream istr(sval);
    istr >> value;
    if(!istr)
	throw Hitem_Error(std::string("void Subs::Hitem::set_value(const std::string&): failed to translate \"") + sval + 
			  std::string("\" into a Position"));
}

// Telescopes
void Subs::Htelescope::print(std::ostream& s) const {
    s << std::setfill(' ') << value;
}

void Subs::Htelescope::write(std::ofstream& s) const {
    write_string(s,comment);
    value.write(s);
}

void Subs::Htelescope::write_ascii(std::ofstream& s) const {
    s << "TEL " << value << " " << comment << std::endl;
    if(!s) throw Hitem_Error("Subs::Htelescope::write_ascii(std::ofstream&): failed to write item");
}

void Subs::Htelescope::read(std::ifstream& s, bool swap_bytes) {
    read_string(s,comment,swap_bytes);
    value.read(s,swap_bytes);
}

void Subs::Htelescope::skip(std::ifstream& s, bool swap_bytes) {
    skip_string(s,swap_bytes);
    Telescope::skip(s, swap_bytes);
}

void Subs::Htelescope::get_value(std::string &sval) const {
    std::ostringstream ostr;
    ostr << value;
    if(!ostr)
	throw Hitem_Error("void Subs::Htelescope::get_value(std::string&) const: error translating a Telescope into a string.");
    sval = ostr.str();
}

std::string Subs::Htelescope::get_string() const {
    std::ostringstream ostr;
    ostr << value;
    if(!ostr)
	throw Hitem_Error("void Subs::Htelescope::get_string() const: error translating a Telescope into a string.");
    return ostr.str();
}

void Subs::Htelescope::set_value(const std::string &sval){
    std::istringstream istr(sval);
    istr >> value;
    if(!istr)
	throw Hitem_Error(std::string("void Subs::Hitem::set_value(const std::string&): failed to translate \"") + sval + 
			  std::string("\" into a Telescope"));
}

// vector<double>
void Subs::Hdvector::print(std::ostream& s) const {
    s << std::setfill(' ');
    for(size_t i=0; i<value.size(); i++)
	s << value[i] << " ";
}

void Subs::Hdvector::write(std::ofstream& s) const {
    write_string(s,comment);
    UINT4 n = UINT4(value.size());
    s.write((char*)&n, sizeof(UINT4));
    REAL8 d;
    for(size_t i=0; i<value.size(); i++){
	d = value[i];
	s.write((char*)&d, sizeof(REAL8));
    }
    if(!s) throw Hitem_Error("Subs::Hdvector::write(std::ofstream&): failed to write item");
}
 
void Subs::Hdvector::write_ascii(std::ofstream& s) const {
    s << "DV " << value.size();
    for(size_t i=0; i<value.size(); i++){
	s << " " << value[i];
    }
    s << " " << comment << std::endl;
    if(!s) throw Hitem_Error("Subs::Hdvector::write_ascii(std::ofstream&): failed to write item");
}

void Subs::Hdvector::read(std::ifstream& s, bool swap_bytes) {
    read_string(s,comment,swap_bytes);
    UINT4 n;
    s.read((char*)&n,sizeof(UINT4));
    if(swap_bytes) n = Subs::byte_swap(n);
    value.resize(n);
    REAL8 d;
    for(size_t i=0; i<value.size(); i++){
	s.read((char*)&d, sizeof(REAL8));
	if(swap_bytes) d = Subs::byte_swap(d);
	value[i] = d;
    }
    if(!s) throw Hitem_Error("Subs::Hdvector::read(std::ifstream&, bool): failed to read item");
}

void Subs::Hdvector::skip(std::ifstream& s, bool swap_bytes) {
    skip_string(s, swap_bytes);
    UINT4 n;
    s.read((char*)&n,sizeof(UINT4));
    if(swap_bytes) n = Subs::byte_swap(n);
    s.ignore(n*sizeof(REAL8));
    if(!s) throw Hitem_Error("Subs::Hdvector::skip(std::ifstream&, bool): failed to skip item");
}

void  Subs::Hdvector::get_value(std::string &sval) const {
    std::ostringstream ostr;
    for(std::vector<double>::size_type i=0; i<value.size(); i++){
	ostr << value[i] << " ";
	if(!ostr)
	    throw Hitem_Error("void Subs::Hdvector::get_value(std::string&) const: error translating a dvector into a string.");
    }
    sval = ostr.str();
}

std::string Subs::Hdvector::get_string() const {
    std::ostringstream ostr;
    for(std::vector<double>::size_type i=0; i<value.size(); i++){
	ostr << value[i] << " ";
	if(!ostr)
	    throw Hitem_Error("void Subs::Hdvector::get_value(std::string&) const: error translating a dvector into a string.");
    }
    return ostr.str();
}

void Subs::Hdvector::set_value(const std::string &sval){
    // Blank added to avoid apparent bug in istringstream
    std::istringstream istr(sval + " ");
    std::vector<double> temp;
    double d;
    while(!istr.eof()){
	istr >> d;
	if(istr.good()){
	    temp.push_back(d);
	}else if(!istr.eof()){
	    throw Hitem_Error("void Subs::Hdvector::set_value(std::string&) const: error translating a std::string into a dvector.");
	}
    }
    value = temp;
}

// Huchar
void Subs::Huchar::print(std::ostream& s) const {
    s << int(value);
}

void Subs::Huchar::write(std::ofstream& s) const {
    write_string(s,comment);
    s.write((char*)&value,sizeof(UCHAR));
    if(!s) throw Hitem_Error("Subs::Huchar::write(std::ofstream&): failed to write item");
}

void Subs::Huchar::write_ascii(std::ofstream& s) const {
    s << "UC " << value << " " << comment << std::endl;
    if(!s) throw Hitem_Error("Subs::Huchar::write_ascii(std::ofstream&): failed to write item");
}

void Subs::Huchar::read(std::ifstream& s, bool swap_bytes) {
    read_string(s,comment,swap_bytes);
    s.read((char*)&value,sizeof(UCHAR));
    if(!s) throw Hitem_Error("Subs::Huchar::read(std::ifstream&): failed to read item");
}

void Subs::Huchar::skip(std::ifstream& s, bool swap_bytes) {
    skip_string(s,swap_bytes);
    s.ignore(sizeof(UCHAR));
    if(!s) throw Hitem_Error("Subs::Huchar::skip(std::ifstream&): failed to skip item");
}

void Subs::Huchar::get_value(std::string &sval) const {
    std::ostringstream ostr;
    ostr << int(value);
    if(!ostr)
	throw Hitem_Error("void Subs::Huchar::get_value(std::string&) const: error translating a char into a string.");
    sval = ostr.str();
}

std::string Subs::Huchar::get_string() const {
    std::ostringstream ostr;
    ostr << int(value);
    if(!ostr)
	throw Hitem_Error("void Subs::Huchar::get_string() const: error translating a char into a string.");
    return ostr.str();
}

void Subs::Huchar::set_value(const std::string &sval){
    std::istringstream istr(sval);
    istr >> value;
    if(!istr)
	throw Hitem_Error(std::string("void Subs::Hitem::set_value(const std::string&): failed to translate \"") + sval + 
			  std::string("\" into an unsigned char"));
}

// vector<int>
void Subs::Hivector::print(std::ostream& s) const {
    s << std::setfill(' ');
    for(size_t i=0; i<value.size(); i++)
	s << value[i] << " ";
}

void Subs::Hivector::write(std::ofstream& s) const {
    write_string(s,comment);
    UINT4 n = UINT4(value.size());
    s.write((char*)&n, sizeof(UINT4));
    INT4 j;
    for(size_t i=0; i<value.size(); i++){
	j = value[i];
	s.write((char*)&j, sizeof(INT4));
    }
    if(!s) throw Hitem_Error("Subs::Hivector::write(std::ofstream&): failed to write item");
}
 
void Subs::Hivector::write_ascii(std::ofstream& s) const {
    s << "IV " << value.size();
    for(size_t i=0; i<value.size(); i++){
	s << " " << value[i];
    }
    s << " " << comment << std::endl;
    if(!s) throw Hitem_Error("Subs::Hivector::write_ascii(std::ofstream&): failed to write item");
}

void Subs::Hivector::read(std::ifstream& s, bool swap_bytes) {
    read_string(s,comment,swap_bytes);
    UINT4 n;
    s.read((char*)&n,sizeof(UINT4));
    if(swap_bytes) n = Subs::byte_swap(n);
    value.resize(n);
    INT4 j;
    for(size_t i=0; i<value.size(); i++){
	s.read((char*)&j, sizeof(INT4));
	if(swap_bytes) j = Subs::byte_swap(j);
	value[i] = j;
    }
    if(!s) throw Hitem_Error("Subs::Hivector::read(std::ifstream&, bool): failed to read item");
}

void Subs::Hivector::skip(std::ifstream& s, bool swap_bytes) {
    skip_string(s, swap_bytes);
    UINT4 n;
    s.read((char*)&n, sizeof(UINT4));
    if(swap_bytes) n = Subs::byte_swap(n);
    s.ignore(n*sizeof(INT4));
    if(!s) throw Hitem_Error("Subs::Hivector::skip(std::ifstream&): failed to skip item");
}

void  Subs::Hivector::get_value(std::string &sval) const {
    std::ostringstream ostr;
    for(std::vector<int>::size_type i=0; i<value.size(); i++){
	ostr << value[i] << " ";
	if(!ostr)
	    throw Hitem_Error("void Subs::Hivector::get_value(std::string&) const: error translating an ivector into a string.");
    }
    sval = ostr.str();
}

std::string Subs::Hivector::get_string() const {
    std::ostringstream ostr;
    for(std::vector<int>::size_type i=0; i<value.size(); i++){
	ostr << value[i] << " ";
	if(!ostr)
	    throw Hitem_Error("void Subs::Hivector::get_value(std::string&) const: error translating an ivector into a string.");
    }
    return ostr.str();
}

void Subs::Hivector::set_value(const std::string &sval){
    // Blank added to avoid apparent bug in istringstream
    std::istringstream istr(sval + " ");
    std::vector<int> temp;
    int j;
    while(!istr.eof()){
	istr >> j;
	if(istr.good()){
	    temp.push_back(j);
	}else if(!istr.eof()){
	    throw Hitem_Error("void Subs::Hivector::set_value(std::string&) const: error translating a std::string into an ivector.");
	}
    }
    value = temp;
}

// vector<float>
void Subs::Hfvector::print(std::ostream& s) const {
    s << std::setfill(' ');
    for(size_t i=0; i<value.size(); i++)
	s << value[i] << " ";
}

void Subs::Hfvector::write(std::ofstream& s) const {
    write_string(s,comment);
    UINT4 n = UINT4(value.size());
    s.write((char*)&n, sizeof(UINT4));
    REAL4 f;
    for(size_t i=0; i<value.size(); i++){
	f = value[i];
	s.write((char*)&f, sizeof(REAL4));
    }
    if(!s) throw Hitem_Error("Subs::Hfvector::write(std::ofstream&): failed to write item");
}
 
void Subs::Hfvector::write_ascii(std::ofstream& s) const {
    s << "FV " << value.size();
    for(size_t i=0; i<value.size(); i++){
	s << " " << value[i];
    }
    s << " " << comment << std::endl;
    if(!s) throw Hitem_Error("Subs::Hfvector::write_ascii(std::ofstream&): failed to write item");
}

void Subs::Hfvector::read(std::ifstream& s, bool swap_bytes) {
    read_string(s,comment,swap_bytes);
    UINT4 n;
    s.read((char*)&n,sizeof(UINT4));
    if(swap_bytes) n = Subs::byte_swap(n);
    value.resize(n);
    REAL4 f;
    for(size_t i=0; i<value.size(); i++){
	s.read((char*)&f, sizeof(REAL4));
	if(swap_bytes) f = Subs::byte_swap(f);
	value[i] = f;
    }
    if(!s) throw Hitem_Error("Subs::Hfvector::read(std::ifstream&, bool): failed to read item");
}

void Subs::Hfvector::skip(std::ifstream& s, bool swap_bytes) {
    skip_string(s, swap_bytes);
    UINT4 n;
    s.read((char*)&n,sizeof(UINT4));
    if(swap_bytes) n = Subs::byte_swap(n);
    s.ignore(n*sizeof(REAL4));
    if(!s) throw Hitem_Error("Subs::Hfvector::skip(std::ifstream&, bool): failed to skip item");
}

void  Subs::Hfvector::get_value(std::string &sval) const {
    std::ostringstream ostr;
    for(std::vector<float>::size_type i=0; i<value.size(); i++){
	ostr << value[i] << " ";
	if(!ostr)
	    throw Hitem_Error("void Subs::Hfvector::get_value(std::string&) const: error translating a fvector into a string.");
    }
    sval = ostr.str();
}

std::string Subs::Hfvector::get_string() const {
    std::ostringstream ostr;
    for(std::vector<float>::size_type i=0; i<value.size(); i++){
	ostr << value[i] << " ";
	if(!ostr)
	    throw Hitem_Error("void Subs::Hfvector::get_value(std::string&) const: error translating a fvector into a string.");
    }
    return ostr.str();
}

void Subs::Hfvector::set_value(const std::string &sval){
    // Blank added to avoid apparent bug in istringstream
    std::istringstream istr(sval + " ");
    std::vector<float> temp;
    float f;
    while(!istr.eof()){
	istr >> f;
	if(istr.good()){
	    temp.push_back(f);
	}else if(!istr.eof()){
	    throw Hitem_Error("void Subs::Hfvector::set_value(std::string&) const: error translating a std::string into a fvector.");
	}
    }
    value = temp;
}


/**
 * Reads an Hitem of any type from an ASCII file. This routine reads an initial code
 * for each type and then constructs a pointer to an Hitem based upon the type
 * code. If you create another supported Hitem type, you need to update this
 * routine. Supported codes are 'D' = double, 'C' = char, 'I' = int, 'UI' = unsigned
 * int, 'LI@ = long int, 'ULI' = unsigned long int, 'F' = float, 'S' = string, 'B' = bool, 
 * 'POS' = position, 'TEL' = telescope, 'USI' = unsigned short int.  The Hitem must occupy a line by itself.
 * \param istr the input file stream
 */

Subs::Hitem* Subs::Hitem::read_ascii_item(std::ifstream& istr){

    // First read the type
    std::string stype;
    istr >> stype;
    if(!istr) throw Hitem_Error("Subs::Hitem::read_ascii_item(std::ifstream&): failed to read type");
    stype = Subs::toupper(stype);
    std::string comment;
    if(stype == "D"){
	double value;
	istr >> value;
	getline(istr, comment);
	return new Hdouble(value, comment);
	
    }else if(stype == "C"){
	char value;
	istr >> value;
	getline(istr, comment);
	return new Hchar(value, comment);
	
    }else if(stype == "I"){
	int value;
	istr >> value;
	getline(istr, comment);
	return new Hint(value, comment);
	
    }else if(stype == "LI"){
	long int value;
	istr >> value;
	getline(istr, comment);
	return new Hlint(value, comment);
	
    }else if(stype == "UI"){
	unsigned int value;
	istr >> value;
	getline(istr, comment);
	return new Huint(value, comment);
	
    }else if(stype == "ULI"){
	unsigned long int value;
	istr >> value;
	getline(istr, comment);
	return new Hulint(value, comment);

    }else if(stype == "USI"){
	UINT2 value;
	istr >> value;
	getline(istr, comment);
	return new Husint(value, comment);
	
    }else if(stype == "F"){
	float value;
	istr >> value;
	getline(istr, comment);
	return new Hfloat(value, comment);
	
    }else if(stype == "S"){
	char c = istr.peek();
	while(istr &&  (c == ' ' || c == '\t')){
	    istr.get(c);
	    c = istr.peek();
	}
	
	std::string value;
	getline(istr, value);
	std::string::size_type pos = value.find('#');
	if(pos == std::string::npos)
	    return new Hstring(value);
	else
	    return new Hstring(value.substr(0,pos), value.substr(pos+1));
	
    }else if(stype == "B"){
	std::string value;
	istr >> value;
	value = Subs::toupper(value);
	getline(istr, comment);
	if(value == "TRUE" || value == "T" || value == "Y" || value == "YES"){
	    return new Hbool(true, comment);
	}else if(value == "FALSE" || value == "F" || value == "N" || value == "NO"){
	    return new Hbool(false, comment);
	}else{
	    throw Hitem_Error("Subs::Hitem::read_ascii_item(std::ifstream&): failed to understand boolean value = " + value);
	}
	
    }else if(stype == "DIR"){
	getline(istr, comment);
	return new Hdirectory(comment);
	
    }else if(stype == "DATE"){
	Date value;
	istr >> value;
	getline(istr, comment);
	return new Hdate(value, comment);
	
    }else if(stype == "T"){
	Time value;
	istr >> value;
	getline(istr, comment);
	return new Htime(value, comment);
	
    }else if(stype == "POS"){
	Position value;
	istr >> value;
	getline(istr, comment);
	return new Hposition(value, comment);
	
    }else if(stype == "TEL"){
	Telescope value;
	istr >> value;
	getline(istr, comment);
	return new Htelescope(value, comment);
	
    }else if(stype == "DV"){
	std::vector<double> value;
	int nvec;
	istr >> nvec;
	double dval;
	for(int i=0; i<nvec; i++){
	    istr >> dval;
	    if(!istr) value.push_back(dval);
	}
	getline(istr, comment);
	return new Hdvector(value, comment);

    }else if(stype == "IV"){
	std::vector<int> value;
	int nvec;
	istr >> nvec;
	int ival;
	for(int i=0; i<nvec; i++){
	    istr >> ival;
	    if(!istr) value.push_back(ival);
	}
	getline(istr, comment);
	return new Hivector(value, comment);

    }else if(stype == "FV"){
	std::vector<float> value;
	int nvec;
	istr >> nvec;
	float fval;
	for(int i=0; i<nvec; i++){
	    istr >> fval;
	    if(!istr) value.push_back(fval);
	}
	getline(istr, comment);
	return new Hfvector(value, comment);
	
    }else if(stype == "UC"){
	unsigned char value;
	istr >> value;
	getline(istr, comment);
	return new Huchar(value, comment);
	
    }else{
	throw Hitem_Error("Subs::Hitem::read_ascii_item(std::ifstream&): failed to understand type = " + stype);
    }
}

/**
 * Reads an Hitem of any type from a binary file. This routine reads an initial code
 * for each type and then constructs a pointer to an Hitem based upon the type
 * code. If you create another supported Hitem type, you need to update this
 * routine.
 * \param istr the input file stream opened for binary IO
 */

Subs::Hitem* Subs::Hitem::read_item(std::ifstream& istr, bool swap_bytes){

    Subs::Hitem::HITEM_TYPE otype = Subs::Hitem::read_type(istr, swap_bytes);

    // Each type must be represented here
    switch(otype) {

	case Hitem::HCHAR:
	    return new Hchar(istr, swap_bytes);

	case Hitem::HDOUBLE:
	    return new Hdouble(istr, swap_bytes);

	case Hitem::HINT:
	    return new Hint(istr, swap_bytes);

	case Hitem::HUINT:
	    return new Huint(istr, swap_bytes);

	case Hitem::HLINT:
	    return new Hlint(istr, swap_bytes);

	case Hitem::HULINT:
	    return new Hulint(istr, swap_bytes);

	case Hitem::HFLOAT:
	    return new Hfloat(istr, swap_bytes);

	case Hitem::HBOOL:
	    return new Hbool(istr, swap_bytes);

	case Hitem::HSTRING:
	    return new Hstring(istr, swap_bytes);

	case Hitem::HDIRECTORY:
	    return new Hdirectory(istr, swap_bytes);

	case Hitem::HDATE:
	    return new Hdate(istr, swap_bytes);

	case Hitem::HTIME:
	    return new Htime(istr, swap_bytes);

	case Hitem::HPOSITION:
	    return new Hposition(istr, swap_bytes);

	case Hitem::HTELESCOPE:
	    return new Htelescope(istr, swap_bytes);

	case Hitem::HDVECTOR:
	    return new Hdvector(istr, swap_bytes);

	case Hitem::HIVECTOR:
	    return new Hivector(istr, swap_bytes);

	case Hitem::HFVECTOR:
	    return new Hfvector(istr, swap_bytes);

	case Hitem::HUCHAR:
	    return new Huchar(istr, swap_bytes);

	case Hitem::HUSINT:
	    return new Husint(istr, swap_bytes);

	default:
	    throw 
		Hitem_Error("Subs::Hitem* Subs::Hitem::read_item(std::ifstream&, bool): unrecognised Subs::Hitem::HITEM_TYPE type = " + Subs::str(otype));
    }
}

/**
 * Skips an Hitem of any type in a file. This routine reads an initial code
 * for each type and then skips the appropriate number of bytes based upon this
 * type.. If you create another supported Hitem type, you need to update this
 * routine.
 * \param istr the input file stream opened for binary IO
 */

void Subs::Hitem::skip_item(std::ifstream& istr, bool swap_bytes){

    Subs::Hitem::HITEM_TYPE otype = Subs::Hitem::read_type(istr, swap_bytes);

    // Each type must be represented here

    switch(otype) {

	case Hitem::HCHAR:
	    Hchar::skip(istr, swap_bytes);
	    return;

	case Hitem::HDOUBLE:
	    Hdouble::skip(istr, swap_bytes);
	    return;

	case Hitem::HINT:
	    Hint::skip(istr, swap_bytes);
	    return;

	case Hitem::HUINT:
	    Huint::skip(istr, swap_bytes);
	    return;

	case Hitem::HLINT:
	    Hlint::skip(istr, swap_bytes);
	    return;

	case Hitem::HULINT:
	    Hulint::skip(istr, swap_bytes);
	    return;

	case Hitem::HFLOAT:
	    Hfloat::skip(istr, swap_bytes);
	    return;

	case Hitem::HBOOL:
	    Hbool::skip(istr, swap_bytes);
	    return;

	case Hitem::HSTRING:
	    Hstring::skip(istr, swap_bytes);
	    return;

	case Hitem::HDIRECTORY:
	    Hdirectory::skip(istr, swap_bytes);
	    return;

	case Hitem::HDATE:
	    Hdate::skip(istr, swap_bytes);
	    return;

	case Hitem::HTIME:
	    Htime::skip(istr, swap_bytes);
	    return;

	case Hitem::HPOSITION:
	    Hposition::skip(istr, swap_bytes);
	    return;

	case Hitem::HTELESCOPE:
	    Htelescope::skip(istr, swap_bytes);
	    return;

	case Hitem::HDVECTOR:
	    Hdvector::skip(istr, swap_bytes);
	    return;

	case Hitem::HIVECTOR:
	    Hivector::skip(istr, swap_bytes);
	    return;

	case Hitem::HFVECTOR:
	    Hfvector::skip(istr, swap_bytes);
	    return;

	case Hitem::HUCHAR:
	    Huchar::skip(istr, swap_bytes);
	    return;

	case Hitem::HUSINT:
	    Husint::skip(istr, swap_bytes);
	    return;

	default:
	    throw 
		Hitem_Error("void Subs::Hitem::skip_item(std::ifstream&): unrecognised Subs::Hitem::HITEM_TYPE type = " + str(otype));
    }
}

/** Writes out the Hitem pointed to by item to a disk file
 * \param ostr output file opened for binary IO
 * \param item pointer to Hitem to write out
 */

void Subs::Hitem::write_item(std::ofstream& ostr, Hitem* item) {
    item->write_type(ostr);
    item->write(ostr);
}

