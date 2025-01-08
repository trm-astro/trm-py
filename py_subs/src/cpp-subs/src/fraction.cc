#include "trm/fraction.h"

/** Constructor of general fraction. 
 * \param numerator the top of the fraction
 * \param denominator the bottom of the fraction
 */
Subs::Fraction::Fraction(int numerator, int denominator) : num(numerator), denom(denominator) {
    if(denom <= 0)
	throw Fraction_Error("Subs::Fraction(int, int): denominator <= 0 in constructor");
    this->reduce();
}

/** Constructor from strings of the form -1/3 or 5/6 etc
 * \param fstring the string to set from
 */
Subs::Fraction::Fraction(const std::string& fstring){
    std::istringstream istr(fstring);
    char c;
    istr >> num;
    if(!istr)
	throw Fraction_Error("Subs::Fraction(std::string&): failed to read string");
    istr >> c >> denom;
    if(!istr){
	denom = 1;
    }else{
	if(denom <= 0)
	    throw Fraction_Error("Subs::Fraction(std::string&): denominator <= 0 in constructor");
    }
    this->reduce();
}  

/** Set equal to an integer.
 * \param number the integer to set the fraction equal to.
 */
Subs::Fraction& Subs::Fraction::operator=(int number){
    num   = number;
    denom = 1;
    return *this;
}

/** Returns negative of a fraction
 */
Subs::Fraction Subs::Fraction::operator-() const {
    return Fraction(-num, denom);
}

/** This function comes back with true or false according to whether two
 * fractions are equal. 
 * \param frac1 the first fraction
 * \param frac2 the second fraction
 * \return true if the two units are equal
 */
bool Subs::operator==(const Subs::Fraction& frac1, const Subs::Fraction& frac2){

    if(frac1.denom > frac2.denom){
	if(frac1.denom % frac2.denom != 0 || frac1.num % frac2.num != 0) return false;
	if(frac1.num / frac2.num == frac1.denom / frac2.denom)
	    return true;
	else
	    return false;
    }else if(frac1.denom < frac2.denom){
	if(frac2.denom % frac1.denom != 0 || frac2.num % frac1.num != 0) return false;
	if(frac2.num / frac1.num == frac2.denom / frac1.denom)
	    return true;
	else
	    return false;
    }else{
	if(frac1.num == frac2.num)
	    return true;
	else
	    return false;
    }
}

/** This function comes back with true or false according to whether a fraction
 * equals an integer
 * \param frac the  fraction
 * \param constant the integer
 * \return true if the two units are equal
 */
bool Subs::operator==(const Subs::Fraction& frac, int number){
    return frac.num == number*frac.denom;
}

/** This function comes back with true or false according to whether a fraction
 * does not equal an integer
 * \param frac the  fraction
 * \param constant the integer
 * \return true if the two units are equal
 */
bool Subs::operator!=(const Subs::Fraction& frac, int number){
    return frac.num != number*frac.denom;
}

/** Checks whether a given fraction is greater than an integer
 * \param frac the fraction to compare
 * \param number the integer to compare it against
 * \return true if fraction is greater than value
 */
bool Subs::operator>(const Subs::Fraction& frac, int number){
    return frac.num > number*frac.denom;
}

/** Checks whether a given fraction is less than an integer
 * \param frac the fraction to compare
 * \param number the integer to compare it against
 * \return true if fraction is less than integer
 */
bool Subs::operator<(const Subs::Fraction& frac, int number){
    return frac.num < number*frac.denom;
}

/** Adds fractions.
 * \param frac1 the first fraction
 * \param frac2 the second fraction
 * \return the combined fraction
 */
Subs::Fraction Subs::operator+(const Subs::Fraction& frac1, const Subs::Fraction& frac2){
    return Fraction(frac1.num*frac2.denom + frac2.num*frac1.denom,  frac1.denom * frac2.denom);
}

/** Subtracts fractions.
 * \param frac1 the first fraction
 * \param frac2 the second fraction
 * \return the subtracted fraction
 */
Subs::Fraction Subs::operator-(const Subs::Fraction& frac1, const Subs::Fraction& frac2){
    return Fraction(frac1.num*frac2.denom - frac2.num*frac1.denom,  frac1.denom * frac2.denom);
}

// Reduces fractions. i.e. it will convert a fraction such as 6/24 into 1/4.
// This is in general a very difficult task so this only divides by primes up to
//
void Subs::Fraction::reduce(){

    if(num == 0){
	denom = 1;
	return;
    }

    const int NPRIME = 20;
    int prime[NPRIME] = {2,3,5,7,11,13,17,19,23,29,31,37,41,43,47,53,59,61,67,71};

    int nred = 1;
    while(nred > 0){
	nred = 0;
	for(int i=0; i<NPRIME; i++){
	    if(num % prime[i] == 0 && denom % prime[i] == 0){
		num   /= prime[i];
		denom /= prime[i];
		nred++;
	    }
	}
    }
}

/** This function returns the name of units suited for plotting with
 * the PGPLOT package. It attempts to recognise simple forms of units,
 * but failing that it will come back with a full string.
 */
std::string Subs::Fraction::get_string() const {
    if(denom == 1)
	return Subs::str(num);
    else
	return Subs::str(num) + "/" + Subs::str(denom);
}

/** Writes a Fraction out in binary format
 * \param ostr the output stream
 */
void Subs::Fraction::write(std::ofstream& ostr) const {
    ostr.write((char*)&num,   sizeof(int));  
    ostr.write((char*)&denom, sizeof(int));
    if(!ostr) throw Fraction_Error("Subs::Fraction::write(std::ofstream&): failure during write");
}

/** Reads a Fraction in binary format
 * \param istr the input stream
 */
void Subs::Fraction::read(std::ifstream& istr, bool swap_bytes) {
    istr.read((char*)&num,   sizeof(int));  
    if(swap_bytes) num = Subs::byte_swap(num);
    istr.read((char*)&denom, sizeof(int));
    if(swap_bytes) denom = Subs::byte_swap(denom);
    if(!istr) throw Fraction_Error("Subs::Fraction::read(std::ifstream&): failure during read");
}

/** Skips a Fraction in binary format
 * \param istr the input stream
 */
void Subs::Fraction::skip(std::ifstream& istr) {
    istr.ignore(sizeof(int));  
    istr.ignore(sizeof(int));
    if(!istr) throw Fraction_Error("Subs::Fraction::skip(std::ifstream&): failure during skip");
}

/** Outputs a fraction as -1/3 or whatever
 * \param ostr the output stream
 * \param frac the fraction
 */
std::ostream& Subs::operator<<(std::ostream& ostr, const Subs::Fraction& frac){
    if(frac.denom == 1)
	ostr << frac.num;
    else
	ostr << frac.num << "/" << frac.denom;
    return ostr;
}
