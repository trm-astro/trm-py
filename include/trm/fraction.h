#ifndef TRM_FRACTION
#define TRM_FRACTION

#include <cmath>
#include <string>
#include "trm/subs.h"
#include "trm/constants.h"

namespace Subs {

    //! A class to represent fractions

    /**
     * The Subs::Fraction class is a helper class for the Subs::Units class. It represents
     * fractions such as 1/3 or -2/13. Internally it stores integers for exact precision.
     */

    class Fraction {
    public:

	//! Default constructor
	Fraction() : num(1), denom(1) {}

	//! General constructor
	Fraction(int numerator, int denominator);

	//! Constructor from an integer
	Fraction(int number) : num(number), denom(1) {}

	//! Constructor from a string
	Fraction(const std::string& fstring);

	//! Set equal to an integer
	Fraction& operator=(int number);

	//! Unary - operator
	Fraction operator-() const;

	friend bool operator==(const Subs::Fraction& frac1, const Subs::Fraction& frac2);

	friend bool operator==(const Subs::Fraction& frac, int number);

	friend bool operator!=(const Subs::Fraction& frac, int number);

	friend bool operator>(const Subs::Fraction& frac, int number);

	friend bool operator<(const Subs::Fraction& frac, int number);

	friend Fraction operator+(const Subs::Fraction& frac1, const Subs::Fraction& frac2);

	friend Fraction operator-(const Subs::Fraction& frac1, const Subs::Fraction& frac2);

	friend std::ostream& operator<<(std::ostream& ostr, const Subs::Fraction& frac);

	//! Returns string representing the fraction
	std::string get_string() const;

	//! Binary write
	void write(std::ofstream& ostr) const;

	//! Binary read
	void read(std::ifstream& istr, bool swap_bytes);

	//! Binary skip
	static void skip(std::ifstream& istr);

	//! Exception class for Units
	class Fraction_Error : public Subs::Subs_Error {
	public:
	    //! Default constructor
	    Fraction_Error() : Subs::Subs_Error() {};
	    //! Constructor from a string
	    Fraction_Error(const std::string& str) : Subs::Subs_Error(str) {}
	};

    private:

	int num; // numerator
	int denom; // denominator

	// Attempts to reduces a fraction as far as possible
	void reduce();

    };

    //! test equality 
    bool operator==(const Fraction& frac1, const Fraction& frac2);

    //! Does a fraction equal an integer?
    bool operator==(const Fraction& frac, int number);

    //! Does a fraction not equal an integer?
    bool operator!=(const Fraction& frac, int number);

    //! Greater than
    bool operator>(const Fraction& frac, int number);

    //! Less than
    bool operator<(const Fraction& frac, int number);

    //! Add fractions
    Fraction operator+(const Fraction& frac1, const Fraction& frac2);

    //! Subtract fractions
    Fraction operator-(const Fraction& frac1, const Fraction& frac2);

    //! ASCII output
    std::ostream& operator<<(std::ostream& ostr, const Fraction& frac);  
};

#endif










