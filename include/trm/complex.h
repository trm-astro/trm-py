#include "trm/subs.h"

namespace Subs {

    /** Class representing complex numbers
     */
    class Complex {
    public:
	
	//! Default constructor. No initialisation
	Complex() : real_(), imag_() {}
	
	//! Sets complex number from a pure real number
	Complex(double real) : real_(real), imag_(0.) {}
	
	//! General constructor
	Complex(double real, double imag) : real_(real), imag_(imag) {}
	
	//! Modulus of a complex number
	double modulus() const {return sqrt(sqr(real_) + sqr(imag_));}
	
	double real() const {return real_;}

	double imag() const {return imag_;}

	//! Multiplication in place
	void operator*=(const Complex& c){
	    real_ = real_*c.real_ - imag_*c.imag_;
	    imag_ = real_*c.imag_ + imag_*c.real_;
	}

	//! Multiplication in place
	void operator*=(double c){
	    real_ = real_*c;
	    imag_ = imag_*c;
	}

	friend Complex operator*(const Complex& c1, const Complex& c2);

    private:
	double real_;
	double imag_;
    };

    Complex operator*(const Complex& c1, const Complex& c2);
};
