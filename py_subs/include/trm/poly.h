#ifndef TRM_POLY
#define TRM_POLY

#include "trm/subs.h"
#include "trm/array1d.h"

namespace Subs {

  //! A class for representing polynomials 
  /** Poly defines a class to be used for representing scales
   * when rebinning. The poly is specified by giving a start and end value
   * that it is scaled between. Let middle=(start+end)/2 and hrange=(end-start)/2
   * then the value is given by
   * sum_{n=0}^{npoly-1} poly coeff(n)*((x-middle)/hrange))^n
   * This avoids any chance of overflow, although there could easily be underflow.
   * There is also an option to take the exponential of the above number to give 
   * a logarithmic scale in the case of npoly=2 
   */

  class Poly : public Array1D<double> {

  public:

    //! Default constructor
    Poly() : Array1D<double>(), norm(true), middle(0.), hrange(1.) {}

    //! General constructor
    Poly(int npoly, bool normal, double xs, double xe);

    //! Another general constructor
    Poly(bool normal, double xs, double xe, const Array1D<double>& coeff);

    //! Another general constructor
    Poly(bool normal, double xs, double xe, const std::vector<double>& coeff);

    //! Constructor of linear poly running from xs at pixel -0.5 to xe at pixel npix-0.5
    Poly(double xs, double xe, int npix);

    //! Copy constructor
    Poly(const Poly& obj);

    //! Constructor of linear pixel scale running 1 to npix
    Poly(int npix);

    //! Compute the value of the poly at x
    double get_value(double x) const;

    //! Compute the derivative of the poly at x
    double get_deriv(double x) const;

    //! Compute the value of x for a given poly value (not always soluble)
    double get_x(double value, double xguess, double acc) const;

    //! Returns whether poly is 'normal' as opposed to exponential in character
    bool is_normal() const { return norm;}

    //! Write out a poly to a binary file
    void write(std::ofstream& s) const;

    //! Skip a poly in a binary file
    static void skip(std::ifstream& s, bool swap_bytes);

    //! Read a poly from a binary file
    void read(std::ifstream& s, bool swap_bytes);

    //! ASCII output
    friend std::ostream& operator<<(std::ostream& s, const Subs::Poly& poly);

    //! ASCII input
    friend std::istream& operator>>(std::istream& s, Subs::Poly& poly);

  private:

    bool norm;
    double middle, hrange;

  };


};

#endif
