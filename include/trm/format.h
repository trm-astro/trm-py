#ifndef TRM_FORMAT
#define TRM_FORMAT

#include <iostream>
#include <iomanip>
#include <sstream>
#include "trm/subs.h"

namespace Subs {

  class Bound_form_d;
  class Bound_form_f;
  class Bound_form_s;

  /** Class to help with output formatting of floating point numbers. 
   * The idea is that this can guarantee the output format regardless of
   * what has been set before. Code using this could look like this:
   *
   * Format form;
   * form.precision(8);
   * form.scientific();
   * double d = 1.23456789;
   * cout << form(d) << '\n';
   *
   */

  class Format {

  public:

    //! Constructor
    explicit Format(int precision=8);

    //! Create and object with a format bound to a double
    Bound_form_d operator()(double d) const;

    //! Create and object with a format bound to a float
    Bound_form_f operator()(float  f) const;

    //! Create and object with a format bound to a string
    Bound_form_s operator()(const std::string& s) const;

    //! Create and object with a format bound to a string, changing width on the fly
    Bound_form_s operator()(const std::string& s, int width);
 
    //! Sets the precision
    void precision(int p);

    //! Sets scientific notation
    void scientific();

    //! Sets fixed format
    void fixed();

    //! Sets best variable output format
    void general();

    //! Sets minimum field width
    void width(int w);

    //! Print decimal point
    void showpoint();

    //! Do not print decimal point unless necessary
    void noshowpoint();

    //! Makes anything uppercase it can
    void uppercase();
    
    //! Makes anything lowercase it can
    void lowercase();

    //! Sets the fill character
    void fill(char c);

    //! Pads after value
    void left();

    //! Pads before value 
    void right();

    //! Pads between number and its sign
    void internal();

    //! Output of double
    friend std::ostream& operator<<(std::ostream& ostr, const Bound_form_d& bf);

    //! Output of float
    friend std::ostream& operator<<(std::ostream& ostr, const Bound_form_f& bf);

    //! Output of string
    friend std::ostream& operator<<(std::ostream& ostr, const Bound_form_s& bf);

  private:

    int precision_;
    int width_;
    std::ios::fmtflags format_;
    bool upper;
    bool showpnt;
    char fill_char;
    std::ios::fmtflags fadjust;
  };
  
  //! Combination of format and a double
  struct Bound_form_d {
    //! the format
    const Format& form;
    //! the double value
    double val;
    //! Constructor
    Bound_form_d(const Format& format, double value) : form(format), val(value) { }
  };

  struct Bound_form_f {
    //! the format
    const Format& form;
    //! the double value
    float val;
    //! Constructor
    Bound_form_f(const Format& format, float value) : form(format), val(value) { }
  };

  struct Bound_form_s {
    //! the format
    const Format& form;
    //! the string
    std::string val;
    //! Constructor
    Bound_form_s(const Format& format, const std::string& value) : form(format), val(value) { }
  };

  std::ostream& operator<<(std::ostream& ostr, const Subs::Bound_form_d& bf);

  std::ostream& operator<<(std::ostream& ostr, const Subs::Bound_form_f& bf);

  std::ostream& operator<<(std::ostream& ostr, const Subs::Bound_form_s& bf);

};

#endif



















