#ifndef TRM_BINARY
#define TRM_BINARY

#include <string>
using namespace std;

//! Binary software namespace

/*
 * This is the namespace for all the \c binary software for things to do with
 * binary stars. This is mainly for use in studies of the evolution of binary
 * stars and provides functions for evaluating the usual stuff such as orbital
 * periods, separations etc, which one tends to need often. Circular orbits are
 * assumed throughout.
 */

namespace Binary{

  //! Default directory for command defaults
  const char BINARY_DIR[] = ".binary";

  //! Environment variable for switching directory for command defaults
  const char BINARY_ENV[] = "BINARY_ENV";

  //! An exception class.

  /** Binary::Binary_Error is the error class for the Binary programs.
   * It is inherited from the standard string class.
   */
  class Binary_Error : public std::string {
  public:

    //! Default constructor
    Binary_Error() : std::string() {}

    //! Constructor storing a message
    Binary_Error(const std::string& err) : std::string(err) {} 
  };

  //! Calculates whether a system of two spheres eclipses
  bool eclipses(double m1, double m2, double r1, double r2, 
		double period, double iangle, double& width);

  //! Computes angular momentum loss due to GR
  double jdotgr(double m1, double m2, double a);

  //! Computes radius of white dwarf according to Eggleton
  double mr_wd_eggleton(double m);

  //! Computes radius of main sequence star
  double mr_main_sequence(double m);

  //! Computes radius of white dwarf according to Nauenberg
  double mr_wd_nauenberg(double m);

  //! Computes Keplerian angular velocity
  double kepler_omega(double m, double r);

  //! Computes orbital angular momentum
  double orbital_ang_mom(double m1, double m2, double a);
  
  //! Computes orbital angular velocity
  double orbital_omega(double m1, double m2, double a);

  //! Computes orbital period
  double orbital_period(double m1, double m2, double a);

  //! Computes orbital separation
  double orbital_separation(double m1, double m2, double period);

  //! Computes scaled circularisation radius
  double rcirc(double q);

  //! Computes derivative of scaled circularisation radius wrt q
  double drcircdq(double q);

  //! Computes minimum radius
  double rmin(double q);

  //! Computes d log R / d log M for main sequence star
  double zeta_main_sequence(double m);

  //! Computes d log R / d log M for Eggleton WD formula
  double zeta_wd_eggleton(double m);

  //! Computes d log R / d log M for Nauenberg WD formula
  double zeta_wd_nauenberg(double m);
  
  //! Computes d zeta / d log M for Eggleton WD formula
  double dzetadm_wd_eggleton(double m);

  //! Computes moment of inertia constant 'k'
  double gyrate(double m);

  //! Computes d log (k) / d log(M) for moment of inertia constant 'k'
  double zeta_gyrate(double m);

  //! Structure for information on elliptical orbit stuff
  struct Einfo {
    double vx;      /**< X velocity km/s **/
    double vy;      /**< Y velocity km/s **/
    double vx_x;    /**< partial derivative of X velocity wrt x **/
    double vx_y;    /**< partial derivative of X velocity wrt y **/
    double vy_x;    /**< partial derivative of Y velocity wrt x **/
    double vy_y;    /**< partial derivative of Y velocity wrt y **/
    double r;       /**< radius, solar */
    double omega;   /**< angular velocity, rads/sec */
    double E;       /**< Eccentric anomoly, radians **/
  };

  //! Computes x- and y-velocity for elliptical orbits
  Einfo vellipse(double m, double e, double a, double theta);

}

#endif
