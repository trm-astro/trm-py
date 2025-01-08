#ifndef TRM_VEC3_H
#define TRM_VEC3_H

#include <iostream>
#include "trm/subs.h"

namespace Subs{

    //! Class for 3-vectors

    /** This class represent  Cartesian 3-vectors. It works in double precision for
     * accuracy, and defines dot and cross-products etc. It is done largely
     * inline for speed.
     */

    class Vec3 {
    public:

	//! Default constructor (zero vector)
	Vec3() : xc(0.), yc(0.), zc(0.) {}

	//! General constructor from three numbers
	Vec3(double x, double y, double z) : xc(x), yc(y), zc(z) {}

	//! General constructor 3 element array
	Vec3(double* v) : xc(v[0]), yc(v[1]), zc(v[2]) {}

	//! Access to x
	double& x(){return xc;}

	//! Access to y
	double& y(){return yc;}

	//! Access to z
	double& z(){return zc;}

	//! Access to x
	const double& x() const {return xc;}

	//! Access to y
	const double& y() const {return yc;}

	//! Access to z
	const double& z() const {return zc;}

	//! Set to 3 numbers
	void set(double x, double y, double z){xc=x; yc=y; zc=z;}
 
	//! Set to a 3 element array
	void set(double* v){xc=v[0]; yc=v[1]; zc=v[2];}

	//! Get to a 3 element array
	void get(double* v) const {v[0]=xc; v[1]=yc; v[2]=zc;}

	//! Normalises a vector
	void unit(){
	    double norm = Subs::sqr(xc)+Subs::sqr(yc)+Subs::sqr(zc);
	    if(norm > 0.){
		norm = sqrt(norm);
		xc /= norm;
		yc /= norm;
		zc /= norm;
	    }else{
		throw Subs_Error("void Subs::Vec3::unit(): null vector");
	    }
	}
    
	//! Returns the length of a vector
	double length() const {
	    return sqrt(Subs::sqr(xc)+Subs::sqr(yc)+Subs::sqr(zc));
	}
    
	//! Returns the length squared of a vector
	double sqr() const {
	    return Subs::sqr(xc)+Subs::sqr(yc)+Subs::sqr(zc);
	}
    
	//! Sets a vector to be a unit vector in the X direction
	void unitx(){xc=1.; yc=zc=0.;}
    
	//! Sets a vector to be a unit vector in the Y direction
	void unity(){yc=1.; xc=zc=0.;}
    
	//! Sets a vector to be a unit vector in the Z direction
	void unitz(){zc=1.; xc=yc=0.;}
    
	//! Returns scalar or dot product of two vectors
	friend double dot(const Vec3& vec1, const Vec3& vec2){
	    return (vec1.xc*vec2.xc+vec1.yc*vec2.yc+vec1.zc*vec2.zc);
	}
    
	//! Returns cross product of two vectors
	friend Vec3 cross(const Vec3& vec1, const Vec3& vec2){
	    Vec3 temp;
	    temp.xc = vec1.yc*vec2.zc-vec1.zc*vec2.yc;
	    temp.yc = vec1.zc*vec2.xc-vec1.xc*vec2.zc;
	    temp.zc = vec1.xc*vec2.yc-vec1.yc*vec2.xc;
	    return temp;
	}
    
	//! Multiplication by a constant in place
	void operator*=(double con){xc *= con; yc *= con; zc *= con;}
    
	//! Division by a constant in place
	void operator/=(double con){xc /= con; yc /= con; zc /= con;}
    
	//! Addition of another Vec3 in place
	void operator+=(const Vec3& vec){xc += vec.xc; yc += vec.yc; zc += vec.zc;}
    
	//! Subtraction of another Vec3 in place
	void operator-=(const Vec3& vec){xc -= vec.xc; yc -= vec.yc; zc -= vec.zc;}
    
	//! Difference between two vectors
	friend Vec3 operator-(const Vec3& vec1, const Vec3& vec2){
	    Vec3 temp = vec1;
	    temp -= vec2;
	    return temp;
	}

	//! Sum of two vectors
	friend Vec3 operator+(const Vec3& vec1, const Vec3& vec2){
	    Vec3 temp = vec1;
	    temp += vec2;
	    return temp;
	}

	//! Multiplication by a constant
	friend Vec3 operator*(double con, const Vec3& vec){Vec3 temp = vec; temp *= con; return temp;}

	//! Division by a constant
	friend Vec3 operator/(const Vec3& vec, double con){Vec3 temp = vec; temp /= con; return temp;}

	//! ASCII output
	friend std::ostream& operator<<(std::ostream& ostr, const Vec3& vec){
	    ostr << vec.xc << " " << vec.yc << " " << vec.zc;
	    return ostr;
	}

	//! ASCII input
	friend std::istream& operator>>(std::istream& istr, Vec3& vec){
	    istr >> vec.xc >> vec.yc >> vec.zc;
	    return istr;
	}

	//! negative of a vector
	Vec3 operator-(){
	    Vec3 temp; 
	    temp.xc = -xc; 
	    temp.yc = -yc; 
	    temp.zc = -zc; 
	    return temp;
	}

    private:
	double xc, yc, zc;
    };

    double dot(const Vec3& vec1, const Vec3& vec2);

    Vec3 cross(const Vec3& vec1, const Vec3& vec2);

    Vec3 operator-(const Vec3& vec1, const Vec3& vec2);

    Vec3 operator+(const Vec3& vec1, const Vec3& vec2);

    Vec3 operator*(double con, const Vec3& vec);

    Vec3 operator/(const Vec3& vec, double con);

    std::ostream& operator<<(std::ostream& ostr, const Vec3& vec);
  
    std::istream& operator>>(std::istream& istr, Vec3& vec);
};

#endif
