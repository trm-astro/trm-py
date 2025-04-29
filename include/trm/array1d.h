#ifndef TRM_SUBS_ARRAY1D
#define TRM_SUBS_ARRAY1D

#include <cmath>
#include <new>
#include <iostream>
#include <string>
#include <fstream>
#include <vector>
#include <algorithm>
#include "trm/subs.h"
#include "plplot.h"

namespace Subs {

  //! Template 1D array class

  /**
   * Subs::Array1D<X> is a template 1D numerical array class based upon
   * the memory handling object Buffer1D. This class should be used if you want
   * to add, subtract, multiply etc vectors (and the data type X can support it).
   * This typically means that you can use Array1D for X=float, double etc but not
   * more complex types. For the latter Buffer1D is the usual choice.
   */

  template <class X> 
  class Array1D : public Buffer1D<X> {
    
  public:

    //! Default constructor
    Array1D() : Buffer1D<X>() {}

    //! Constructor of an npix element vector
    Array1D(int npix) : Buffer1D<X>(npix) {}

    //! Constructs an npix element vector but with memory up to nmem
    Array1D(int npix, int nmem) : Buffer1D<X>(npix, nmem) {}

    //! Copy constructor
    Array1D(const Array1D& obj) : Buffer1D<X>(obj) {}

    //! Constructs an Array from a vector 
    template <class Y> 
    Array1D(const std::vector<Y>& vec);

    //! Constructor from an ordinary array
    template <class Y> 
    Array1D(int npix, const Y* vec);

    //! Constructor from a function of the array index
    template <class Y> 
    Array1D(int npix, const Y& func);

    //! Constructor from a file
    Array1D(const std::string& file) : Buffer1D<X>(file) {};

    //! Assign to a constant
    Array1D<X>& operator=(const X& con);

    //! Addition of a constant, in place
    void operator+=(const X& con);

    //! Subtraction of a constant, in place
    void operator-=(const X& con);

    //! Division by a constant, in place
    void operator/=(const X& con);

    //! Multiplication by a constant, in place
    void operator*=(const X& con);

    //! Addition of another array, in place
    template <class Y>
    void operator+=(const Array1D<Y>& vec);

    //! Subtraction of another array, in place
    template <class Y>
    void operator-=(const Array1D<Y>& vec);

    //! Multiplication by another array, in place
    template <class Y>
    void operator*=(const Array1D<Y>& vec);

    //! Division by another array, in place
    template <class Y>
    void operator/=(const Array1D<Y>& vec);

    //! Returns maximum value
    X max() const;

    //! Returns minimum value
    X min() const;

    //! Takes cosine of array
    void cos();

    //! Takes sine of array
    void sin();

    //! Determines whether values are monotonic
    bool monotonic() const;

    //! Locate a value in an ordered array
    void hunt(const X& x, int& jhi) const;

    //! Locate a value in an ordered array
    unsigned long locate(const X& x) const;

    //! Sorts an array into ascending order and returns a key to the original order
    Buffer1D<unsigned long int> sort();

    //! Return percentile (pcent from 0 to 100)
    X centile(double pcent);

    //! Return value of k-th smallest element (scrambles element order!)
    X select(int k);

    //! Returns median (scrambles element order!)
    X median();

    //! Returns sum
    X sum() const;

    //! Returns mean
    X mean() const;

    //! Returns length in Euclidean sense
    X length() const;

  };

  //! Error class
  class Array1D_Error : public Subs_Error {
  public:
    Array1D_Error() : Subs_Error("") {};
    Array1D_Error(const std::string& str) : Subs_Error(str) {};
  };
    
  /** Constructor of an Array1D from the STL vector class. It must be
   * possible to assign the two types to each other.
   * \param vec the vector to construct from. The Array1D will have the same number
   * of elements with the same values as the vector.
   */
  template <class X>
  template <class Y> 
  Array1D<X>::Array1D(const std::vector<Y>& vec) : Buffer1D<X>(vec.size()) {
    for(int i=0; i<this->size(); i++)
      this->buff[i] = vec[i];
  }

  /** Constructor of an Array1D from a standard C-like array. It must be
   * possible to assign the two types to each other.
   * \param vec the vector to construct from. The Array1D will have the same number
   * of elements with the same values as the vector.
   */
  template <class X>
  template <class Y> 
  Array1D<X>::Array1D(int npix, const Y* vec) : Buffer1D<X>(npix) {
    for(int i=0; i<this->size(); i++)
      this->buff[i] = vec[i];
  }

  /** Constructor of an Array1D from a function of the pixel index, with the first
   * pixel = 0.
   * \param func the function to construct from. The function must support a call of the
   * form func(int i) which returns the value at pixel i. Here is an example of a function (as
   * a function object class) that adds a simple offset
   * class Func {
   * public:
   * Func(int off) : offset(off) {}
   * int operator()(int i) const {return i+offset;}
   * private:
   * int offset;
   *};
   */
  template <class X>
  template <class Y> 
  Array1D<X>::Array1D(int npix, const Y& func) : Buffer1D<X>(npix) {
    for(int i=0; i<this->size(); i++)
      this->buff[i] = func(i);
  }

  // Operations with constants  
  template <class X>
  void Array1D<X>::operator+=(const X& con){
    for(int i=0; i<this->size(); i++)
      this->buff[i] += con;
  }

  template <class X>
  void Array1D<X>::operator-=(const X& con){
    for(int i=0; i<this->size(); i++)
      this->buff[i] -= con;
  }

  template <class X>
  void Array1D<X>::operator*=(const X& con){
    for(int i=0; i<this->size(); i++)
      this->buff[i] *= con;
  }

  template <class X>
  void Array1D<X>::operator/=(const X& con){
    for(int i=0; i<this->size(); i++)
      this->buff[i] /= con;
  }

  // Returns maximum value
  template <class X>
  X Array1D<X>::max() const {
    if(this->size() == 0)
      throw Array1D_Error("X Array1D<X>::max(): null array, cannot take maximum");
    
    X vmax = this->buff[0];
    for(int i=1; i<this->size(); i++)
      vmax = vmax < this->buff[i] ? this->buff[i] : vmax;
    return vmax;
  }

  // Returns minimum value
  template <class X>
  X Array1D<X>::min() const {
    if(this->size() == 0)
      throw Array1D_Error("X Array1D<X>::min(): null array, cannot take minimum");
    
    X vmin = this->buff[0];
    for(int i=1; i<this->size(); i++)
      vmin = vmin > this->buff[i] ? this->buff[i] : vmin;
    return vmin;
  }

  /** This function returns the value of the k-th smallest element of an Array1D. It uses
   * a routine that scrambles the order for speed. Copy the Array1D first if the 
   * order is important to you. The operation can be done by sorting, but this is faster
   * if you only want a single value.
   */
  template <class X>
  X Array1D<X>::select(int k) {
    return Subs::select(this->buff, this->size(), k);
  }

  /** This function returns the value of a percentile of an Array1D. It uses
   * a routine that scrambles the order for speed. Copy the Array1D first if the 
   * order is important to you. The operation can be done by sorting, but this is faster
   * if you only want a single value.
   */
  template <class X>
  X Array1D<X>::centile(double pcent) {
      int k = int(pcent/100*this->size());
      k = k < 0 ? 0 : k;
      k = k < this->size() ? k : this->size() - 1;
      return Subs::select(this->buff, this->size(), k);
  }

  /** This function computes the median. It uses 'select' which scrambles the order. 
   * Copy the Array1D first if the order is important to you. The median is clearly defined
   * for odd numbers of elements; for even numbers this routine averages the two middle values
   * (thus requiring two calls to 'select'and therefore slower than a similar odd number case.
   */
  template <class X>
  X Array1D<X>::median() {
    if(this->size() % 2 == 0){
      return (Subs::select(this->buff, this->size(), this->size()/2-1) +
	      Subs::select(this->buff, this->size(), this->size()/2))/2;
    }else{
      return Subs::select(this->buff, this->size(), this->size()/2);
    }
  }

  /** This function computes the mean. It returns zero if there are no elements
   * 
   */
  template <class X>
  X Array1D<X>::mean() const {
    X sum = 0;
    if(this->size()){
      for(int i=0; i<this->size(); i++)
	sum += this->buff[i];
      sum /= this->size();
    }
    return sum;
  }

  /** This function computes the sum. It returns zero if there are no elements
   * 
   */
  template <class X>
  X Array1D<X>::sum() const {
      X sum = 0;
      if(this->size()){
	  for(int i=0; i<this->size(); i++)
	      sum += this->buff[i];
      }
      return sum;
  }

  /** This function computes the length as the square root of the sum of squares
   */
  template <class X>
  X Array1D<X>::length() const {
    X sum = 0;
    for(int i=0; i<this->size(); i++)
      sum += Subs::sqr(this->buff[i]);
    return sqrt(sum);
  }

  // Operations with other arrays
  template <class X>
  template <class Y>
  void Array1D<X>::operator+=(const Array1D<Y>& vec){
    if(this->size() != vec.size())
      throw Array1D_Error("void Array1D<X>::operator+=(const Array1D<Y>& vec): incompatible numbers of elements, " + Subs::str(this->size()) + 
			  " versus " + Subs::str(vec.size()) );
    for(int i=0; i<this->size(); i++)
      this->buff[i] += vec.buff[i];
  }

  template <class X>
  template <class Y>
  void Array1D<X>::operator-=(const Array1D<Y>& vec){
    if(this->size() != vec.size())
      throw Array1D_Error("void Array1D<X>::operator-=(const Array1D<Y>& vec): incompatible numbers of elements, " + Subs::str(this->size()) + 
			  " versus " + Subs::str(vec.size()));
    for(int i=0; i<this->size(); i++)
      this->buff[i] -= vec.buff[i];
  }

  template <class X>
  template <class Y>
  void Array1D<X>::operator*=(const Array1D<Y>& vec){
    if(this->size() != vec.size())
      throw Array1D_Error("void Array1D<X>::operator*=(const Array1D<Y>& vec): incompatible numbers of elements, " + Subs::str(this->size()) + 
			  " versus " + Subs::str(vec.size()) );
    for(int i=0; i<this->size(); i++)
      this->buff[i] *= vec.buff[i];
  }

  template <class X>
  template <class Y>
  void Array1D<X>::operator/=(const Array1D<Y>& vec){
    if(this->size() != vec.size())
      throw Array1D_Error("void Array1D<X>::operator/=(const Array1D<Y>& vec): incompatible numbers of elements, " + Subs::str(this->size()) + 
			  " versus " + Subs::str(vec.size()) );
    for(int i=0; i<this->size(); i++)
      this->buff[i] /= vec.buff[i];
  }

  template <class X>
  void Array1D<X>::cos(){
    for(int i=0; i<this->size(); i++) 
      this->buff[i] = std::cos(this->buff[i]);
  }

  template <class X>
  void Array1D<X>::sin(){
    for(int i=0; i<this->size(); i++) 
      this->buff[i] = std::sin(this->buff[i]);
  }

  /** Finds position of x assuming values are ordered. If they are not ordered, the routine will not
   * check, and wrong results will ensue. This one is useful if you have some idea of the likely
   * index value.
   * \param x the value to locate
   * \param jhi input roughly the index expected; returned as the value such that x lies between the jhi-1 and jhi'th
   * elements
   */
  template <class X>
  void Array1D<X>::hunt(const X& x, int& jhi) const {
    unsigned long int uljhi = jhi;
    Subs::hunt(this->buff, this->size(), x, uljhi);
    jhi = int(uljhi);
  }

  /** Finds position of x assuming values are ordered. If they are not ordered, the routine will not
   * check, and wrong results will ensue. This one is useful if you no idea the likely index value.
   * \param x the value to locate
   * \param jhi returned as the value such that x lies between the jhi-1 and jhi'th elements
   */
  template <class X>
  unsigned long Array1D<X>::locate(const X& x) const {
    return Subs::locate(this->buff, this->size(), x);
  }

  /** Sorts an array into ascending order
   */
  template <class X>
  Buffer1D<unsigned long int> Array1D<X>::sort(){
    Buffer1D<unsigned long int> key(this->size());
    heaprank(this->buff, key.ptr(), this->size());
    Array1D<X> temp(*this);
    for(int i=0; i<this->size(); i++)
      this->buff[i] = temp[key[i]];
    return key;
  }
  
  /** Tests whether values increase or decrease monotonically and so whether hunt can be used.
   */
  template <class X>
  bool Array1D<X>::monotonic() const {
    if(this->size() < 3){
      return true;
    }else{
      bool up = this->buff[0] < this->buff[this->size()-1];
      for(int i=1; i<this->size(); i++){
	bool up_now = this->buff[i-1] < this->buff[i];
	if( (up && !up_now) || (!up && up_now)) return false;
      }
      return true;
    }
  }
  
  // Non-member functions

  // Operations with other arrays
  template <class X, class Y>
  Array1D<X> operator-(const Array1D<X>& v1, const Array1D<Y>& v2){
    if(v1.size() != v2.size())
      throw Array1D_Error("void operator-=(const Array1D<X>&, const Array1D<Y>&): incompatible numbers of elements, " + Subs::str(v1.size()) 
			  + " versus " + Subs::str(v2.size()) );
    Array1D<X> temp = v1;
    temp -= v2;
    return temp;
  }

  // Operations with other arrays
  template <class X, class Y>
  Array1D<X> operator+(const Array1D<X>& v1, const Array1D<Y>& v2){
    if(v1.size() != v2.size())
      throw Array1D_Error("void operator+=(const Array1D<X>&, const Array1D<Y>&): incompatible numbers of elements, " + Subs::str(v1.size()) 
			  + " versus " + Subs::str(v2.size()) );
    Array1D<X> temp = v1;
    temp += v2;
    return temp;
  }

  // Plotting routines

  //! Plots an Array1D versus pixel number
  template <class X>
  void pgline(const Array1D<X>& vec){
    if(vec.size()){
      float x[vec.size()];
      float y[vec.size()];
      for(int i=0; i<vec.size(); i++){
	x[i] = float(i+1);
	y[i] = float(vec[i]);
      }
      cpgline(vec.size(),x,y);
    }
  }


  //! Plots two Array1D objects against each other as a line
  template <class X, class Y> 
  void pgline(const Array1D<X>& x, const Array1D<Y>& y){
    if(x.size() > 0){
      if(x.size() == y.size()){
	float xt[x.size()];
	float yt[y.size()];
	for(int i=0; i<x.size(); i++){
	  xt[i] = float(x[i]);
	  yt[i] = float(y[i]);
	}
	cpgline(x.size(),xt,yt);
      }else{
	throw Array1D_Error("void pgline(const Array1D<X>&, const Array1D<Y>): incompatible numbers of elements, " + Subs::str(x.size()) + 
			    " versus " + Subs::str(y.size()) );
      }
    }
  }

  //! Plots two Array1D objects against each other in bar form.
  template <class X, class Y> 
  void pgbin(const Array1D<X>& x, const Array1D<Y>& y){
    if(x.size() > 0){
      if(x.size() == y.size()){
	float xt[x.size()];
	float yt[y.size()];
	for(int i=0; i<x.size(); i++){
	  xt[i] = float(x[i]);
	  yt[i] = float(y[i]);
	}
	cpgbin(x.size(),xt,yt,1);
      }else{
	throw Array1D_Error("void pgbin(const Array1D<X>&, const Array1D<Y>): incompatible numbers of elements, " + Subs::str(x.size()) + 
			    " versus " + Subs::str(y.size()) );
      }
    }
  }

  //! Plots two Array1D objects against each other as a set of points 
  /** 
   * \param x X values
   * \param y Y values
   * \param symbol PGPLOT symbol
   */
  template <class X, class Y> 
  void pgpt(const Array1D<X>& x, const Array1D<Y>& y, int symbol){
    if(x.size() > 0){
      if(x.size() == y.size()){
	float xt[x.size()];
	float yt[y.size()];
	for(int i=0; i<x.size(); i++){
	  xt[i] = x[i];
	  yt[i] = y[i];
	}
	cpgpt(x.size(),xt,yt,symbol);
      }else{
	throw Array1D_Error("void pgpt(const Array1D<X>&, const Array1D<Y>, int): incompatible numbers of elements, " + Subs::str(x.size()) + 
			    " versus " + Subs::str(y.size()) );
      }
    }
  }

  //! Returns the maximum value
  template <class X> 
  X max(const Array1D<X>& vec){
    return vec.max();
  }

  //! Returns the minimum value
  template <class X> 
  X min(const Array1D<X>& vec){
    return vec.min();
  }

  //! Subtracts a constant to create a new Array1D
  template <class X, class Y>
  Array1D<X> operator-(const Array1D<X>& vec, const Y& con){
    Array1D<X> temp = vec;
    temp -= con;
    return temp;
  }

  //! Pre-multiplies by a constant
  template <class X, class Y>
  Array1D<X> operator*(const Y& con, const Array1D<X>& vec){
    Array1D<X> temp = vec;
    temp *= con;
    return temp;
  }

  //! Divides by a constant
  template <class X, class Y>
  Array1D<X> operator/(const Array1D<X>& vec, const Y& con){
    Array1D<X> temp = vec;
    temp /= con;
    return temp;
  }

  //! Mulitplies two Array1Ds, element by element
  template <class X, class Y>
  Array1D<X> operator*(const Array1D<X>& vec1, const Array1D<X>& vec2){
    Array1D<X> temp = vec1;
    temp *= vec2;
    return temp;
  }

  //! Divides two Array1Ds, element by element
  template <class X, class Y>
  Array1D<X> operator/(const Array1D<X>& vec1, const Array1D<X>& vec2){
    Array1D<X> temp = vec1;
    temp /= vec2;
    return temp;
  }

  //! Takes cosine of array
  template <class X>
  Array1D<X> cos(const Array1D<X>& vec){
    Array1D<X> temp = vec;
    temp.cos();
    return temp;
  }

  //! Takes sine of array
  template <class X>
  Array1D<X> sin(const Array1D<X>& vec){
    Array1D<X> temp = vec;
    temp.sin();
    return temp;
  }

  /** Sets an Array1D to a constant
   */
  template <class X> 
  Array1D<X>& Array1D<X>::operator=(const X& con){
    Buffer1D<X>::operator=(con);
    return *this;
  }

};

#endif











