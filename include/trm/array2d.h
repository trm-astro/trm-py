#ifndef TRM_SUBS_ARRAY2D
#define TRM_SUBS_ARRAY2D

#include <cmath>
#include <new>
#include <string>
#include <iostream>
#include "cpgplot.h"
#include "trm/buffer2d.h"

namespace Subs {
  
  //! Error class for Array2D
  class Array2D_Error : public Subs_Error {
  public:
    Array2D_Error() : Subs_Error("") {};
    Array2D_Error(const std::string& str) : Subs_Error(str) {};
  };

  // Pre-declaration of friends
  template <class X>
  class Array2D;

  template <class X>
  Array2D<X> operator-(const X& con, const Array2D<X>& img);
    
  template <class X>
  std::ostream& operator<<(std::ostream& ost, const Array2D<X>& img);
  
  //! Array2D, a 2D array class
  
  /** This class is for representing numerical 2D arrays and matrices. Use it if
   * you want to be able to add, subtract etc with arrays, otherwise use Buffer2D
   * from which it is inherited.
   */
  
  template <class X> 
  class Array2D : public Buffer2D<X> {
    
  public:
    
    //! Default constructor
    Array2D() : Buffer2D<X>() {}

    //! Constructs an ny rows by nx columns matrix
    Array2D(int ny, int nx) : Buffer2D<X>(ny,nx) {}

    //! Assignment to a constant
    Array2D& operator=(const X& con);

    //! Assignment to another Array2D
    Array2D& operator=(const Array2D<X>& obj);

    // functions

    //! Returns minimum value
    X min() const;

    //! Returns maximum value
    X max() const;

    //! Returns mean value
    X mean() const;

    //! Returns RMS value
    X rms() const;

    //! Returns median value
    X median() const;

    //! Returns a percentiles 
    X centile(float frac) const;

    //! Returns two percentiles 
    void centile(float frac1, float frac2, X& t1, X& t2) const;

    //! Returns sum over all elements
    X sum() const;

    //! Takes square root of all elements
    void sqrt(); 

    //! Transform using a step function at a certain threshold
    void step(const X& thresh);
    
    //! Addition of a constant in place
    void operator+=(const X& con);

    //! Subtraction of a constant in place
    void operator-=(const X& con);

    //! Multiplication by a constant in place
    void operator*=(const X& con);

    //! Division by a constant in place
    void operator/=(const X& con);
    
    //! Addition of an Array2D in place
    void operator+=(const Array2D<X>& img);

    //! Subtraction of an Array2D in place
    void operator-=(const Array2D<X>& img);

    //! Multiplication by an Array2D in place
    void operator*=(const Array2D<X>& img);

    //! Division by an Array2D in place
    void operator/=(const Array2D<X>& img);
    
    //! Binary output
    void write(std::ostream& ostr) const;

    //! Binary input
    void read(std::istream& istr);
    
    //! Array2D = constant - Array2D
    friend Array2D<X> operator-<>(const X& con, const Array2D<X>& img);
    
    //! ASCII output
    friend std::ostream& operator<<<>(std::ostream& ost, const Array2D<X>& img);
    
  };

  template <class X>
  Array2D<X>& Array2D<X>::operator=(const X& con){
    this->Buffer2D<X>::operator=(con);
    return *this;
  }

  template <class X>
  Array2D<X>& Array2D<X>::operator=(const Array2D<X>& obj){
    this->Buffer2D<X>::operator=(obj);
    return *this;
  }

  template <class X>
  X Array2D<X>::max() const {
    if(this->get_ny()){
      X t = this->ptr()[0][0];
      for(int iy=0; iy<this->get_ny(); iy++)
	for(int ix=0; ix<this->get_nx(); ix++)
	  if(t < this->ptr()[iy][ix]) t = this->ptr()[iy][ix];
      return t;
    }else{
      throw Array2D_Error("X Array2D<X>::max(): operation undefined on null array");
    }
  }

  template <class X>
  X Array2D<X>::min() const {
    if(this->get_ny()){
      X t = this->ptr()[0][0];
      for(int iy=0; iy<this->get_ny(); iy++)
	for(int ix=0; ix<this->get_nx(); ix++)
	  if(t > this->ptr()[iy][ix]) t = this->ptr()[iy][ix];
      return t;
    }else{
      throw Array2D_Error("X Array2D<X>::min(): operation undefined on null array");
    }
  }

  /**
   * The mean is computed in double precision to maintain accuracy over large
   * numbers of pixels.
   */
  template <class X>
  X Array2D<X>::mean() const {
    double t = 0;
    for(int iy=0; iy<this->get_ny(); iy++)
      for(int ix=0; ix<this->get_nx(); ix++)
	t += double(this->ptr()[iy][ix]);
    return X(t/this->size());
  }

  template <class X>
  X Array2D<X>::rms() const {
    double t = 0, m = mean();
    for(int iy=0; iy<this->get_ny(); iy++)
      for(int ix=0; ix<this->get_nx(); ix++)
	t += Subs::sqr(this->ptr()[iy][ix]-m);
    return std::sqrt(t/(std::max(this->size()-1,int(1))));
  }

  template <class X>
  X Array2D<X>::centile(float frac) const {
    if(this->get_ny()){

      unsigned long int nelem = this->size();

      X *ptr = new X[nelem];
      this->get(ptr);

      X t;
      if(frac <= 0.f){
	t = select(ptr, nelem, 0);
      }else if(frac >= 1.f){
	t = select(ptr, nelem, nelem-1);
      }else{
	double x = (nelem-1)*frac;
	unsigned long int k = (unsigned long int)floor(x);
	if(k == nelem-1){
	  t = select(ptr, nelem, k);
	}else{
	  t = (k+1-x)*select(ptr, nelem, k) + (x-k)*select(ptr, nelem, k+1);
	}
      }

      delete[] ptr;
      return t;

    }else{
      throw Array2D_Error("X Array2D<X>::centile(float): operation undefined on null array");
    }
  }

  template <class X>
  void Array2D<X>::centile(float frac1, float frac2, X& t1, X& t2) const {
    if(this->nrow()){

      unsigned long int nelem = this->size();

      X *ptr = new X[nelem];
      this->get(ptr);

      if(frac1 <= 0.f){
	t1 = select(ptr, nelem, 0);
      }else if(frac1 >= 1.f){
	t1 = select(ptr, nelem, nelem-1);
      }else{
	double x1 = (nelem-1)*frac1;
	unsigned long int k1 = (unsigned long int)floor(x1);
	if(k1 == nelem-1){
	  t1 = select(ptr, nelem, k1);
	}else{
	  t1 = (k1+1-x1)*select(ptr, nelem, k1) + (x1-k1)*select(ptr, nelem, k1+1);
	}
      }

      if(frac2 <= 0.f){
	t2 = select(ptr, nelem, 0);
      }else if(frac2 >= 1.f){
	t2 = select(ptr, nelem, nelem-1);
      }else{
	double x2 = (nelem-1)*frac2;
	unsigned long int k2 = (unsigned long int)floor(x2);
	if(k2 == nelem-1){
	  t2 = select(ptr, nelem, k2);
	}else{
	  t2 = (k2+1-x2)*select(ptr, nelem, k2) + (x2-k2)*select(ptr, nelem, k2+1);
	}
      }

      delete[] ptr;

    }else{
      throw Array2D_Error("X Array2D<X>::centile(float, float, X&, X&): operation undefined on null array");
    }
  }

  template <class X>
  X Array2D<X>::median() const {
    if(this->get_ny()){
      unsigned long int nelem = this->size();

      X *ptr = new X[nelem];
      X med;
      this->get(ptr);
      if(nelem % 2 == 0){
	med = (Subs::select(ptr, nelem, nelem/2-1) + Subs::select(ptr, nelem, nelem/2))/2;
      }else{
	med = Subs::select(ptr, nelem, nelem/2);
      }
      delete[] ptr;
      return med;
    }else{
      throw Array2D_Error("X Array2D<X>::median(): operation undefined on null array");
    }
  }

  template <class X>
  X Array2D<X>::sum() const {
    if(this->get_ny()){
      X t = 0;
      for(int iy=0; iy<this->get_ny(); iy++)
	for(int ix=0; ix<this->get_nx(); ix++)
	  t += this->ptr()[iy][ix];
      return t;
    }else{
      throw Array2D_Error("X Array2D<X>::sum(): operation undefined on null array");
    }
  }

  template <class X>
  X max(const Array2D<X>& img){
    return img.max();
  }

  template <class X>
  X min(const Array2D<X>& img){
    return img.min();
  }

  template <class X>
  X sum(const Array2D<X>& img){
    return img.sum();
  }

    template <class X>
    void Array2D<X>::sqrt(){
	for(int iy=0; iy<this->get_ny(); iy++)
	    for(int ix=0; ix<this->get_nx(); ix++)
		this->ptr()[iy][ix] = std::sqrt(this->ptr()[iy][ix]);
    }

    /** Converts all elements of the array to 0 or 1 depending on whether they are
     * above or below the specified threshold.
     * \param thresh the threshold level. Elements greater than this will be 1,
     * elements less than or equal to it will become 0.
     */
    template <class X>
    void Array2D<X>::step(const X& thresh) {
	for(int iy=0; iy<this->get_ny(); iy++)
	    for(int ix=0; ix<this->get_nx(); ix++)
		this->ptr()[iy][ix] = this->ptr()[iy][ix] > thresh ? X(1) : X(0);
    }

  // Arithematic

  template <class X>
  void Array2D<X>::operator+=(const X& con){
    for(int iy=0; iy<this->get_ny(); iy++)
      for(int ix=0; ix<this->get_nx(); ix++)
	this->ptr()[iy][ix] += con;
  }

  template <class X>
  void Array2D<X>::operator-=(const X& con){
    for(int iy=0; iy<this->get_ny(); iy++)
      for(int ix=0; ix<this->get_nx(); ix++)
	this->ptr()[iy][ix] -= con;
  }

  template <class X>
  void Array2D<X>::operator*=(const X& con){
    for(int iy=0; iy<this->get_ny(); iy++)
      for(int ix=0; ix<this->get_nx(); ix++)
	this->ptr()[iy][ix] *= con;
  }

  template <class X>
  void Array2D<X>::operator/=(const X& con){
    for(int iy=0; iy<this->get_ny(); iy++)
      for(int ix=0; ix<this->get_nx(); ix++)
	this->ptr()[iy][ix] /= con;
  }

  template <class X>
  void Array2D<X>::operator+=(const Array2D<X>& img){
    if(this->get_ny() != img.get_ny() || this->get_nx() != img.get_nx()){
      throw Array2D_Error("void Array2D<X>::operator+=(Array2D<X>&): mis-matching array dimensions; ["
			  + Subs::str(this->get_ny()) + "][" + Subs::str(this->get_nx()) + "] vs ["
			  + Subs::str(img.get_ny()) + "][" + Subs::str(img.get_nx()) + "]");
    }else{
      for(int iy=0; iy<this->get_ny(); iy++)
	for(int ix=0; ix<this->get_nx(); ix++)
	  this->ptr()[iy][ix] += img.ptr()[iy][ix];
    }
  }

  template <class X>
  void Array2D<X>::operator-=(const Array2D<X>& img){
    if(this->get_ny() != img.get_ny() || this->get_nx() != img.get_nx()){
      throw Array2D_Error("void Array2D<X>::operator-=(Array2D<X>&): mis-matching array dimensions; ["
			  + Subs::str(this->get_ny()) + "][" + Subs::str(this->get_nx()) + "] vs ["
			  + Subs::str(img.get_ny()) + "][" + Subs::str(img.get_nx()) + "]");
    }else{
      for(int iy=0; iy<this->get_ny(); iy++)
	for(int ix=0; ix<this->get_nx(); ix++)
	  this->ptr()[iy][ix] -= img.ptr()[iy][ix];
    }
  }

  template <class X>
  void Array2D<X>::operator*=(const Array2D<X>& img){
    if(this->get_ny() != img.get_ny() || this->get_nx() != img.get_nx()){
      throw Array2D_Error("void Array2D<X>::operator*=(Array2D<X>&): mis-matching array dimensions; ["
			  + Subs::str(this->get_ny()) + "][" + Subs::str(this->get_nx()) + "] vs ["
			  + Subs::str(img.get_ny()) + "][" + Subs::str(img.get_nx()) + "]");
    }else{
      for(int iy=0; iy<this->get_ny(); iy++)
	for(int ix=0; ix<this->get_nx(); ix++)
	  this->ptr()[iy][ix] *= img.ptr()[iy][ix];
    }
  }

  template <class X>
  void Array2D<X>::operator/=(const Array2D<X>& img){
    if(this->get_ny() != img.get_ny() || this->get_nx() != img.get_nx()){
      throw Array2D_Error("void Array2D<X>::operator/=(Array2D<X>&): mis-matching array dimensions; ["
			  + Subs::str(this->get_ny()) + "][" + Subs::str(this->get_nx()) + "] vs ["
			  + Subs::str(img.get_ny()) + "][" + Subs::str(img.get_nx()) + "]");
    }else{
      for(int iy=0; iy<this->get_ny(); iy++)
	for(int ix=0; ix<this->get_nx(); ix++)
	  this->ptr()[iy][ix] /= img.ptr()[iy][ix];
    }
  }

  // non-member arithematic functions
  template <class X>
  Array2D<X> operator+(const Array2D<X>& img1, const Array2D<X>& img2){
    Array2D<X> temp = img1;
    temp += img2;
    return temp;
  }

  template <class X>
  Array2D<X> operator+(const Array2D<X>& img, const X& con){
    Array2D<X> temp = img;
    temp += con;
    return temp;
  }

  template <class X>
  Array2D<X> operator+(const X& con, const Array2D<X>& img){
    Array2D<X> temp = img;
    temp += con;
    return temp;
  }

  template <class X>
  Array2D<X> operator-(const Array2D<X>& img1, const Array2D<X>& img2){
    Array2D<X> temp = img1;
    temp -= img2;
    return temp;
  }

  template <class X>
  Array2D<X> operator-(const Array2D<X>& img, const X& con){
    Array2D<X> temp = img;
    temp -= con;
    return temp;
  }

  template <class X>
  Array2D<X> operator-(const X& con, const Array2D<X>& img){
    Array2D<X> temp(img.get_ny(),img.get_nx());
    if(img.get_ny()){
      for(int iy=0; iy<img.get_ny(); iy++)
	for(int ix=0; ix<img.get_nx(); ix++)
	  temp.ptr()[iy][ix] = con - img.ptr()[iy][ix];
    }else{
      throw Array2D_Error("X Array2D<X>::operator-(X& con, Array2D<X>): operation undefined on null array");
    }
    return temp;
  }

  template <class X>
  Array2D<X> operator*(const Array2D<X>& img, const X& con){
    Array2D<X> temp = img;
    temp *= con;
    return temp;
  }

  template <class X>
  Array2D<X> operator*(const X& con, const Array2D<X>& img){
    Array2D<X> temp = img;
    temp *= con;
    return temp;
  }

  template <class X>
  Array2D<X> operator*(const Array2D<X>& img1, const Array2D<X>& img2){
    Array2D<X> temp = img1;
    temp *= img2;
    return temp;
  }

  template <class X>
  Array2D<X> operator/(const Array2D<X>& img, const X& con){
    Array2D<X> temp = img;
    temp /= con;
    return temp;
  }

  template <class X>
  Array2D<X> operator/(const Array2D<X>& img1, const Array2D<X>& img2){
    Array2D<X> temp = img1;
    temp /= img2;
    return temp;
  }

  template <class X>
  std::ostream& operator<<(std::ostream& ost, const Array2D<X>& img){
    if(img.get_ny()){
      if(img.get_ny() < 5){
	for(int iy=0; iy<img.get_ny(); iy++){
	  ost << "  ";
	  for(int ix=0; ix<img.get_nx(); ix++)
	    ost <<  img.ptr()[iy][ix] << " ";
	  ost << std::endl;
	}
      }else{
	for(int iy=0; iy<img.get_ny(); iy++)
	  for(int ix=0; ix<img.get_nx(); ix++)
	    ost << "[" << iy << "," << ix << "] = " << img.ptr()[iy][ix] << std::endl;
      }
      return ost;
    }else{
      throw Array2D_Error("std::ostream operator<<(std::ostream&, Array2D<X>): operation undefined on null array");
    }
  }

  template <class X>
  void pggray(const Array2D<X>& img, float a1, float a2, float* tr){
    if(img.get_ny()){
      float *ptr = new float[img.size()];
      if(ptr){
	img.get(ptr);
	cpggray(ptr,img.get_nx(),img.get_ny(),1,img.get_nx(),1,img.get_ny(),a1,a2,tr);
	delete [] ptr;
      }else{
	throw Array2D_Error("pggray(const Array2D<X>&, float, float, float*): failed to allocate buffer");
      } 
    }
  }

  template <class X>
  void pgimag(const Array2D<X>& img, float a1, float a2, float tr[]){
    if(img.get_ny()){
      float a[img.size()];
      img.get(a);
      cpgimag(a,img.get_nx(),img.get_ny(),1,img.get_nx(),1,img.get_ny(),a1,a2,tr);
    }
  }

  template <class X>
  void pgcont(const Array2D<X>& img, float c[], int nc, float tr[]){
    if(img.get_ny()){
      float a[img.size()];
      img.get(a);
      cpgcont(a,img.get_nx(),img.get_ny(),1,img.get_nx(),1,img.get_ny(),c,nc,tr);
    }
  }

  template <class X>
  void Array2D<X>::write(std::ostream& ostr) const {
    int ny = this->get_ny(), nx = this->get_nx();
    ostr.write((char*)&ny,sizeof(ny));
    ostr.write((char*)&nx,sizeof(nx));
    for(int iy=0; iy<ny; iy++)
      ostr.write((char*)this->ptr()[iy],sizeof(X[nx]));
  }

  template <class X>
  void Array2D<X>::read(std::istream& istr) {
    int ny, nx;
    istr.read((char*)&ny,sizeof(ny));
    istr.read((char*)&nx,sizeof(nx));
    this->resize(ny,nx);
    for(int iy=0; iy<ny; iy++)
      istr.read((char*)this->ptr()[iy],sizeof(X[nx]));
  }

};

#endif











