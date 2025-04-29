#ifndef TRM_BUFFER2D
#define TRM_BUFFER2D

#include "trm/subs.h"

namespace Subs {

  //! Error class
  class Buffer2D_Error : public Subs_Error {
  public:
    Buffer2D_Error() : Subs_Error("") {};
    Buffer2D_Error(const std::string& str) : Subs_Error(str) {};
  };

  //! Buffer class for handling 2D memory.
  /** Simple class designed to supply a safe 2D array, 'safe' in
   * that it is deallocated whenever the object goes out of scope.
   * It creates pointers to the array which can then be used
   * in the usual way as a C-style 2D array. The big advantage is that
   * the pointers are automatically deleted when the Buffer2D goes out of scope. 
   * \sa Array2D for handling arithematic on arrays.
   */
  template <class X>
  class Buffer2D {
  public:

    //! Default constructor
    Buffer2D() : nx_(0), ny_(0), buff(NULL) {}

    //! Constructor grabbing space
    Buffer2D(int ny, int nx);

    //! Copy constructor
    Buffer2D(const Buffer2D& obj);

    //! Assignment
    void operator=(const X& x);

    //! Assignment
    Buffer2D& operator=(const Buffer2D& obj);

    //! Destructor
    ~Buffer2D(){
      dealloc();
    }

    //! Resize a buffer (loses data)
    void resize(int ny, int nx);

    //! Returns the number of X pixels (= columns)
    int get_nx() const {return nx_;}

    //! Returns the number of columns (= X pixels)
    int ncol() const {return nx_;}

    //! Returns the number of Y pixels (= rows)
    int get_ny() const {return ny_;}

    //! Returns the number of rows (= Y pixels)
    int nrow() const {return ny_;}

    //! Returns total number of elements
    int size() const {return nx_*ny_;}

    //! Conversion to 2D C-type array
    operator X** const () const {return buff;}

    //! Conversion to 2D C-type array
    operator X**() {return buff;}

    //! Returns 2D C-type array
    template <class Y> 
    void get(Y** ar) const;

    //! Returns 1D C-type array
    template <class Y> 
    void get(Y* ar) const;

    //! Sets data from a 1D C-type array
    template <class Y>
    void set(Y* ar);

    //! Indexed access to the array
    X* operator[](int iy) {return buff[iy];}

    //! Indexed access to the array
    const X* operator[](int iy) const {return buff[iy];}

  protected:

    // The next two functions allow fairly direct access to the
    // underlying pointer to member functions of this and derived
    // classes.
 
    //! Returns pointer to 2D C-type array
    X** const ptr() const {return buff;}

    //! Conversion to 2D C-type array
    X** ptr() {return buff;}

  private:

    // dimensions 
    int nx_, ny_;

    // The pointer (to an array of pointers)
    X** buff;
    
    // sets up memory buffers correctly (for first time)
    void alloc(int ny, int nx);

    // deletes memory buffers
    void dealloc();

  };

  template <class X> 
  void Buffer2D<X>::alloc(int ny, int nx) {
    if(nx == 0 && ny == 0){
      nx_ = ny_ = 0;
      buff = NULL;
    }else{
      nx_ = nx;
      ny_ = ny;
      if(nx < 1 || ny < 1)
	  throw Buffer2D_Error("Subs::Buffer2D<>::alloc(int, int): nx,ny = " + Subs::str(nx) + ", " + Subs::str(ny));
      
      if((buff = new(std::nothrow) X* [ny]) == NULL){
	buff = NULL;
	throw Buffer2D_Error("Subs::Buffer2D<>::alloc(int, int): failed to allocate memory (1)");
      }
      for(int iy=0; iy<ny; iy++){
	if((buff[iy] = new(std::nothrow) X [nx]) == NULL){
	  for(int jy=0; jy<iy; jy++)
	    delete[] buff[jy];
	  delete[] buff;	    
	  buff = NULL;
	  throw Buffer2D_Error("Subs::Buffer2D<>::alloc(int, int): failed to allocate memory (2)");
	}
      }
    }
  }

  template <class X> 
  void Buffer2D<X>::dealloc() {
    if(buff != NULL){
      for(int iy=0; iy<ny_; iy++)
	delete[] buff[iy];
      delete[] buff;
    }
  }

  /** This constructor gets space for exactly ny by nx points. The order
   * or arguments is designed to be the same as for normal C-style arrays which
   * is the reverse of FORTRAN-style. 'ny' is also the number of rows and 'nx'
   * the number of columns if this is being regarded as a matrix.
   * \param ny the number of pixels in Y (rows)
   * \param nx the number of pixels in X (columns)
   */
  template <class X> 
  Buffer2D<X>::Buffer2D(int ny, int nx) : nx_(nx), ny_(ny) {
    try{
      alloc(ny_, nx_);
    }
    catch(const Buffer2D_Error& e){
      throw Buffer2D_Error("Subs::Buffer2D<>::Buffer2D(int, int): " + e);
    }
  }

  /** Copy constructor to make an element by element copy of an object
   */
  template <class X> 
  Buffer2D<X>::Buffer2D(const Buffer2D<X>& obj) : nx_(obj.nx_), ny_(obj.ny_) { 
   
    if(obj.buff == NULL && (nx_ !=0 || ny_ != 0))
      throw Buffer2D_Error("Buffer2D<X>::Buffer2D(const Buffer2D<X>&): null pointer error!");

    // allocate memory
    try{
      alloc(ny_, nx_);
    }
    catch(const Buffer2D_Error& e){
      throw Buffer2D_Error("Subs::Buffer2D<>::Buffer2D(const Buffer2D<>&): " + e);
    }

    // copy over data
    for(int iy=0; iy<ny_; iy++)
      for(int ix=0; ix<nx_; ix++)
	buff[iy][ix] = obj.buff[iy][ix];
  }

  template <class X> 
  Buffer2D<X>& Buffer2D<X>::operator=(const Buffer2D<X>& obj){

    if(this == &obj) return *this;

    // reallocate memory if necessary
    if(nx_ != obj.nx_ || ny_ != obj.ny_){
      dealloc();
      try{
	alloc(obj.ny_, obj.nx_);
      }
      catch(const Buffer2D_Error& e){
	throw Buffer2D_Error("Subs::Buffer2D<>::operator=(const Buffer2D<>&): " + e);
      }
    }
	
    // copy over data
    for(int iy=0; iy<ny_; iy++)
      for(int ix=0; ix<nx_; ix++)
	buff[iy][ix] = obj.buff[iy][ix];

    return *this;
  }

  /** Resizes a Buffer2D. This checks that any change is being made. If it is
   * it deallocates and then reallocates the memory. All existing data will be lost
   * in general.
   * \param ny the number of rows
   * \param nx the number of columns
   */
  template <class X> 
  void Buffer2D<X>::resize(int ny, int nx){
    if(nx_ != nx || ny_ != ny){
      dealloc();
      try{
	alloc(ny, nx);
      }
      catch(const Buffer2D_Error& e){
	throw Buffer2D_Error("Subs::Buffer2D<>::resize(int, int): " + e);
      }
    }
  }

  /** Sets all values in the buffer to a constant
   * \param x the constant value to set
   */
  template <class X>
  void Buffer2D<X>::operator=(const X& x){
    for(int iy=0; iy<ny_; iy++)
      for(int ix=0; ix<nx_; ix++)
	buff[iy][ix] = x;
  }

  /** Copies all the array values into an ordinary C style 2D array
   * the memory for which which must have been allocated beforehand.
   */
  template <class X> 
  template <class Y>
  void Buffer2D<X>::get(Y** arr) const {
    for(int iy=0; iy<ny_; iy++)
      for(int ix=0; ix<nx_; ix++)
	arr[iy][ix] = Y(buff[iy][ix]);
  }

  /** Copies all the array values into an ordinary C style 1D array
   * the memory for which which must have been allocated beforehand.
   */
  template <class X> 
  template <class Y>
  void Buffer2D<X>::get(Y* ar) const {
    for(int iy=0, k=0; iy<ny_; iy++)
      for(int ix=0; ix<nx_; ix++, k++)
	ar[k] = Y(buff[iy][ix]);
  }

  /** Copies all the values from an ordinary C style 1D array.
   */
  template <class X> 
  template <class Y>
  void Buffer2D<X>::set(Y* ar) {
    for(int iy=0, k=0; iy<ny_; iy++)
      for(int ix=0; ix<nx_; ix++, k++)
	buff[iy][ix] = X(ar[k]);
  }

  //! ASCII output
  template <class X>
  std::ostream& operator<<(std::ostream& s, const Buffer2D<X>& obj){
    for(int iy=0; iy<obj.get_ny(); iy++)
      for(int ix=0; ix<obj.get_nx(); ix++)
	s << ix+1 << " " << iy+1 << " " << obj[iy][ix] << "\n";
    return s;
  }

}

#endif
