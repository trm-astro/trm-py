#ifndef TRM_HITEM
#define TRM_HITEM
 
#include <fstream>
#include <string>
#include <vector>
#include "trm/subs.h"
#include "trm/date.h"
#include "trm/time.h"
#include "trm/position.h"
#include "trm/telescope.h"

namespace Subs {

    //! Base class for storage of heterogeneous data types
    
    /**
     * Hitem is an abstract base class for the storage of
     * heterogeneous data types such as int, float etc in a consistent
     * manner with a comment field. These can be stored into homogeneous
     * containers using pointers to Hitems, i.e. Hitem*.
     */
    
    class Hitem{
	
    public:
	
	// The numbers below are fixed.
	//! Used for reading/writing to disk. 
	enum HITEM_TYPE {
	    HDOUBLE    = 0, 
	    HCHAR      = 1, 
	    HINT       = 2, 
	    HUINT      = 3, 
	    HLINT      = 4, 
	    HULINT     = 5, 
	    HFLOAT     = 6, 
	    HSTRING    = 7, 
	    HBOOL      = 8, 
	    HDIRECTORY = 9, 
	    HDATE      = 10, 
	    HTIME      = 11, 
	    HPOSITION  = 12, 
	    HDVECTOR   = 13, 
	    HUCHAR     = 14, 
	    HTELESCOPE = 15,
	    HUSINT     = 16,
	    HIVECTOR   = 17, 
	    HFVECTOR   = 18, 
	};
	
	//! Printing mode
	static int pmode;
	
	//! field width for value output    
	static size_t value_width; 

	//! field width for type output    
	static size_t type_width; 
	
	//! precision for output of doubles 
	static size_t double_precision; 
	
	//! precision for output of floats
	static size_t float_precision; 
	
	//! Default constructor
	Hitem() : comment("") {};
	
	//! Constructor from a comment string.
	Hitem(const std::string& com) : comment(com) {};
	
	//! Destructor
	virtual ~Hitem(){};
	
	//! Returns the comment field
	const std::string& get_comment() const {return comment;};
	
	//! Set the comment value
	void set_comment(const std::string& com){
		// Check on if(this) removed as it is removed at compile time
		comment = com;
	};
	
	//! Name of type to use in output
	virtual std::string type() const = 0;
	
	//! Code for use in I/O
	virtual HITEM_TYPE type_code() const = 0;
	
	//! Is the object a directory?
	virtual bool is_a_dir() const = 0;
	
	//! Get value as a char
	virtual void get_value(char &cval) const;
	
	//! Get value as a double
	virtual void get_value(double &dval) const;
	
	//! Get value as an int
	virtual void get_value(int &ival) const;
	
	//! Get value as an unsigned int
	virtual void get_value(unsigned int &ival) const;
	
	//! Get value as a long int
	virtual void get_value(long int &ival) const;
	
	//! Get value as an unsigned long int
	virtual void get_value(unsigned long int &ival) const;
	
	//! Get value as a float int
	virtual void get_value(float  &fval) const;
	
	//! Get value as a boolean
	virtual void get_value(bool &bval) const;
	
	//! Get value as a string
	virtual void get_value(std::string &sval) const = 0;
	
	//! Get value as a Date
	virtual void get_value(Date   &dval) const;
	
	//! Get value as a Time
	virtual void get_value(Time   &tval) const;
	
	//! Get value as a Position
	virtual void get_value(Position &pval) const;
	
	//! Get value as a Telescope
	virtual void get_value(Telescope &pval) const;
	
	//! Get value as a vector of doubles
	virtual void get_value(std::vector<double> &val) const;
	
	//! Get value as an unsigned char
	virtual void get_value(unsigned char& val) const;
	
	//! Get value as a char
	virtual char get_char() const;
	
	//! Get value as an unsigned char
	virtual unsigned char get_uchar() const;
	
	//! Get value as a double
	virtual double get_double() const;
	
	//! Get value as an int
	virtual int get_int() const;
	
	//! Get value as an unsigned int
	virtual unsigned int get_uint() const;
	
	//! Get value as a long int
	virtual long int get_lint() const;
	
	//! Get value as an unsigned long int
	virtual unsigned long int get_ulint() const;
	
	//! Get value as a float
	virtual float  get_float() const;
	
	//! Get value as a boolean
	virtual bool   get_bool() const;
	
	//! Get value as a string
	virtual std::string get_string() const = 0;
	
	//! Get value as a Date
	virtual const Date& get_date() const;
	
	//! Get value as a Time
	virtual const Time& get_time() const;
	
	//! Get value as a Position
	virtual const Position& get_position() const;
	
	//! Get value as a Telescope
	virtual const Telescope& get_telescope() const;
	
	//! Get value as a vector of doubles
	virtual const std::vector<double>& get_dvector() const;

	//! Get value as a vector of integers
	virtual const std::vector<int>& get_ivector() const;

	//! Get value as a vector of floats
	virtual const std::vector<float>& get_fvector() const;
	
	//! Set value from a char
	virtual void set_value(const char &cval);
	
	//! Set value from an unsigned char
	virtual void set_value(const unsigned char &cval);
	
	//! Set value from a double
	virtual void set_value(const double &dval);
	
	//! Set value from an int
	virtual void set_value(const int &ival);
	
	//! Set value from an unsigned int
	virtual void set_value(const unsigned int &ival);
	
	//! Set value from a long int
	virtual void set_value(const long int &ival);
	
	//! Set value from an unsigned long int
	virtual void set_value(const unsigned long int &ival);
	
	//! Set value from a float
	virtual void set_value(const float  &fval);
	
	//! Set value from a boolean
	virtual void set_value(const bool  &bval);
	
	//! Set value from a string
	virtual void set_value(const std::string &sval) = 0;
	
	//! Set value from a date
	virtual void set_value(const Date   &dval);
	
	//! Set value from a Time
	virtual void set_value(const Time   &tval);
	
	//! Set value from a Position
	virtual void set_value(const Position &pval);
	
	//! Set value from a Telescope
	virtual void set_value(const Telescope& tval);
	
	//! Set value from a vector of doubles
	virtual void set_value(const std::vector<double> &val);
	
	//! Print the value
	virtual void print(std::ostream& s) const = 0;
	
	//! Binary output
	virtual void write(std::ofstream& ostr) const = 0;
	
	//! Binary input
	virtual void read(std::ifstream& istr, bool swap_bytes) = 0;
	
	//! Binary output of the type
	virtual void write_type(std::ofstream& ostr) const;
	
	//! Generate a copy of the item
	virtual Hitem* copy() const  = 0;
	
	//! Write an Hitem to an ASCII file
	virtual void write_ascii(std::ofstream& ostr) const = 0;
	
	//! Binary input of the type
	static HITEM_TYPE read_type(std::ifstream& istr, bool swap_bytes);

	//! Read any Hitem from a disk file
	static Hitem* read_item(std::ifstream& istr, bool swap_bytes);

	//! Read any Hitem from an ASCII file
	static Hitem* read_ascii_item(std::ifstream& istr);

	//! Skip any Hitem in a disk file
	static void skip_item(std::ifstream& istr, bool swap_byes);

	//! Writes an Hitem to disk
	static void write_item(std::ofstream& ostr, Hitem* item);

	//! Set output format
	static void set_default_output();
    
	//! Error class
	class Hitem_Error : public Subs_Error { 
	public:
	    //! Default constructor
	    Hitem_Error () : Subs_Error() {}

	    //! Constructor from a string
	    Hitem_Error (const std::string& str) : Subs_Error(str) {}
	};
    
    protected:
    
	std::string comment;
    
    };
  
  
    //! Storage of chars

    class Hchar : public Hitem {
    
    public:

	//! Default constructor
	Hchar() : Hitem(), value(' ') {} ;

	//! Constructor from a value plus a comment
	Hchar(char cval, const std::string& com="") : Hitem(com), value(cval) {};

	//! Constructor from file input
	Hchar(std::ifstream& istr, bool swap_bytes){read(istr, swap_bytes);}

	//! Destructor
	~Hchar(){};  
    
	//! Returns name to use when referring to this type
	std::string type() const {return "char";};

	//! Type code for use in disk I/O
	HITEM_TYPE type_code() const {return HCHAR;}

	//! For tests of whether an Hitem is a directory or not
	bool is_a_dir() const {return false;};

	//! Gets value as a char
	void  get_value(char &cval) const {cval = value;};

	//! Gets value as a char
	char  get_char() const {return value;};

	//! Sets value    
	void  set_value(const char &cval){value = cval;};

	void get_value(std::string& sval) const;

	std::string get_string() const;

	void set_value(const std::string& sval);

	//! Prints value
	void print(std::ostream& s) const;

	//! Binary output
	void   write(std::ofstream& ostr) const;

	//! ASCII output
	void   write_ascii(std::ofstream& ostr) const;

	//! Binary input
	void   read(std::ifstream& istr, bool swap_bytes);

	//! Skip in a binary file
	static void skip(std::ifstream& istr, bool swap_bytes);

	//! Copy an character header item and return a pointer to it.
	Hitem* copy() const {return new Hchar(value,comment);}

    private:
    
	CHAR value;
    
    };

    //! Storage of doubles
  
    class Hdouble : public Hitem {
    
    public:
    
	//! Default constructor
	Hdouble() : Hitem(), value(0.) {} ;

	//! Constructor from a value plus a comment
	Hdouble(double dval, const std::string& com="") : Hitem(com), value(dval) {};

	//! Constructor from file input
	Hdouble(std::ifstream& istr, bool swap_bytes){read(istr, swap_bytes);}

	//! Destructor
	~Hdouble(){};  
        
	//! Returns name to use when referring to this type
	std::string type() const {return "double";};

	//! Type code for use in disk I/O
	HITEM_TYPE type_code() const {return HDOUBLE;}

	//! For tests of whether an Hitem is a directory or not
	bool is_a_dir() const  {return false;};

	//! Gets a value as a double    
	void   get_value(double &dval) const {dval = value;};

	//! Gets a value as a float
	void   get_value(float &fval) const  {fval = float(value);};

	//! Gets a value as a double    
	double get_double() const {return value;};

	//! Gets a value as a float
	float  get_float() const {return float(value);};

	//! Sets the value from an int
	void   set_value(const int &ival){value = double(ival);};

	//! Sets the value from an unsigned int
	void   set_value(const unsigned int  &ival){value = double(ival);};

	//! Sets the value from a long int
	void   set_value(const long int &ival){value = double(ival);};

	//! Sets the value from an unsigned long int
	void   set_value(const unsigned long int &ival){value = double(ival);};

	//! Sets the value from a double
	void   set_value(const double &dval){value = dval;};

	//! Sets the value from a float
	void   set_value(const float  &fval){value = double(fval);};

	void get_value(std::string& sval) const;

	std::string get_string() const;

	void set_value(const std::string& sval);
    
	//! ASCII print
	void   print(std::ostream& s) const;

	//! Binary input
	void   write(std::ofstream& ostr) const;

	//! ASCII output
	void   write_ascii(std::ofstream& ostr) const;

	//! Binary output
	void   read(std::ifstream& istr, bool swap_bytes);

	//! Skip in a binary disk file
	static void skip(std::ifstream& istr, bool swap_bytes);

	//! Copy an character header item and return a pointer to it.
	Hitem* copy() const {return new Hdouble(value,comment);}

    private:
    
	REAL8 value;
    
    };
  
    //! Storage of ints
  
    class Hint : public Hitem {
    
    public:
    
	//! Default constructor
	Hint() : Hitem(), value(0) {};

	//! Constructor from a value plus a comment
	Hint(int ival, const std::string& com="") : Hitem(com), value(ival) {};

	//! Constructor from file input
	Hint(std::ifstream& istr, bool swap_bytes){read(istr, swap_bytes);}

	//! Destructor
	~Hint(){};
    
	//! Gets value as an int
	void   get_value(int &ival) const  {ival = value;};
    
	//! Gets value as an int
	void   get_value(double &dval) const {dval = double(value);};
    
	//! Gets value as an int
	void   get_value(float &fval) const {fval = float(value);};
    
	//! Gets value as an int
	int    get_int() const {return value;};
    
	//! Gets value as a double
	double get_double() const {return double(value);};
    
	//! Gets value as a float
	float  get_float() const {return float(value);};

	void get_value(std::string& sval) const;

	std::string get_string() const;

	void set_value(const std::string& sval);

	//! Returns name to use when referring to this type
	std::string type() const {return "int";}

	//! Type code for use in disk I/O
	HITEM_TYPE type_code() const {return HINT;}

	//! For tests of whether an Hitem is a directory or not
	bool is_a_dir() const  {return false;}
    
	void   set_value(const int &ival){value = ival;};
    
	void   print(std::ostream& s) const;

	void   write(std::ofstream& ostr) const;

	//! ASCII output
	void   write_ascii(std::ofstream& ostr) const;

	void   read(std::ifstream& istr, bool swap_bytes);

	static void skip(std::ifstream& istr, bool swap_bytes);

	//! Copy an integer header item and return a pointer to it.
	Hitem* copy() const {return new Hint(value,comment);}
    
    private:
    
	INT4 value;
    
    };
  
    //! Storage of unsigned ints
  
    class Huint : public Hitem {
    
    public:
    
	//! Default constructor
	Huint() : Hitem(), value(0) {};

	//! Constructor from a value plus a comment
	Huint(unsigned int ival, const std::string& com="") : Hitem(com), value(ival) {};

	//! Constructor from file input
	Huint(std::ifstream& istr, bool swap_bytes){read(istr, swap_bytes);}

	//! Destructor
	~Huint(){};
    
	void   get_value(unsigned int &ival) const  {ival = value;};
	void   get_value(double &dval) const {dval = double(value);};
	void   get_value(float &fval) const {fval = float(value);};
	unsigned int get_uint() const {return value;};
	double get_double() const {return double(value);};
	float  get_float() const {return float(value);};

	void get_value(std::string& sval) const;
	std::string get_string() const;
	void set_value(const std::string& sval);
    
	//! Returns name to use when referring to this type
	std::string type() const {return "unsigned int";}

	//! Type code for use in disk I/O
	HITEM_TYPE type_code() const {return HUINT;}

	//! For tests of whether an Hitem is a directory or not
	bool is_a_dir() const  {return false;}
    
	void   set_value(const unsigned int &ival){value = ival;};
    
	void   print(std::ostream& s) const;

	void   write(std::ofstream& ostr) const;

	//! ASCII output
	void   write_ascii(std::ofstream& ostr) const;

	void   read(std::ifstream& istr, bool swap_bytes);

	static void skip(std::ifstream& istr, bool swap_bytes);

	//! Copy an unsigned integer header item and return a pointer to it.
	Hitem* copy() const {return new Huint(value,comment);}

    
    private:
    
	UINT4 value;
    
    };

    //! Storage of floats
  
    class Hfloat : public Hitem {
    
    public:
    
	//! Default constructor
	Hfloat() : Hitem(), value(0.) {}; 

	//! Constructor from a value plus a comment
	Hfloat(float fval, const std::string& com="") : Hitem(com), value(fval) {};

	//! Constructor from file input
	Hfloat(std::ifstream& istr, bool swap_bytes){read(istr, swap_bytes);}

	//! Destructor
	~Hfloat(){};
	void   get_value(double &dval) const {dval = double(value);};
	void   get_value(float &fval) const {fval  = value;};
	double get_double() const {return double(value);};
	float  get_float() const {return value;};

	void get_value(std::string& sval) const;
	std::string get_string() const;
	void set_value(const std::string& sval);

	//! Returns name to use when referring to this type
	std::string type() const {return "float";};

	//! Type code for use in disk I/O
	HITEM_TYPE type_code() const {return HFLOAT;}

	//! For tests of whether an Hitem is a directory or not
	bool is_a_dir() const  {return false;}
    
	void   set_value(const int &ival){value = float(ival);};
	void   set_value(const unsigned int &ival){value = float(ival);};
	void   set_value(const long int &ival){value = float(ival);};
	void   set_value(const unsigned long int &ival){value = float(ival);};
	void   set_value(const double &dval){value = float(dval);};
	void   set_value(const float &fval){value = fval;};
    
	void   print(std::ostream& s) const;

	void   write(std::ofstream& ostr) const;

	//! ASCII output
	void   write_ascii(std::ofstream& ostr) const;

	void   read(std::ifstream& istr, bool swap_bytes);

	static void skip(std::ifstream& istr, bool swap_bytes);
    
	//! Copy a float integer header item and return a pointer to it.
	Hitem* copy() const {return new Hfloat(value,comment);}

    private:
    
	REAL4 value;
    
    };

    //! Storage of bools
  
    class Hbool : public Hitem {
    
    public:
    
	//! Default constructor
	Hbool() : Hitem(), value(0) {}; 

	//! Constructor from a value plus a comment
	Hbool(bool fval, const std::string& com="") : Hitem(com), value(fval) {};

	//! Constructor from file input
	Hbool(std::ifstream& istr, bool swap_bytes){read(istr, swap_bytes);}

	//! Destructor
	~Hbool(){};

	void   get_value(bool &bval) const {bval  = value;};

	bool  get_bool() const {return value;};

	void get_value(std::string& sval) const;

	std::string get_string() const;

	void set_value(const std::string& sval);
    
	//! Returns name to use when referring to this type
	std::string type() const {return "bool";};

	//! Type code for use in disk I/O
	HITEM_TYPE type_code() const {return HBOOL;}

	//! For tests of whether an Hitem is a directory or not
	bool is_a_dir() const  {return false;}
    
	void   set_value(const bool &bval){value = bval;};
    
	void   print(std::ostream& s) const;

	void   write(std::ofstream& ostr) const;

	//! ASCII output
	void   write_ascii(std::ofstream& ostr) const;

	void   read(std::ifstream& istr, bool swap_bytes);

	static void skip(std::ifstream& istr, bool swap_bytes);
    
	//! Copy a bool integer header item and return a pointer to it.
	Hitem* copy() const {return new Hbool(value,comment);}

    private:
    
	bool value;
    
    };
  
    //! Storage of strings
  
    class Hstring : public Hitem {
    
    public:
    
	//! Default constructor
	Hstring() : Hitem(), value("") {} 

	//! Constructor from a value plus a comment
	Hstring(const std::string& sval, const std::string& com = "") : Hitem(com), value(sval) {}

	//! Constructor from file input
	Hstring(std::ifstream& istr, bool swap_bytes){read(istr, swap_bytes);}

	//! Destructor
	~Hstring(){}
    
	void  get_value(std::string &sval) const {sval = value;}
	std::string get_string() const {return value;}
    
	//! Returns name to use when referring to this type
	std::string type() const {return "string";}

	//! Type code for use in disk I/O
	HITEM_TYPE type_code() const {return HSTRING;}

	//! For tests of whether an Hitem is a directory or not
	bool is_a_dir() const  {return false;}
    
	void   set_value(const std::string &sval){value = sval;};
    
	void   print(std::ostream& s) const;

	void   write(std::ofstream& ostr) const;

	//! ASCII output
	void   write_ascii(std::ofstream& ostr) const;

	void   read(std::ifstream& istr, bool swap_bytes);

	static void skip(std::ifstream& istr, bool swap_bytes);
    
	//! Copy a string integer header item and return a pointer to it.
	Hitem* copy() const {return new Hstring(value,comment);}

    private:
    
	std::string value;
    
    };
  
    //! Storage of directories
  
    class Hdirectory : public Hitem {
    
    public:
    
	//! Default constructor
	Hdirectory() : Hitem() {}

	//! Constructor from a value plus a comment
	Hdirectory(const std::string& com) : Hitem(com) {}

	//! Constructor from file input
	Hdirectory(std::ifstream& istr, bool swap_bytes){read(istr, swap_bytes);}

	//! Destructor
	~Hdirectory(){}
    
	//! Returns name to use when referring to this type
	std::string type() const {return "directory";} 

	//! Type code for use in disk I/O
	HITEM_TYPE type_code() const {return HDIRECTORY;}

	//! For tests of whether an Hitem is a directory or not
	bool is_a_dir() const  {return true;}

	void  get_value(std::string &sval) const;

	std::string get_string() const;

	void   set_value(const std::string &sval);
    
	void   print(std::ostream& s) const;

	void   write(std::ofstream& ostr) const;

	//! ASCII output
	void   write_ascii(std::ofstream& ostr) const;

	void   read(std::ifstream& istr, bool swap_bytes);

	static void skip(std::ifstream& istr, bool swap_bytes);
    
	//! Copy a directory integer header item and return a pointer to it.
	Hitem* copy() const {return new Hdirectory(comment);}

    };
  
    //! Storage of Dates
  
    class Hdate : public Hitem {
    
    public:
    
	//! Default constructor
	Hdate() : Hitem(), value() {}

	//! Constructor from a value plus a comment
	Hdate(const Date& dval, const std::string& com = "") : Hitem(com), value(dval) {}

	//! Constructor from file input
	Hdate(std::ifstream& istr, bool swap_bytes){read(istr, swap_bytes);}

	//! Destructor
	~Hdate(){}
    
	void get_value(Date &dval) const {dval = value;}
	const Date& get_date() const {return value;}

	void get_value(std::string& sval) const;
	std::string get_string() const;
	void set_value(const std::string& sval);
    
	//! Returns name to use when referring to this type
	std::string type() const {return "date";}

	//! Type code for use in disk I/O
	HITEM_TYPE type_code() const {return HDATE;}

	//! For tests of whether an Hitem is a directory or not
	bool is_a_dir() const {return false;}
    
	void  set_value(const Date &dval){value = dval;}
	void  set_value(const Time &tval){value = Date(tval);}
    
	void  print(std::ostream& s) const;

	void  write(std::ofstream& ostr) const;

	//! ASCII output
	void   write_ascii(std::ofstream& ostr) const;

	void  read(std::ifstream& istr, bool swap_bytes);

	static void skip(std::ifstream& istr, bool swap_bytes);
    
	//! Copy a Date integer header item and return a pointer to it.
	Hitem* copy() const {return new Hdate(value,comment);}

    private:
    
	Date value;
    
    };
  
    //! Storage of Times
  
    class Htime : public Hitem {
    
    public:
    
	//! Default constructor
	Htime() : Hitem(), value() {} 

	//! Constructor from a value plus a comment
	Htime(const Time& tval, const std::string& com = "") : Hitem(com), value(tval) {}

	//! Constructor from file input
	Htime(std::ifstream& istr, bool swap_bytes){read(istr, swap_bytes);}

	//! Destructor
	~Htime(){}
    
	void get_value(Time &tval) const {tval = value;}

	const Time& get_time() const {return value;}

	double get_double() const {return value.mjd();}

	void get_value(std::string& sval) const;

	std::string get_string() const;

	void set_value(const std::string& sval);
    
	//! Returns name to use when referring to this type
	std::string type() const {return "time";}

	//! Type code for use in disk I/O
	HITEM_TYPE type_code() const {return HTIME;}

	//! For tests of whether an Hitem is a directory or not
	bool is_a_dir() const {return false;}
    
	void  set_value(const Time &tval){value = tval;}
    
	void  print(std::ostream& s) const;

	void  write(std::ofstream& ostr) const;

	//! ASCII output
	void   write_ascii(std::ofstream& ostr) const;

	void  read(std::ifstream& istr, bool swap_bytes);

	static void skip(std::ifstream& istr, bool swap_bytes);
    
	//! Copy a Time integer header item and return a pointer to it.
	Hitem* copy() const {return new Htime(value,comment);}

    private:
    
	Time value;
    
    };
  
    //! Storage of Positions
  
    class Hposition : public Hitem {
    
    public:
    
	//! Default constructor
	Hposition() : Hitem(), value() {} 

	//! Constructor from a value plus a comment
	Hposition(const Position& pval, const std::string& com = "") : Hitem(com), value(pval) {}

	//! Constructor from file input
	Hposition(std::ifstream& istr, bool swap_bytes){read(istr, swap_bytes);}

	//! Destructor
	~Hposition(){}
    
	void get_value(Position &pval) const {pval = value;}
	const Position& get_position() const {return value;}

	void get_value(std::string& sval) const;
	std::string get_string() const;
	void set_value(const std::string& sval);
    
	//! Returns name to use when referring to this type
	std::string type() const {return "position";}

	//! Type code for use in disk I/O
	HITEM_TYPE type_code() const {return HPOSITION;}

	//! For tests of whether an Hitem is a directory or not
	bool is_a_dir() const {return false;}
    
	void  set_value(const Position &pval){value = pval;}
    
	void  print(std::ostream& s) const;

	void  write(std::ofstream& ostr) const;

	//! ASCII output
	void   write_ascii(std::ofstream& ostr) const;

	void  read(std::ifstream& istr, bool swap_bytes);

	static void skip(std::ifstream& istr, bool swap_bytes);
    
	//! Copy a Position integer header item and return a pointer to it.
	Hitem* copy() const {return new Hposition(value,comment);}

    private:
    
	Position value;
    
    };

    //! Storage of Telescopes
  
    class Htelescope : public Hitem {
    
    public:
    
	//! Default constructor
	Htelescope() : Hitem(), value() {} 

	//! Constructor from a value plus a comment
	Htelescope(const Telescope& tval, const std::string& com = "") : Hitem(com), value(tval) {}

	//! Constructor from file input
	Htelescope(std::ifstream& istr, bool swap_bytes){read(istr, swap_bytes);}

	//! Destructor
	~Htelescope(){}
    
	void get_value(Telescope &tval) const {tval = value;}
	const Telescope& get_telescope() const {return value;}

	void get_value(std::string& sval) const;
	std::string get_string() const;
	void set_value(const std::string& sval);
    
	//! Returns name to use when referring to this type
	std::string type() const {return "telescope";}

	//! Type code for use in disk I/O
	HITEM_TYPE type_code() const {return HTELESCOPE;}

	//! For tests of whether an Hitem is a directory or not
	bool is_a_dir() const {return false;}
    
	void  set_value(const Telescope &tval){value = tval;}
    
	void  print(std::ostream& s) const;

	void  write(std::ofstream& ostr) const;

	//! ASCII output
	void   write_ascii(std::ofstream& ostr) const;

	void  read(std::ifstream& istr, bool swap_bytes);

	static void skip(std::ifstream& istr, bool swap_bytes);
    
	//! Copy a Telescope integer header item and return a pointer to it.
	Hitem* copy() const {return new Htelescope(value,comment);}

    private:
    
	Telescope value;
    
    };
  
    //! Storage of vectors of doubles
  
    class Hdvector : public Hitem {
    
    public:
    
	//! Default constructor
	Hdvector() : Hitem(), value() {} 

	//! Constructor from a value plus a comment
	Hdvector(const std::vector<double>& val, const std::string& com = "") : Hitem(com), value(val) {}

	//! Constructor from file input
	Hdvector(std::ifstream& istr, bool swap_bytes){read(istr, swap_bytes);}

	//! Destructor
	~Hdvector(){}
    
	void get_value(std::vector<double>& val) const {val = value;}

	const std::vector<double>& get_dvector() const {return value;}

	void  get_value(std::string &sval) const;

	std::string get_string() const;

	void  set_value(const std::string &sval);

	//! For tests of whether an Hitem is a directory or not
	bool is_a_dir() const {return false;}

	//! Type code for use in disk I/O
	HITEM_TYPE type_code() const {return HDVECTOR;}

	//! Returns name to use when referring to this type
	std::string type() const {return "dvector";}
    
	void  set_value(const std::vector<double>& val){value = val;}
    
	void  print(std::ostream& s) const;

	void  write(std::ofstream& ostr) const;

	//! ASCII output
	void   write_ascii(std::ofstream& ostr) const;

	void  read(std::ifstream& istr, bool swap_bytes);

	static void skip(std::ifstream& istr, bool swap_bytes);
    
	//! Copy a vector of double header item and return a pointer to it.
	Hitem* copy() const {return new Hdvector(value,comment);}

    private:
    
	std::vector<double> value;

    };

    std::ostream& operator<<(std::ostream& s, const Subs::Hitem* obj);

    //! Storage of unsigned chars

    class Huchar : public Hitem {
    
    public:

	//! Default constructor
	Huchar() : Hitem(), value(0) {} ;

	//! Constructor from a value plus a comment
	Huchar(unsigned char cval, const std::string& com="") : Hitem(com), value(cval) {};

	//! Constructor from file input
	Huchar(std::ifstream& istr, bool swap_bytes){read(istr, swap_bytes);}

	//! Destructor
	~Huchar(){};  
    
	//! Returns name to use when referring to this type
	std::string type() const {return "unsigned char";};

	//! Type code for use in disk I/O
	HITEM_TYPE type_code() const {return HUCHAR;}

	//! For tests of whether an Hitem is a directory or not
	bool is_a_dir() const {return false;};

	//! Gets value as an unsigned char
	void  get_value(unsigned char &cval) const {cval = value;};

	//! Gets value as an unsigned char
	unsigned char get_uchar() const {return value;};

	//! Sets value    
	void  set_value(const unsigned char &cval){value = cval;};

	void get_value(std::string &sval) const;

	std::string get_string() const;

	void set_value(const std::string& sval);

	//! Prints value
	void   print(std::ostream& s) const;

	//! Binary output
	void   write(std::ofstream& ostr) const;

	//! ASCII output
	void   write_ascii(std::ofstream& ostr) const;

	//! Binary input
	void   read(std::ifstream& istr, bool swap_bytes);

	//! Skip in a binary file
	static void skip(std::ifstream& istr, bool swap_bytes);

	//! Copy an character header item and return a pointer to it.
	Hitem* copy() const {return new Huchar(value,comment);}

    private:
    
	UCHAR value;
    
    };


    //! Storage of long ints
  
    class Hlint : public Hitem {
    
    public:
    
	//! Default constructor
	Hlint() : Hitem(), value(0) {};

	//! Constructor from a value plus a comment
	Hlint(int ival, const std::string& com="") : Hitem(com), value(ival) {};

	//! Constructor from file input
	Hlint(std::ifstream& istr, bool swap_bytes){read(istr, swap_bytes);}

	//! Destructor
	~Hlint(){};
    
	//! Gets value as a long int
	void   get_value(long int &ival) const  {ival = (long int)(value);};

	//! Gets value as an int
	void   get_value(int &ival) const  {ival = int(value);};
    
	//! Gets value as a double
	void   get_value(double &dval) const {dval = double(value);};
    
	//! Gets value as a float
	void   get_value(float &fval) const {fval = float(value);};
    
	//! Gets value as a long int
	long int get_lint() const {return (long int)(value);};
    
	//! Gets value as a double
	double get_double() const {return double(value);};
    
	//! Gets value as a float
	float  get_float() const {return float(value);};

	void get_value(std::string& sval) const;

	std::string get_string() const;

	void set_value(const std::string& sval);

	//! Returns name to use when referring to this type
	std::string type() const {return "long int";}

	//! Type code for use in disk I/O
	HITEM_TYPE type_code() const {return HLINT;}

	//! For tests of whether an Hitem is a directory or not
	bool is_a_dir() const  {return false;}
    
	void   set_value(const int &ival){value = ival;};
    
	void   print(std::ostream& s) const;

	void   write(std::ofstream& ostr) const;

	//! ASCII output
	void   write_ascii(std::ofstream& ostr) const;

	void   read(std::ifstream& istr, bool swap_bytes);

	static void skip(std::ifstream& istr, bool swap_bytes);

	//! Copy an integer header item and return a pointer to it.
	Hitem* copy() const {return new Hlint(value,comment);}
    
    private:
    
	// store as 32-bit
	INT4 value;
    
    };
  
    //! Storage of unsigned ints
  
    class Hulint : public Hitem {
    
    public:
    
	//! Default constructor
	Hulint() : Hitem(), value(0) {};

	//! Constructor from a value plus a comment
	Hulint(unsigned int ival, const std::string& com="") : Hitem(com), value(ival) {};

	//! Constructor from file input
	Hulint(std::ifstream& istr, bool swap_bytes){read(istr, swap_bytes);}

	//! Destructor
	~Hulint(){};
    
	void   get_value(unsigned long int &ival) const  {ival = (unsigned long int)(value);};

	void   get_value(double &dval) const {dval = double(value);};

	void   get_value(float &fval) const {fval = float(value);};

	unsigned long int get_ulint() const {return (unsigned long int)(value);};

	double get_double() const {return double(value);};

	float  get_float() const {return float(value);};

	void get_value(std::string& sval) const;

	std::string get_string() const;

	void set_value(const std::string& sval);
    
	//! Returns name to use when referring to this type
	std::string type() const {return "unsigned long int";}

	//! Type code for use in disk I/O
	HITEM_TYPE type_code() const {return HULINT;}

	//! For tests of whether an Hitem is a directory or not
	bool is_a_dir() const  {return false;}
    
	void   set_value(const unsigned long int &ival){value = ival;};
    
	void   print(std::ostream& s) const;

	void   write(std::ofstream& ostr) const;

	//! ASCII output
	void   write_ascii(std::ofstream& ostr) const;

	void   read(std::ifstream& istr, bool swap_bytes);

	static void skip(std::ifstream& istr, bool swap_bytes);

	//! Copy an unsigned integer header item and return a pointer to it.
	Hitem* copy() const {return new Hulint(value,comment);}

    
    private:
    
	UINT4 value;
    
    };

    //! Storage of unsigned short ints
  
    class Husint : public Hitem {
    
    public:
    
	//! Default constructor
	Husint() : Hitem(), value(0) {};

	//! Constructor from a value plus a comment
	Husint(unsigned short int ival, const std::string& com="") : Hitem(com), value(ival) {};

	//! Constructor from file input
	Husint(std::ifstream& istr, bool swap_bytes){read(istr, swap_bytes);}

	//! Destructor
	~Husint(){};
    
	void   get_value(unsigned int &ival) const  {ival = value;};
	void   get_value(double &dval) const {dval = double(value);};
	void   get_value(float &fval) const {fval = float(value);};
	unsigned int get_uint() const {return value;};
	double get_double() const {return double(value);};
	float  get_float() const {return float(value);};

	void get_value(std::string& sval) const;
	std::string get_string() const;
	void set_value(const std::string& sval);
    
	//! Returns name to use when referring to this type
	std::string type() const {return "unsigned short int";}

	//! Type code for use in disk I/O
	HITEM_TYPE type_code() const {return HUSINT;}

	//! For tests of whether an Hitem is a directory or not
	bool is_a_dir() const  {return false;}
    
	void   set_value(const unsigned int &ival){value = ival;};
    
	void   print(std::ostream& s) const;

	void   write(std::ofstream& ostr) const;

	//! ASCII output
	void   write_ascii(std::ofstream& ostr) const;

	void   read(std::ifstream& istr, bool swap_bytes);

	static void skip(std::ifstream& istr, bool swap_bytes);

	//! Copy an unsigned integer header item and return a pointer to it.
	Hitem* copy() const {return new Husint(value,comment);}

    
    private:
    
	UINT2 value;
    
    };


    //! Storage of vectors of integers
  
    class Hivector : public Hitem {
    
    public:
    
	//! Default constructor
	Hivector() : Hitem(), value() {} 

	//! Constructor from a value plus a comment
	Hivector(const std::vector<int>& val, const std::string& com = "") : Hitem(com), value(val) {}

	//! Constructor from file input
	Hivector(std::ifstream& istr, bool swap_bytes){read(istr, swap_bytes);}

	//! Destructor
	~Hivector(){}
    
	void get_value(std::vector<int>& val) const {val = value;}

	const std::vector<int>& get_ivector() const {return value;}

	void  get_value(std::string &sval) const;

	std::string get_string() const;

	void  set_value(const std::string &sval);

	//! For tests of whether an Hitem is a directory or not
	bool is_a_dir() const {return false;}

	//! Type code for use in disk I/O
	HITEM_TYPE type_code() const {return HIVECTOR;}

	//! Returns name to use when referring to this type
	std::string type() const {return "ivector";}
    
	void  set_value(const std::vector<int>& val){value = val;}
    
	void  print(std::ostream& s) const;

	void  write(std::ofstream& ostr) const;

	//! ASCII output
	void   write_ascii(std::ofstream& ostr) const;

	void  read(std::ifstream& istr, bool swap_bytes);

	static void skip(std::ifstream& istr, bool swap_bytes);
    
	//! Copy a vector of integers header item and return a pointer to it.
	Hitem* copy() const {return new Hivector(value,comment);}

    private:
    
	std::vector<int> value;

    };


    //! Storage of vectors of floats
  
    class Hfvector : public Hitem {
    
    public:
    
	//! Default constructor
	Hfvector() : Hitem(), value() {} 

	//! Constructor from a value plus a comment
	Hfvector(const std::vector<float>& val, const std::string& com = "") : Hitem(com), value(val) {}

	//! Constructor from file input
	Hfvector(std::ifstream& istr, bool swap_bytes){read(istr, swap_bytes);}

	//! Destructor
	~Hfvector(){}
    
	void get_value(std::vector<float>& val) const {val = value;}

	const std::vector<float>& get_fvector() const {return value;}

	void  get_value(std::string &sval) const;

	std::string get_string() const;

	void  set_value(const std::string &sval);

	//! For tests of whether an Hitem is a directory or not
	bool is_a_dir() const {return false;}

	//! Type code for use in disk I/O
	HITEM_TYPE type_code() const {return HFVECTOR;}

	//! Returns name to use when referring to this type
	std::string type() const {return "fvector";}
    
	void  set_value(const std::vector<float>& val){value = val;}
    
	void  print(std::ostream& s) const;

	void  write(std::ofstream& ostr) const;

	//! ASCII output
	void   write_ascii(std::ofstream& ostr) const;

	void  read(std::ifstream& istr, bool swap_bytes);

	static void skip(std::ifstream& istr, bool swap_bytes);
    
	//! Copy a vector of float header item and return a pointer to it.
	Hitem* copy() const {return new Hfvector(value,comment);}

    private:
    
	std::vector<float> value;

    };

};


  
#endif










