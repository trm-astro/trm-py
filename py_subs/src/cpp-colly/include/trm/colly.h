#ifndef TRM_MOLLY
#define TRM_MOLLY

//
// Header file for C++ molly routines
//

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include "trm/subs.h"
#include "trm/header.h"
#include "trm/position.h"

namespace Colly {

    //! Name of environment variable which can be set to specify location of default files
    const char COLLY_ENV[]         = "COLLY_ENV";
  
    //! Standard name of directory for default files if environment variable not set.
    const char COLLY_DIR[]         = ".colly";

    bool read_molly_head(std::ifstream& ist, Subs::Header& head, std::vector<std::string>& original,
			 bool warn=false);

    void write_molly_head(std::ofstream& ost, Subs::Header head, 
			  const std::vector<std::string>& original, bool old=true);

    void rejig_molly_head(Subs::Header& head, bool warn=false);
  
    void skip_molly_data(std::ifstream& ist, Subs::Header& head);
  
    bool skip_bytes(std::ifstream& ist, size_t nbytes);
  
    class Colly_Error : public std::string {
    public:

	//! Default constructor
	Colly_Error() :  std::string() {}

	//! Constructor from a string
	Colly_Error(const std::string& s) : std::string(s) {}

    };
  
    // Information kept on each spectrum including
    // its position, date, name and group (initialised to 0)

    struct Info {
	Info() : group(0) {}
	Subs::Position pos;
	Subs::Date date;
	std::string name;
	size_t group;
    };

    //! Recursive group finder
    void find_groups(std::vector<Info>& info, const Subs::Position& pos, 
		     const size_t& ngroup, const double& thresh);

    // minimal structure to store a molly spectrum.
    // Default constructor, destructor, copy constructor
    // and a read routine to set data. Stores header
    // and data in a charactr array.
  
    struct molly {

	molly() : head(), original(), nbytes(0), buff(0) {}
	molly(const molly& mol) : 
	    head(mol.head), original(mol.original), nbytes(mol.nbytes){
	    buff = new char [nbytes];
	    for(size_t i=0; i<nbytes; i++) buff[i] = mol.buff[i];
	}  
	~molly(){if(buff) delete[] buff;}
  
	void read(std::ifstream& ist, size_t ndata){    
	    if(buff && ndata != nbytes){
		delete[] buff;
		buff = 0;
	    }
	    if(!buff)
		buff = new char [ndata];
	    nbytes = ndata;
	    if(!(ist.read(buff,nbytes)))
		throw std::string("Failed to read bytes in molly::read");
	}
    
	Subs::Header head;
	std::vector<std::string> original;
	size_t nbytes;
	char* buff;
    };
}

#endif








