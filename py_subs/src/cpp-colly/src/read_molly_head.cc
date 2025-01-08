#include <fstream>
#include <vector>
#include <string>
#include "trm/subs.h"
#include "trm/header.h"
#include "trm/colly.h"

/**
 * read_molly_head reads in the headers and arc poly
 * (if set) of a molly spectrum. It assumes that the stream is
 * just positioned at the start of the header and leves it just
 * as the data record starts. This routine reads in the header items
 * in a virtually unaltered format (exception: any "Epoch"s are corrected
 * to "Equinox" and "DEC" to "Dec"). If you want a nicer format,
 * use rejig_molly_head, but note that write_molly_head expects
 * unvarnished headers. read_molly_head also returns a string array 
 * which shows the original storage order of the molly headers and can be used by
 * write_molly_head in order to write out an almost-untouched version of the header. 
 *
 * A directory called Xtra is created which stores a few extras, 
 * such as the number of pixels and the start and end wavelengths.
 *
 * \param ist the input file stream which must have been opened in binary mode
 * \param head the header to read into. It will be cleared at the start of the function
 * \param original vector of string names in the original
 * order in which the header items were loaded (since they are stored
 * alphabetically in 'head')
 * \param warn true to print warnings of header items not found
 * although expected; false to suppress these warnings
 * \return The routine comes back with false if the end of file is
 * encountered. It throws Colly_Error objects if serious
 * errors are encountered. These contain messages that can be
 * printed.
 */

namespace Colly {
  std::string convert_molly_to_header(std::string name);
}

bool Colly::read_molly_head(std::ifstream& ist, Subs::Header& head, 
			    std::vector<std::string>& original, bool warn){
  
    int  fcode, npix, narc, nchar, ndoub, nintr, nfloat;
    char units[17]; 
    int j, nstart;

    head.clear(); // clear the header
    original.clear(); // clear the std::string std::vector

    if(!(ist.read((char *)&nstart, sizeof(int)))) return false;
    if(nstart != 44)
	throw Colly_Error("Colly::read_molly_head: First 4 bytes != 44");
    if(!(ist.read((char *)&fcode,sizeof(int)))) 
	throw Colly_Error("Colly::read_molly_head: Error reading fcode");
    if(!(ist.read(units,16)))
	throw Colly_Error("Colly::read_molly_head: Error reading units");
    units[16] = 0;
    if(!(ist.read((char *)&npix,sizeof(int))))
	throw Colly_Error("Colly::read_molly_head: Error reading npix");
    if(!(ist.read((char *)&narc,sizeof(int))))
	throw Colly_Error("Colly::read_molly_head: Error reading narc");

    try{
	head.set("Xtra",       new Subs::Hdirectory("Extra information on spectrum"));
	head.set("Xtra.FCODE", new Subs::Hint(fcode,"molly format code"));
	head.set("Xtra.UNITS", new Subs::Hstring(units, "Flux units"));
	head.set("Xtra.NPIX",  new Subs::Hint(npix,"Number of pixels"));
	head.set("Xtra.NARC",  new Subs::Hint(narc, "Number of arc coefficients"));
    }
    catch(std::string s){
	if(warn) std::cerr << s << std::endl;
    } 

#ifdef DEBUG
    std::cerr << "FCODE = " << fcode << ", UNITS = " << units
	      << ", NPIX = " << npix << std::endl;
#endif

    if(!(ist.read((char *)&nchar,sizeof(int))))
	throw Colly::Colly_Error("Colly::read_molly_head: Error reading nchar");
    if(!(ist.read((char *)&ndoub,sizeof(int))))
	throw Colly::Colly_Error("Colly::read_molly_head: Error reading ndoub");
    if(!(ist.read((char *)&nintr,sizeof(int))))
	throw Colly::Colly_Error("Colly::read_molly_head: Error reading nintr");
    if(!(ist.read((char *)&nfloat,sizeof(int))))
	throw Colly::Colly_Error("Colly::read_molly_head: Error reading nfloat");
    if(!(skip_bytes(ist,2*sizeof(int))))
	throw Colly::Colly_Error("Colly::read_molly_head: Error skipping bytes after first record");

#ifdef DEBUG
    std::cerr << "NCHAR = " << nchar << ", NDOUB = " << ndoub
	      << ", NINTR = " << nintr << ", NFLOAT = " << nfloat << std::endl;
#endif

    // Read in header item names, convert them to the right
    // format (no whitespace)

    char hname[17];
    hname[16] = 0;

    std::vector<std::string> cname(nchar);
    for(j=0; j<nchar; j++){
	if(!(ist.read(hname,16)))
	    throw Colly::Colly_Error("Colly::read_molly_head: Error reading char header name");
	cname[j] = convert_molly_to_header(hname);
	original.push_back(cname[j]);
    }

    // Read double with a few corrections

    std::vector<std::string> dname(ndoub);
    for(j=0; j<ndoub; j++){
	if(!(ist.read(hname,16)))
	    throw Colly::Colly_Error("Colly::read_molly_head: Error reading double header name");
	dname[j] = convert_molly_to_header(hname);

	if(dname[j] == "DEC") dname[j] = "Dec";
	if(dname[j] == "EQUINOX" || dname[j] == "Epoch") dname[j] = "Equinox";
	original.push_back(dname[j]);
    }

    std::vector<std::string> iname(nintr);
    for(j=0; j<nintr; j++){
	if(!(ist.read(hname,16)))
	    throw Colly::Colly_Error("Colly::read_molly_head: Error reading int header name");
	iname[j] = convert_molly_to_header(hname);
	original.push_back(iname[j]);
    }

    std::vector<std::string> fname(nfloat);
    for(j=0; j<nfloat; j++){
	if(!(ist.read(hname,16)))
	    throw Colly::Colly_Error("Colly::read_molly_head: Error reading float header name");
	fname[j] = convert_molly_to_header(hname);
	original.push_back(fname[j]);
    }

    if(!(skip_bytes(ist,2*sizeof(int))))
	throw Colly::Colly_Error("Colly::read_molly_head: Error skipping bytes after header names");

#ifdef DEBUG
    std::cerr << "Read header item names" << std::endl;
#endif

    // Now read and set the values, chopping
    // trailing blanks from strings.

    char citem[33];
    citem[32] = 0;
    for(j=0; j<nchar; j++){
	if(!(ist.read(citem,32)))
	    throw Colly::Colly_Error("Colly::read_molly_head: Error reading char item values");
	for(int k=31; k>=0; k--)
	    if(citem[k] == ' ')
		citem[k] = 0;
	    else
		break;
	head.set(cname[j], new Subs::Hstring(citem));
    }

    double d;
    for(j=0; j<ndoub; j++){
	if(!(ist.read((char *)&d,sizeof(double))))
	    throw Colly::Colly_Error("Colly::read_molly_head: Error reading double item values");
	head.set(dname[j], new Subs::Hdouble(d));
    }

    int i;
    for(j=0; j<nintr; j++){
	if(!(ist.read((char *)&i,sizeof(int))))
	    throw Colly::Colly_Error("Colly::read_molly_head: Error reading int item values");
	head.set(iname[j], new Subs::Hint(i));
    }

    float f;
    for(j=0; j<nfloat; j++){
	if(!(ist.read((char *)&f,sizeof(float))))
	    throw Colly::Colly_Error("Colly::read_molly_head: Error reading float item values");
	head.set(fname[j], new Subs::Hfloat(f));
    } 

    if(!(skip_bytes(ist,2*sizeof(int))))
	throw Colly::Colly_Error("Colly::read_molly_head: Error skipping bytes after reading header item values");

#ifdef DEBUG
    std::cerr << "Read header item values" << std::endl;
#endif

    if(narc != 0){
	double arc[abs(narc)];
	if(!(ist.read((char *)&arc,abs(narc)*sizeof(double))))
	    throw Colly::Colly_Error("Colly::read_molly_head: Error reading arc coefficients");
	double w1 = 0., w2 = 0.;
	for(i=0; i<abs(narc); i++){
	    w1 += arc[i]*pow(0.5/double(npix),i);
	    w2 += arc[i]*pow((npix+0.5)/double(npix),i);;
	}
	if(w1 > w2) std::swap(w1,w2);
	if(narc < 0){
	    w1 = exp(w1);
	    w2 = exp(w2);
	}

	double disp = (w2-w1)/npix;
	head.set("Xtra.ARC", 
		 new Subs::Hdvector(std::vector<double>(arc,arc+abs(narc)), 
				    "Arc coefficients"));
	head.set("Xtra.WLO", 
		 new Subs::Hdouble(w1, "Lower wavelength (A, at pixel 0.5)"));
	head.set("Xtra.WHI", 
		 new Subs::Hdouble(w2,"Upper wavelength (A, at pixel npix+0.5)"));
	head.set("Xtra.DISP",
		 new Subs::Hdouble(disp,"Dispersion (A/pix)"));

#ifdef DEBUG
	std::cerr << "Read arc coefficients" << std::endl;
#endif

    }

    // skip bytes at end of record

    if(!(skip_bytes(ist,sizeof(int))))
	throw Colly::Colly_Error("Colly::read_molly_head: Error skipping bytes after reading arc coefficients");

    // A few corrections ...

    // Correct for UTC >= 24.

    try{
	double utc = head["UTC"]->get_double();
	if(utc >= 24.){
	    Subs::Header::Hnode* dhit = head.find("Day");
	    Subs::Header::Hnode* mhit = head.find("Month");
	    Subs::Header::Hnode* yhit = head.find("Year");
	    int day    = dhit->value->get_int();
	    int month  = mhit->value->get_int();
	    int year   = yhit->value->get_int();
	    utc -= 24.;
	    Subs::Time time(day,Subs::Date::Month(month),year,utc);
	    time.add_day(1);
	    dhit->value->set_value(time.day());
	    mhit->value->set_value(int(time.month()));
	    yhit->value->set_value(time.year());
	    head["UTC"]->set_value(utc);
	}
    }
    catch(const std::string& s){
	if(warn) std::cerr << s << std::endl;
    }

    return true;
}

namespace Colly {
  std::string convert_molly_to_header(std::string name){
    std::string::size_type n1 = name.find_first_not_of(" \t");
    if(n1 == std::string::npos)
      throw Colly::Colly_Error("Blank name in std::string convert_molly_to_header(std::string)");
    std::string::size_type n2 = name.find_last_not_of(" \t");
    name = name.substr(n1,n2-n1+1);
    std::string::size_type p;
    while((p = name.find(" ")) != std::string::npos)
      name.replace(p,1,"_");
    while((p = name.find("\t")) != std::string::npos)
      name.replace(p,1,"_");
    return name;
  }
}

  






