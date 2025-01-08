/*

!!begin
!!author  T.R. Marsh
!!created 25 June 2001
!!revised 07 Jan 2008
!!descr   groups molly spectra by position
!!root    group
!!index   group
!!class   Programs
!!css     style.css
!!head1   group - groups molly spectra by position

!!emph{group} reads in molly files and then groups spectra
by their position. It then reports the number and position of each
group along with its most common name and the first and last dates
on which the spectra were taken. It can also select spectra
covering specified wavelengths only.

It is not uncommon to encounter problems with group so it is  designed
to give fairly full error information when encountering problems. A
common one is a message about the first 4 bytes != 44. This is an
indication that the file being read is not a molly file. A common way
for this to happen is for it in fact to be a molly file from a machine
of different architecture. You may then need to convert the file.

!!head2 Invocation

group delta warn [-w wave] file1 file2 file3 ...

!!head2 Arguments

!!table

!!arg{delta}{Angle to define groups in arcminutes.}

!!arg{warn}{"true" or "false" to switch on warnings when reading the
headers or not. Although these can be useful, they can also be irritating.
You should probably run it once with warn=true to see that things
are working approximately correctly.}

!!arg{wave}{A wavelength that must be covered by a spectrum for it
to be considered. <= 0 to ignore}

!!arg{flist}{List of molly files.}

!!table

!!end

*/

#include <cstdlib>
#include <iomanip>
#include <fstream>
#include <vector>
#include <map>

#include "trm/subs.h"
#include "trm/constants.h"
#include "trm/input.h"
#include "trm/header.h"
#include "trm/position.h"
#include "trm/date.h"
#include "trm/colly.h"

// Information kept on each group 
//
// names (with numbers of each), a reference position
// and a date range,

struct Ginfo{
    std::map<std::string,size_t> names;
    Subs::Position repos;
    Subs::Date dmin, dmax;
};

int main(int argc, char* argv[]){

    const double CEPOCH = 2000.;
    Subs::Date::print_method = 2;

    std::cout << "PRECESSION BROKEN" << std::endl;

    try{
    
	// Construct Input object
	Subs::Input input(argc, argv, Colly::COLLY_ENV, Colly::COLLY_DIR);

	// Define inputs
	input.sign_in("delta",    Subs::Input::GLOBAL,  Subs::Input::PROMPT);
	input.sign_in("warn",     Subs::Input::LOCAL,  Subs::Input::PROMPT);
	input.sign_in("flist",    Subs::Input::GLOBAL, Subs::Input::PROMPT);

	float delta;
	input.get_value("delta", delta,  5.f, 0.f, 10000.f, "angle in arcminutes to define groups");
	bool warn;
	input.get_value("warn", warn,  true, "issue warnings when reading headers or not");
	double wave;
	input.get_value("wave", wave,  0., -10000., 1.e10, "wavelength that must be covered (<= 0 to ignore)");
	std::string flist;
	input.get_value("flist",  flist, "list", "list of molly files");


	// Read file names
	std::ifstream liststr(flist.c_str());
	if(!liststr)
	    throw std::string("Failed to open ") + flist;    
	std::vector<std::string> file;
	std::string fname;
	while(liststr >> fname)
	    file.push_back(fname);
	if(file.size() == 0)
	    throw std::string("No file names specified!");

	Subs::Header head; 
	Subs::Header::Hnode *namep, *w1p, *w2p = NULL, *dayp, *monthp = NULL;
	Subs::Header::Hnode *yearp = NULL, *rap = NULL, *decp = NULL, *eqp = NULL;
	std::vector<Colly::Info> info;
	std::vector<std::string> original;
	Colly::Info tinf;
	size_t count, nspec = 0;

	for(size_t j=0; j<file.size(); j++){
	    if(warn)
		std::cerr << "Now reading file " << j+1 << " = " << file[j] << std::endl;

	    std::ifstream f(file[j].c_str(), std::ios::in | std::ios::binary);
	    count = 0;
	    try{
		while(Colly::read_molly_head(f,head,original,warn)){
		    count++;
		    if(warn)
			std::cerr << "Read spectrum " << count << std::endl;
	  
#ifdef DEBUG
		    std::cerr << head << std::endl;
#endif
	  
		    nspec++;
		    if(nspec % 1000 == 0)
			std::cerr << "Read " << nspec << " spectra." << std::endl;
	  
		    if(
			(dayp   = head.find("Day"))->has_data()   &&
			(monthp = head.find("Month"))->has_data() &&
			(yearp  = head.find("Year"))->has_data()  &&
			(rap    = head.find("RA"))->has_data()    &&
			(decp   = head.find("Dec"))->has_data()   &&
			(eqp    = head.find("Equinox"))->has_data()){

			tinf.date = Subs::Date(dayp->value->get_int(),
					       Subs::Date::Month(monthp->value->get_int()),
					       yearp->value->get_int());

			tinf.pos  = Subs::Position(rap->value->get_double(),
						   decp->value->get_double(),
						   eqp->value->get_double());
	    
			if((namep = head.find("Object"))->has_data() && 
			   namep->value->type() == "string")
			    tinf.name = namep->value->get_string();
			else
			    tinf.name = "NULL";
	      
			if(wave > 0.){
			    if((w1p = head.find("Xtra.WLO"))->has_data() &&
			       (w2p = head.find("Xtra.WHI"))->has_data()){
		      
				if(w1p->value->get_double() < wave && w2p->value->get_double() < wave)
				    info.push_back(tinf);
		      
			    }else if(warn){
				std::cerr << "Spectrum " << count << " of " << file[j] 
					  << " has no wavelength limits." << std::endl;
			    }
		  
			}else{ 
			    info.push_back(tinf);
			}
		    }else{
			std::cerr << "Spectrum " << count << " of " << file[j] 
				  << " has no positional and/or timing information." << std::endl;
			if(warn) std::cerr << head << std::endl;
		    }
		    Colly::skip_molly_data(f,head);
		}
	    }
	    catch(const Colly::Colly_Error& error){
		std::cerr << error << std::endl;
		std::cerr << head << std::endl;
	    }
	    catch(const std::string& message){
		std::cerr << message << std::endl;
	    }
	    catch(...){
		throw;
	    }
	    f.close();
	}
	if(info.size()){
	    std::cout << "Finished file input. " << nspec
		      << " spectra read of which " << info.size() 
		      << " were valid." << std::endl;
	}else{
	    std::cout << "Finished file input. No valid spectra found." << std::endl;
	    exit(EXIT_FAILURE);
	}

	// Now precess all coordinates to a common epoch
    
//	for(size_t i=0; i<info.size(); i++)
//	    info[i].pos.precess(CEPOCH);

//	std::cout << "Precessed to a common epoch." << std::endl;

	// Allocate to groups
	size_t ngroup = 0;
	double thresh = cos(Constants::TWOPI*delta/60./360.);
    
	std::cout << "Now allocating groups." << std::endl;

	// Recursive group allocator

	for(size_t i=0; i<info.size(); i++)
	    if(!info[i].group)
		find_groups(info,info[i].pos,++ngroup,thresh);

	std::cout << ngroup << " target groups were found." << std::endl; 
    
	// Now compute bits of info on groups, store in 
	// vector of Ginfo objects. Also create a map based
	// on keys made from RA and Dec so that an ordered
	// output is given at the end.

	if(ngroup){
	    std::vector<Ginfo> ginfo(ngroup);
	    std::map<std::string,size_t> order;

	    typedef std::map<std::string,size_t>::const_iterator CI;
      
	    // Count names

	    for(size_t i=0; i<info.size(); i++)
		ginfo[info[i].group-1].names[info[i].name]++;

	    // Find reference position and initialise date range

	    for(size_t i=0; i<ngroup; i++){
		for(size_t j=0; j<info.size(); j++){
		    if(info[j].group-1 == i){
			ginfo[i].repos = info[j].pos;
			ginfo[i].dmin  = ginfo[i].dmax = info[j].date;
			order[info[j].pos.ra_dec()] = i;
			break;
		    }
		}

		// Finish computation of date range

		for(size_t j=0; j<info.size(); j++){
		    if(info[j].group-1 == i){
			if(ginfo[i].dmin > info[j].date){
			    ginfo[i].dmin = info[j].date;
			}else if(ginfo[i].dmax < info[j].date){
			    ginfo[i].dmax = info[j].date;
			}
		    }
		}
	    }

	    // Finally output results
      
	    size_t ng = 1, ind;
	    std::string name, name1, name2, name3;
	    size_t n, n1, n2, n3, ntot;
	    for(CI finp=order.begin(); finp!=order.end(); finp++){

		ind = finp->second;

		// Find 3 most common names

		n1 = n2 = n3 = ntot = 0;
		for(CI p=ginfo[ind].names.begin(); p!=ginfo[ind].names.end(); p++){
		    ntot += p->second;
		    if(p->second > n1){
			n3    = n2;
			name3 = name2;
			n2    = n1;
			name2 = name1;
			n1    = p->second;
			name1 = p->first;
		    }else if(p->second > n2){
			n3    = n2;
			name3 = name2;
			n2    = p->second;
			name2 = p->first;
		    }else if(p->second > n3){
			n3    = p->second;
			name3 = p->first;
		    }
		}

		// select one avoiding certain uninformative ones

		if((name1 == "ARC"  || name1 == "THAR" || name1 == "CUAR" ||
		    name1 == "CUNE" || name1 == "CuNe" || name1 == "CuAr") &&
		   n2){
		    if((name2 == "ARC"  || name2 == "THAR" || name2 == "CUAR" ||
			name2 == "CUNE" || name2 == "CuNe" || name2 == "CuAr") &&
		       n3){
			name = name3;
			n    = n3;
		    }else{
			name = name2;
			n    = n2;
		    }
		}else{
		    name = name1;
		    n    = n1;
		} 

		std::cout.setf(std::ios::right);
		std::cout << std::setfill(' ') << std::setw(5) << ng++ << ") " 
			  << name; 
		int l = 24-name.length();
		if(l > 0){
		    std::cout.width(l);
		    std::cout << " ";
		}
		std::cout << " [" << n << "/" << ntot << "]";
		int ndig = int(floor(log10(double(n)))+floor(log10(double(ntot))));
		l = 6-ndig;
		if(l > 0){
		    std::cout.width(l);
		    std::cout << " ";
		}
		std::cout << " " << ginfo[ind].repos << ", Dates: " 
			  << ginfo[ind].dmin << " to " << ginfo[ind].dmax << std::endl;
	    }
	}
    }
    catch(std::string error){
	std::cerr << error << std::endl;
    }
}







