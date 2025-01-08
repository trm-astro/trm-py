/*

!!begin
!!author  T.R. Marsh
!!created 20 August 2001
!!revised 07 Jan 2008
!!descr   corrects positions of molly spectra.
!!root    correct
!!index   correct
!!class   Programs
!!css     style.css
!!head1   correct - corrects positions of molly spectra

!!emph{correct} reads in a series of molly spectra which it then
groups in the same way as !!ref{group.html}{group}. It also reads in a
file of positions which it then uses to change the positions of groups
which match. This provides a way to correct the often-erroneous
positions found in archives of spectra and thus more secure selection
later on (as well as better heliocentric correction). Every file read
in is duplicated with the addition of .zzz at the end of the file
name. It is then up to the user to check these are OK and then move
them over the top of the old ones.

!!head2 Invocation

correct corfile delta warn list

!!head2 Arguments

!!table

!!arg{corfile}{File of positions to act as group centres each of which
comes along with a new position. See below for the format.}

!!arg{delta}{Angle to define groups in arcminutes.}

!!arg{warn}{"true" or "false" to switch on warnings when reading the
headers or not. Although these can be useful, they can also be irritating.
You should probably run it once with warn=true to see that things
are working approximately correctly.}

!!arg{list}{List of molly file names. e.g. use "find . -name "*.mol" > list" to
create a list of all molly files from all sub-driectories
of the present working directory.}

!!table

!!head2 Format of the correction file

The correction file should have lines such as!!break
# My favourite target!!break
10 02 32.3 -00 45 53 2000 10 02 34.4 -00 44 53 2000!!break

i.e. comment lines start with #, while a changing line consists of
the old RA, Dec, Equinox, followed by the new. Do not use colons between
the numbers.

!!end

*/

#include <cstdlib>
#include <iomanip>
#include <fstream>
#include <vector>
#include "trm/subs.h"
#include "trm/input.h"
#include "trm/header.h"
#include "trm/position.h"
#include "trm/date.h"
#include "trm/colly.h"

int main(int argc, char* argv[]){

    const double CEQUINOX = 2000.;
    Subs::Date::print_method = 2;

    std::cout << "THIS IS CURRENTLY BROKEN; NEEDS PRECESSION ADDED" << std::endl;
    try{
    
	// Construct Input object
	Subs::Input input(argc, argv, Colly::COLLY_ENV, Colly::COLLY_DIR);

	// Define inputs
	input.sign_in("corfile",  Subs::Input::LOCAL,  Subs::Input::PROMPT);
	input.sign_in("delta",    Subs::Input::GLOBAL, Subs::Input::PROMPT);
	input.sign_in("warn",     Subs::Input::LOCAL,  Subs::Input::PROMPT);
	input.sign_in("flist",    Subs::Input::GLOBAL, Subs::Input::PROMPT);

	std::string scorfile;
	input.get_value("corfile",  scorfile, "correct", "file of position corrections");
	float delta;
	input.get_value("delta", delta,  5.f, 0.f, 10000.f, "angle in arcminutes to define groups");
	bool warn;
	input.get_value("warn", warn,  true, "issue warnings when reading headers or not");
	std::string flist;
	input.get_value("flist",  flist, "list", "list of molly files");

	// Read in positions, skipping lines starting with #
	// switch to common equinox

	std::ifstream corfile(scorfile.c_str());
	if(!corfile)
	    throw std::string("Failed to open ") + scorfile;
	std::vector<Subs::Position> opos, npos;
	Subs::Position op, np;
	char c;
	while(corfile){
	    c = corfile.peek();
	    if(corfile){
		if(c != '#'){
		    if(corfile >> op >> np){
			std::cerr << op << " ---> " << np << std::endl;
//			op.precess(CEQUINOX);
			opos.push_back(op);
			npos.push_back(np);
		    }
		}
		while(corfile.get(c) && c != '\n');
	    }
	}
	corfile.close();
	if(opos.size() == 0){
	    throw std::string("No positions read.");
	}else{
	    std::cout << opos.size() << " positions read.\n";
	}

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
	Subs::Header::Hnode *dayp, *monthp = NULL;
	Subs::Header::Hnode *yearp = NULL, *rap = NULL, *decp = NULL, *eqp = NULL;
	std::vector<Colly::Info> info;
	std::vector<std::string> original;
	Colly::Info tinf;
	size_t count, nspec = 0;

	// Read in all molly files

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
	    

			info.push_back(tinf);

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

	// Now precess all coordinates to a common equinox
    
//	for(size_t i=0; i<info.size(); i++)
//	    info[i].pos.precess(CEQUINOX);

//	std::cout << "Precessed to a common equinox." << std::endl;

	// Allocate to groups

	size_t ngroup = 0;
	double thresh = cos(Constants::TWOPI*delta/60./360.);
    
	std::cout << "Now allocating groups." << std::endl;

	// Recursive group allocator

	for(size_t i=0; i<info.size(); i++)
	    if(!info[i].group)
		find_groups(info,info[i].pos,++ngroup,thresh);

	std::cout << ngroup << " target groups were found." << std::endl; 
    
	// Now work out which groups correspond to which 
	// positions.

	if(ngroup){
	    std::vector<int> gindex(ngroup,-1);
      
	    // Identify positions for each group if possible

	    int mpos = 0;
	    for(size_t j=0; j<ngroup; j++){
		for(size_t i=0; i<info.size(); i++){
		    if(info[i].group-1 == j){
			for(size_t k=0; k<opos.size(); k++){
			    if(dot(info[i].pos,opos[k]) > thresh){
				if(gindex[j] == -1){
				    gindex[j] = k;
				    mpos++;
				}else if(gindex[j] != int(k)){
				    std::string message = "Multiple positions found for group ";
				    message += j+1;
				    throw message;
				}
			    }
			}
		    }
		}
	    }

	    std::cout << mpos << " groups have new positions.\n";

	    // Now read in every file all over again, but this time change
	    // the position where necessary and dump out to a new file.

	    int nok = 0, changed = 0;
	    std::string newfile;
	    Colly::molly spectrum;
      
	    for(size_t j=0; j<file.size(); j++){
	
		newfile = file[j] + ".zzz";  
		if(warn){
		    std::cerr << "Now reading file " << j+1 << " = " << file[j] 
			      << " and writing file = " << newfile << std::endl;
		}
	  
	  
		std::ifstream fin(file[j].c_str(), std::ios::in | std::ios::binary);
		std::ifstream ftest(newfile.c_str(), std::ios::in | std::ios::binary);
		if(ftest)
		    throw std::string("File " + newfile + " already exists.");

		std::ofstream fout(newfile.c_str(), std::ios::out | std::ios::binary);
		if(!fout)
		    throw std::string("Failed to open " + newfile + " for output.");

		count = 0;
		try{
		    while(Colly::read_molly_head(fin,spectrum.head,spectrum.original,warn)){
			count++;
			if(warn)
			    std::cerr << "Read spectrum " << count << std::endl;
	    
#ifdef DEBUG
			std::cerr << spectrum.head << std::endl;
#endif
	    
			size_t nbytes;
			int fcode = spectrum.head["Xtra.FCODE"]->get_int();
			int npix  = spectrum.head["Xtra.NPIX"]->get_int();
			switch(fcode){
			    case 1: case 4:
				nbytes = npix*sizeof(float)+8;
				break;
			    case 2: case 5:
				nbytes = 2*npix*sizeof(float)+8;
				break;
			    case 3:
				nbytes = 3*npix*sizeof(float)+8;
				break;
			    default:
				std::string message = "Unrecognised molly format code = " + fcode;
				throw message;
			}
			spectrum.read(fin,nbytes);
	    
			if(spectrum.head.find("Day")->has_data()   &&
			   spectrum.head.find("Month")->has_data() &&
			   spectrum.head.find("Year")->has_data()  &&
			   spectrum.head.find("RA")->has_data()    &&
			   spectrum.head.find("Dec")->has_data()   &&
			   spectrum.head.find("Equinox")->has_data()){
	      
			    if(info[nok].group && gindex[info[nok].group-1] != -1){
				changed++;
				spectrum.head["RA"]->set_value(npos[gindex[info[nok].group-1]].ra());
				spectrum.head["Dec"]->set_value(npos[gindex[info[nok].group-1]].dec());
				spectrum.head["Equinox"]->set_value(npos[gindex[info[nok].group-1]].epoch());
			    }
	      
			    nok++;
			}
	    
			// write spectrum out
	    
			Colly::write_molly_head(fout,spectrum.head,spectrum.original);
			if(!(fout.write(spectrum.buff,spectrum.nbytes)))
			    throw std::string("Failed during writing of " + newfile);
	    
		    }
		}
		catch(const Colly::Colly_Error& error){
		    std::cerr << error << std::endl;
		    std::cerr << spectrum.head << std::endl;
		}
		catch(const std::string& message){
		    std::cerr << message << std::endl;
		}
		catch(...){
		    throw;
		}
		fin.close();
		fout.close();
		std::cout << "Written " << newfile << std::endl;
	    }
	    std::cout << "Changed a total of " << changed << " spectra." << std::endl;
	}
    }
    catch(const std::string& error){
	std::cerr << error << std::endl;
    }
}








