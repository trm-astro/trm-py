/*

!!begin
!!author  T.R. Marsh
!!created 20 August 2001
!!revised 07 Jan 2008
!!descr   grabs molly spectra.
!!root    grab
!!index   grab
!!class   Programs
!!css     style.css
!!head1   grab - grabs molly spectra

!!emph{grab} grabs molly spectra from a series of files according to
their position and dumps them to a separate file.

!!head2 Invocation

grab newfile position delta warn wave flist

!!head2 Arguments

!!table
!!arg{newfile}{Name of new file to contain the grabbed spectra. It must
not already exist.}

!!arg{position}{RA, Dec, Equinox. This should be specified as a string as in
"10 12 32.4 -52 34 45.0 2000". This needs both the quotes and the sign, even
when positive.}

!!arg{delta}{Tolerance around position for spectra to be included (arcminutes).}

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
#include <string>
#include <vector>
#include "trm/subs.h"
#include "trm/input.h"
#include "trm/constants.h"
#include "trm/header.h"
#include "trm/position.h"
#include "trm/date.h"
#include "trm/colly.h"

int main(int argc, char* argv[]){

    const double CEQUINOX = 2000.;
    Subs::Date::print_method = 2;

    std::cout << "CURRENTLY BROKEN PRECESSION" << std::endl;

    try{

	// Construct Input object
	Subs::Input input(argc, argv, Colly::COLLY_ENV, Colly::COLLY_DIR);

	// Define inputs
	input.sign_in("newfile",  Subs::Input::LOCAL,  Subs::Input::PROMPT);
	input.sign_in("position", Subs::Input::GLOBAL, Subs::Input::PROMPT);
	input.sign_in("delta",    Subs::Input::GLOBAL, Subs::Input::PROMPT);
	input.sign_in("warn",     Subs::Input::LOCAL,  Subs::Input::PROMPT);
	input.sign_in("wave",     Subs::Input::LOCAL,  Subs::Input::PROMPT);
	input.sign_in("flist",    Subs::Input::GLOBAL, Subs::Input::PROMPT);

	std::string snewfile;
	input.get_value("newfile",  snewfile, "newfile", "name of file for grabbed spectra");
	std::string sposition;
	input.get_value("position",  sposition, "10 12 00.23 -01 12 20.0 2000", "position to locate spectra");
	Subs::Position pos(sposition);
	float delta;
	input.get_value("delta", delta,  5.f, 0.f, 10000.f, "angle in arcminutes to define groups");
	bool warn;
	input.get_value("warn", warn,  true, "issue warnings when reading headers or not");
	double wave;
	input.get_value("wave", wave,  0., -10000., 1.e10, "wavelength that must be covered (<= 0 to ignore)");
	std::string flist;
	input.get_value("flist",  flist, "list", "list of molly files");

	std::cerr << "Searching around position: " << pos << std::endl;

//	pos.precess(CEQUINOX);
	const double thresh = cos(Constants::TWOPI*delta/60./360.);

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

	// Open output file

	std::ifstream ftest(snewfile.c_str(), std::ios::in | std::ios::binary);
	if(ftest)
	    throw std::string("File " + snewfile + " already exists.");

	std::ofstream fout(snewfile.c_str(), std::ios::out | std::ios::binary);
	if(!fout)
	    throw std::string("Failed to open ") + snewfile + std::string(". Does it exist already?");

	Subs::Header::Hnode *w1p, *w2p = NULL, *rap, *decp = NULL, *eqp = NULL;
	std::vector<std::string> original;
	Subs::Position tpos;
	int count, nspec = 0;
	Colly::molly spectrum;

	// Read in all molly files

	for(size_t j=0; j<file.size(); j++){
	    if(warn)
		std::cerr << "Now reading file " << j+1 << " = " << file[j] << std::endl;

	    std::ifstream fin(file[j].c_str(), std::ios::in | std::ios::binary);
	    count = 0;
	    try{
		while(Colly::read_molly_head(fin,spectrum.head,spectrum.original,warn)){
		    count++;
		    if(warn)
			std::cerr << "Read spectrum " << count << std::endl;
	  
#ifdef DEBUG
		    std::cerr << spectrum.head << std::endl;
#endif
	  
		    nspec++;
		    if(nspec % 1000 == 0)
			std::cerr << "Read " << nspec << " spectra." << std::endl;
	  
		    if(
			(rap    = spectrum.head.find("RA"))->has_data()  &&
			(decp   = spectrum.head.find("Dec"))->has_data() &&
			(eqp    = spectrum.head.find("Equinox"))->has_data()){

			tpos  = Subs::Position(rap->value->get_double(),
					       decp->value->get_double(),
					       eqp->value->get_double());
//			tpos.precess(CEQUINOX);

			if(dot(tpos,pos) > thresh){
			    bool ok_to_read = true;
			    if(wave > 0.){
				if((w1p = spectrum.head.find("Xtra.WLO"))->has_data() &&
				   (w2p = spectrum.head.find("Xtra.WHI"))->has_data()){
				    if(w1p->value->get_double() > wave || w2p->value->get_double() < wave)
					ok_to_read = false;

				}else if(warn){
				    std::cerr << "Spectrum " << count << " of " << file[j] 
					      << " has no wavelength limits." << std::endl;
				    ok_to_read = false;
				}
			    }
			    if(ok_to_read){
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

				// write it out
		
				Colly::write_molly_head(fout,spectrum.head,spectrum.original);
				if(!(fout.write(spectrum.buff,spectrum.nbytes)))
				    throw std::string("Failed during writing of ") + std::string(argv[1]);
				if(warn)
				    std::cerr << "Dumped spectrum " << count << " of file = " 
					      << file[j] << std::endl;
			    }else{
				Colly::skip_molly_data(fin,spectrum.head);
			    }
			}else{
			    Colly::skip_molly_data(fin,spectrum.head);
			}
		    }else{
			std::cerr << "Spectrum " << count << " of " << file[j] 
				  << " has no positional and/or timing information." << std::endl;
			if(warn) std::cerr << spectrum.head << std::endl;
			Colly::skip_molly_data(fin,spectrum.head);
		    }
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
	}
	fout.close();
    }
    catch(const Subs::Position::Position_Error& err){
	std::cerr << err << std::endl;
    }
    catch(const std::string& error){
	std::cerr << error << std::endl;
    }
}








