/*

!!begin
!!author  T.R. Marsh
!!created 25 June 2001
!!revised 07 Jan 2008
!!descr   reorders molly spectra by position, time or phase
!!root    reorder
!!index   reorder
!!class   Programs
!!css     style.css
!!head1   reorder - reorders molly spectra by position, time or phase

!!emph{reorder} reads in molly files and then reorders spectra
by their position (RA, Dec), time (RJD) or orbital phase (Orbital phase). 
For position it does so on a key 
composed of RA then Dec precessed to a common equinox. If identical 
keys are encountered it preserves the order that it finds the spectra 
in the file on output. All spectra are loaded prior to output so it is 
possible to overwrite the input file safely. If values are not located, then
they default to zero and will appear as a block in the output spectra, i.e.
no spectrum will be lost.

!!head2 Invocation

reorder infile outfile key (clobber)

!!head2 Arguments

!!table
!!arg{infile}{Input molly file}
!!arg{outfile}{Output molly file}
!!arg{key}{'p', 't' or 'o' for position, time or orbital phase}
!!arg{clobber}{Enter "true" to allow existing files to be overwritten
(including the input file). Off by default}
!!table

!!end

*/

#include <cstdlib>
#include <iomanip>
#include <fstream>
#include <vector>
#include <map>

#include "trm/subs.h"
#include "trm/input.h"
#include "trm/header.h"
#include "trm/position.h"
#include "trm/colly.h"

int main(int argc, char* argv[]){

    const double CEPOCH = 2000.;
    Subs::Date::print_method = 2;

    std::cout << "PRECESSION BROKEN" << std::endl;

    try{
    
	// Construct Input object
	Subs::Input input(argc, argv, Colly::COLLY_ENV, Colly::COLLY_DIR);

	// Define inputs
	input.sign_in("infile",   Subs::Input::LOCAL,  Subs::Input::PROMPT);
	input.sign_in("outfile",  Subs::Input::LOCAL,  Subs::Input::PROMPT);
	input.sign_in("key",      Subs::Input::LOCAL,  Subs::Input::PROMPT);
	input.sign_in("clobber",  Subs::Input::LOCAL,  Subs::Input::NOPROMPT);

	std::string infile;
	input.get_value("infile", infile, "infile", "input molly file");
	std::string outfile;
	input.get_value("outfile", outfile, "outfile", "output reordered molly file");
	char key;
	input.get_value("key", key, 'p', "ptoPTO", "sort key p(osition), t(ime), o(rbital phase)");
	key = toupper(key);
	bool clobber;
	input.get_value("clobber", clobber, false, "allow output to ovwerwrite existing files?");

	Subs::Header::Hnode *rap, *decp = NULL, *eqp = NULL, *tp = NULL, *op;

	std::multimap<std::string,Colly::molly> specs;
	std::multimap<double,Colly::molly> specd;
	Subs::Position pos;
	Colly::molly spectrum;
	std::string skey;
	double dkey = 0;

	std::ifstream fin(infile.c_str(), std::ios::in | std::ios::binary);
	if(!fin)
	    throw std::string("Could not open input file = " + infile);

	size_t count = 0;
	try{
	    while(Colly::read_molly_head(fin,spectrum.head,spectrum.original)){
		count++;   
		if(key == 'P'){
		    if(
			(rap    = spectrum.head.find("RA"))->has_data()  &&      
			(decp   = spectrum.head.find("Dec"))->has_data() &&     
			(eqp    = spectrum.head.find("Equinox"))->has_data()){ 
	      
			pos  = Subs::Position(rap->value->get_double(),
					      decp->value->get_double(),
					      eqp->value->get_double());
//			pos.precess(CEPOCH);	    
		    }else{
			pos = Subs::Position(0., 0., CEPOCH);
		    }
		    skey = pos.ra_dec();
		}else if(key == 'T'){
		    if((tp    = spectrum.head.find("RJD"))->has_data()){
			dkey = tp->value->get_double();
		    }else{
			dkey = 0.;
		    }
		}else if(key == 'O'){
		    if((op = spectrum.head.find("Orbital phase"))->has_data()){
			dkey = tp->value->get_double();
		    }else{
			dkey = 0.;
		    }
		} 
	
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
		if(key == 'P'){
		    specs.insert(std::make_pair(skey,spectrum));
		}else if(key == 'T' || key == 'O'){
		    specd.insert(std::make_pair(dkey,spectrum));
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
	fin.close();
    
	if(count){
	    std::cout << "Finished file input. " << count
		      << " spectra read." << std::endl;
	}else{
	    std::cerr << "No spectra read!" << std::endl;
	    abort();
	}
    
	// Now output in order
    
	std::ofstream fout;
	if(clobber){
	    fout.open(outfile.c_str(), std::ios::out | std::ios::binary);
	}else{
	    std::ifstream ftest(outfile.c_str(), std::ios::in | std::ios::binary);
	    if(ftest)
		throw std::string("File = " + outfile + "  already exists.");
	    fout.open(outfile.c_str(), std::ios::out | std::ios::binary);
	}
	if(!fout)
	    throw std::string("Could not open output file = " + outfile);
    
	Subs::Header head;
	if(key == 'P'){
	    typedef std::multimap<std::string,Colly::molly>::const_iterator IS;
	    for(IS j=specs.begin(); j!= specs.end(); j++){
		head = j->second.head;
		Colly::write_molly_head(fout,head,j->second.original);
		if(!(fout.write(j->second.buff,j->second.nbytes)))
		    throw std::string("Failed during write of molly data");
	    }
	}else if(key == 'T' || key == 'O'){
	    typedef std::multimap<double,Colly::molly>::const_iterator ID;
	    for(ID j=specd.begin(); j!= specd.end(); j++){
		head = j->second.head;
		Colly::write_molly_head(fout,head,j->second.original);
		if(!(fout.write(j->second.buff,j->second.nbytes)))
		    throw std::string("Failed during write of molly data");
	    }
	}
	fout.close();
	std::cout << "Successfully written file = " << outfile << std::endl;
    }
    catch(const std::string& error){
	std::cerr << error << std::endl;
    }
}














