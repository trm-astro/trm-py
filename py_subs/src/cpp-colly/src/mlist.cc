/*

!!begin
!!author  T.R. Marsh
!!created 18 June 2001
!!revised 09 July 2007
!!descr   lists molly headers
!!root    mlist
!!index   mlist
!!class   Programs
!!css     style.css
!!head1   mlist - lists molly headers.

!!emph{mlist} reads in molly files and lists their headers.
Several levels  of listing are possible.

!!head2 Invocation

mlist pmode file1 file2 file3 ...

!!head2 Arguments

!!table
!!arg{pmode}{Print level for each item. 0 just prints the name,
1 prints the name and value, 2 prints the name, value and type,
3 prints the name, value, type and comment. 1 is usually the best
choice.}
!!arg{file1 etc}{Molly file names. These can be *.mol for example}
!!table

!!end

*/

#include <stdlib.h>
#include <fstream>

#include "trm/subs.h"
#include "trm/input.h"
#include "trm/header.h"
#include "trm/colly.h"

int main(int argc, char* argv[]){

    try{
    
	// Construct Input object
	Subs::Input input(argc, argv, Colly::COLLY_ENV, Colly::COLLY_DIR);

	// Define inputs
	input.sign_in("pmode",    Subs::Input::LOCAL,  Subs::Input::PROMPT);
	input.sign_in("flist",    Subs::Input::GLOBAL, Subs::Input::PROMPT);

	int pmode;
	input.get_value("pmode",  pmode, 1, 0, 3, "print mode (0,1,2 or 3)");
	std::string flist;
	input.get_value("flist",  flist, "list", "list of molly files");

	Subs::Hitem::set_default_output();
	Subs::Header head; 
	std::vector<std::string> original;

	// Read file names
	std::ifstream liststr(flist.c_str());
	if(!liststr)
	    throw std::string("Failed to open ") + flist;    
	std::string fname;
	while(liststr >> fname){
	    std::ifstream fin(fname.c_str(), std::ios::in | std::ios::binary);
	    if(fin){
		int count = 0;
		while(Colly::read_molly_head(fin,head,original)){
		    Colly::rejig_molly_head(head);
		    try{
			std::cout << "\nFile: " << fname << ", spectrum = " << ++count << std::endl;
			std::cout << head << std::endl;
		    }
		    catch(const std::string& message){
			std::cerr << message << std::endl;
		    }
		    Colly::skip_molly_data(fin,head);
		}
		fin.close();
	    }else{
		std::cerr << "Skipped " << fname << " which could not be opened for input." << std::endl;
	    }
	}
    }
    catch(const Colly::Colly_Error& error){
	std::cerr << error << std::endl;
    }
    catch(const std::string& message){
	std::cerr << message << std::endl;
    }
}





