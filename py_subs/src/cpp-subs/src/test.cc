#include "trm/header.h"
#include "trm/input.h"

int main(int argc, char* argv[]){
   
    try {

	Subs::Header head;
	head.set("Xtra",              new Subs::Hdirectory("Molly data"));
	head.set("Xtra.FCODE",        new Subs::Hint(2, "Molly format code"));
	head.set("Xtra.UNITS",        new Subs::Hstring("COUNTS          ", "Units of fluxes"));
	head.set("Xtra.NARC",         new Subs::Hint(0,"Number of arc coefficients"));
	head.set("Record",            new Subs::Hint(23, "Record number"));
	head.set("Exposure",          new Subs::Hfloat(20.5f, "Exposure time (sec)"));
	head.set("Object",            new Subs::Hstring("IP Peg", "target name"));
	head.set("Timing",            new Subs::Hdirectory("Timing information"));
	head.set("Timing.deadtime",   new Subs::Hfloat(5.5f, "dead time (sec)"));
	head.set("Timing.frametime",  new Subs::Hfloat(25.5f, "frame time (msec)"));
	head.set("Timing.CCD",        new Subs::Hdirectory("CCD specific"));
	head.set("Timing.CCD.xleft",  new Subs::Hint(10, "Left X pixel"));
	head.set("Timing.CCD.xright", new Subs::Hint(610, "Right X pixel"));
	head.set("Timing.CCD.ystart", new Subs::Hint(13, "Lower Y pixel"));
	head.set("Timing.CCD.nx",     new Subs::Hint(100, "X dimension"));
	head.set("Timing.CCD.ny",     new Subs::Hint(100, "Y dimension"));
	head.set("Window",            new Subs::Hdirectory("Window information"));
	head.set("Window.size",       new Subs::Hint(1000, "number of pixels"));


	Subs::Header::start_string = "# ";
//	std::cout << "\n\nhead = \n" << head;

//	head.move_to_top("Exposure");
//	head.move_to_top("Timing.CCD.ystart");

	std::cout << "\n\nhead = \n" << head;

//	std::cout << head["Exposure"]->get_float() << std::endl;

	head.erase("Xtra");

	std::cout << "\n\nhead = \n" << head;

    }
    catch(const std::string& err){
	std::cerr << err << std::endl;
    }
}

