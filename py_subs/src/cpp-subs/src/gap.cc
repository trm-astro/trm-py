/*

!!sphinx

gap -- computes gap probability
===============================

This carries out Monte Carlo runs to work out how often one expects to see a
gap of a given size in a uniform distribution

Arguments:


 width : width of gap, 0.00001 to 0.99999
 npoints : number of values, 2 or more
 nmonte : number of MC trials
 seed : seed integer

!!sphinx

*/

#include <climits>
#include <string>
#include "trm/subs.h"
#include "trm/input.h"

int main(int argc, char* argv[]){
  try{

    // Construct Input object
    Subs::Input input(argc, argv, "SUBS_ENV", ".subs");

    // Define inputs
    input.sign_in("width",   Subs::Input::LOCAL, Subs::Input::PROMPT);
    input.sign_in("npoints", Subs::Input::LOCAL, Subs::Input::PROMPT);
    input.sign_in("nmonte",  Subs::Input::LOCAL, Subs::Input::PROMPT);
    input.sign_in("seed",    Subs::Input::LOCAL, Subs::Input::PROMPT);

    float width;
    input.get_value("width", width, 0.1f, 1.e-5f, 0.9999999f, "width of gap (relative to total width)");
    int npoints;
    input.get_value("npoints", npoints, 40, 2, 1000000, "number of points per monte carlo run");
    int nmonte;
    input.get_value("nmonte", nmonte, 10000, 1, 1000000, "number of monte carlo runs");
    int seed;
    input.get_value("seed", seed, 78932, INT_MIN, INT_MAX, "seed integer for monte carlo runs");
    if(seed > 0) seed = -seed;

    float *vals = new float[npoints], gap;
    int naslarge = 0;
    for(int i=0; i<nmonte; i++){

      for(int j=0; j<npoints; j++)
	vals[j] = Subs::ran2(seed);

      Subs::quicksort(vals, npoints);

      gap = -1.;
      for(int j=0; j<npoints-1; j++)
	if(vals[j+1] - vals[j] > gap)
	  gap = vals[j+1] - vals[j];

      if(gap >= width) naslarge++;
    }

    std::cout << "Largest gap as large as " << width << " on " << 100.*naslarge/nmonte << "% of the " 
	 << nmonte << " trials." << std::endl;
  }
  catch(const std::string& message){
    std::cerr << message << std::endl;
  }
}
