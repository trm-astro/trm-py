#include <iostream>
#include "trm/subs.h"
#include "trm/constants.h"
#include "cpgplot.h"

int main(){

  const int NDATA=100;
  float x[NDATA], y[NDATA], var[NDATA];
  const float FWHM = 10;
  for(int i=0; i<NDATA; i++)
    x[i] = i+1;

  double pos;
  float epos;
  long int seed = -127897;
  for(int j=0; j<100000; j++){ 
    for(int i=0; i<NDATA; i++){	
      var[i]  = 1;
      y[i]    = exp(-Subs::sqr((i-55.)/5.)/2.) + sqrt(var[i])*Subs::gauss2(seed);
    }
    try{
      Subs::centroid(y, var, 0, NDATA-1, FWHM, 50., true, pos, epos);
      std::cout << "trial " << j+1 << ", position = " << pos << std::endl; 
    }
    catch(const Subs::Subs_Error& err){
      std::cerr << "Error message = " << err << std::endl;
      std::cerr << "j = " << j << std::endl;
    } 
  }
}
