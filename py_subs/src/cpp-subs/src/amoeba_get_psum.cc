#include <vector>
#include "trm/subs.h"

void Subs::amoeba_get_psum(const std::vector<std::pair<std::vector<double>, double> >& params, std::vector<double>& psum){
  for(int j=0; j<psum.size(); j++){
    double sum = 0.;
    for(int i=0; i<params.size(); i++)
      sum += params[i].first[j];
    psum[j] = sum;
  }
}



