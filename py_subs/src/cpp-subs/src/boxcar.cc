#include "trm/subs.h"

/* Boxcar filter routine (Stu Littlefair)
 * \param data the data to be filtered
 * \param filt the output filtered data
 * \param width the width in bins of the filter (odd)
 */
void Subs::boxcar(const Buffer1D<double>& data, Buffer1D<double>& filt, int width){

  if(width == 1){
    filt = data;
    return;
  }

  filt.resize(data.size());
  int np     = data.size();

  for(int j=0; j<np; j++){
    double sum=0.;
    int npix=0;
    for(int i=std::max(j-width,0); i<std::min(j+width+1,np); i++){
      sum += data[i];
      npix++;
    }
    filt[j]=sum/npix;
  }

}


