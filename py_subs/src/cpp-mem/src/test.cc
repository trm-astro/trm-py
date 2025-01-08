
#include "../base/subs.h"
#include "../base/graphics.h"
#include "trm/memsys.h"

const int NDATA = 200;
float dat[NDATA];

void op(float model[], float data[] , const int npt);
void tr(float data[] , float model[], const int npt);

int main(){

  const int MXPIX = 100000;
  float buff[MXPIX];
  Mem::Gbl::st = buff;
  
  Mem::memcore(MXPIX,NDATA,NDATA);

  float model[NDATA], data[NDATA], x[NDATA], e[NDATA], image[NDATA];
  for(int i=0; i<NDATA; i++){
    x[i] = i+1;
    model[i] = 1.;
  }
  for(int i=0; i<NDATA; i++){
    model[i] += 2.*exp(-0.2*sqr(i-56.));
    model[i] += 4.*exp(-0.5*sqr(i-106.));
    model[i] += 4.*exp(-0.1*sqr(i-146.));
  }

  op(model,data,NDATA);
  long int seed = -8269;
  float sigma = 0.01;
  for(int i=0; i<NDATA; i++){
    data[i] += sigma*gauss2(seed);
    e[i]     = sigma;
    image[i] = 0.5;
  }

  // move stuff into vectors

  for(int i=0; i<NDATA; i++){
    Mem::Gbl::st[Mem::Gbl::kb[0]+i]  = image[i];
    Mem::Gbl::st[Mem::Gbl::kb[20]+i] = data[i];
    Mem::Gbl::st[Mem::Gbl::kb[21]+i] = 2./(e[i]*e[i]*NDATA);
  }
  
  float caim = 1., rmax = 0.2, c,test,acc,cnew,s,rnew,snew,sumf;
  for(int it=0; it<20; it++){
    std::cout << "\nIteration " << it+1 << std::endl;
    Mem::memprm(10,20,caim,rmax,1.,acc,c,test,cnew,s,rnew,snew,sumf);
  }

  // move stuff out for plot

  for(int i=0; i<NDATA; i++){
    image[i] = Mem::Gbl::st[Mem::Gbl::kb[0]+i];
  }

  try{
    graphics plot("/xs");
    plot.sci(4);
    plot.scf(2);
    plot.sch(1.5);
    plot.env(0.,float(NDATA+1),0.,6.,0,0);
    plot.sci(2);
    plot.lab("X","Y");
    plot.sci(3);
    plot.line(NDATA,x,model);
    plot.sci(2);
    plot.line(NDATA,x,image);
    plot.sci(1);
    plot.bin(NDATA,x,data,true);
    plot.clos();
  }

  catch(graphics::no_device o){
    std::cerr << "graphics error. No device open";
    std::string s = o.get_mess();
    if(s != ""){
      std::cerr << " in " << s << "." << std::endl;
    }else{
      std::cerr << "." << std::endl;
    }
  }

  catch(graphics::graphics_error o){
    std::cerr << "graphics error.";
    std::string s = o.get_mess();
    if(s != ""){
      std::cerr << s << std::endl;
    }else{
      std::cerr << std::endl;
    }
  }

}

void Mem::opus(const int j, const int k){

  std::cout << "    OPUS " << j+1 << " ---> " << k+1 << std::endl;
  op(Mem::Gbl::st+Mem::Gbl::kb[j],
     Mem::Gbl::st+Mem::Gbl::kb[k], NDATA);
}

const int NSP = 15;
float wgt[NSP] = {0.5,1.,1.5,2.,3.,4.,5.,6.,5.,4.,3.,2.,1.5,1.,0.5};

void op(float model[], float data[], const int npt){ 
  int k;
  float sum, wnorm = 0.;

  for(int i=0; i<NSP; i++) 
    wnorm += wgt[i];
  for(int i=0; i<npt; i++){
    sum = 0.;
    k = i - NSP/2;
    for(int j=0; j<NSP; j++, k++){
      if(k<0){
	sum += wgt[j]*model[0];
      }else if(k >= npt){
	sum += wgt[j]*model[npt-1];
      }else{
	sum += wgt[j]*model[k];
      }
    }
    data[i] = sum/wnorm;
  }
}


void tr(float data[], float model[], const int npt){

  int k;
  float wnorm = 0.;

  for(int i=0; i<NSP; i++)
    wnorm += wgt[i];

  for(int i=0; i<npt; i++)
    model[i] = 0.;

  for(int i=0; i<npt; i++){
    k = i - NSP/2;
    for(int j=0; j<NSP; j++, k++){
      if(k<0){
	model[i] += wgt[j]*data[0];
      }else if(k >= npt){
	model[i] += wgt[j]*data[npt-1];
      }else{
	model[i] += wgt[j]*data[k];
      }
    }
    model[i] /= wnorm;
  }
}


