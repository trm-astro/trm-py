// test program to compate normal equations with singular value decomposition
#include "trm/subs.h"
#include "trm/array1d.h"
#include "trm/buffer2d.h"
#include "trm/plot.h"
#include "cpgplot.h"

int main(){

  const int NPOLY = 10;
  const int NDATA = 2000;

  Subs::Array1D<double> x(NDATA), y(NDATA);
  Subs::Array1D<float>  e(NDATA);
  float  *ed = new float[NDATA];
  double *yd = new double[NDATA];

  int seed = -363787;
  
  std::cout << "NPOLY = " << NPOLY << std::endl;

  // generate test data with a ramp
  float x1 = 100, x2 = 1050, step = 4;
  for(int i=0; i<NDATA; i++){

    x[i] = i + 1;
    y[i] = 1 + 0.001*i + 0.11*Subs::gauss2(seed);
    if(x[i] > 500 && x[i] < 550)
      y[i] += sin(6.28*(x[i]-500)/50);
    ed[i] = e[i] = 1.;
    yd[i] = y[i];
  }

  Subs::Buffer2D<double> a(NDATA, NPOLY), v;
  Subs::Array1D<double> w, b(NDATA);

  double **func = new double*[NDATA];
  for(int i=0; i<NDATA; i++)
    func[i] = new double[NPOLY];

  // Fill in matrix a
  for(int i=0; i<NDATA; i++){
    double val = 1/Subs::sqr(e[i]), xtrans = (x[i]-(x[0]+x[NDATA-1])/2.)/((x[NDATA-1]-x[0])/2);
    b[i] = val*y[i];
    func[i][0] = a[i][0] = 1;
    func[i][1] = a[i][1] = xtrans;
    int j;
    for(j=2; j<NPOLY; j++){
      func[i][j] = a[i][j] = 2*xtrans*a[i][j-1] - a[i][j-2];
    } 
  }

  Subs::Plot plot("/xs");

  // plot data
  cpgenv(-10.,NDATA+50, -5, 10, 0, 0);
  pgpt(x,y,5);

  Subs::Array1D<float> fit(NDATA);

  // usual normal equations approach
  double *coeff = new double[NPOLY], **covar = new double*[NPOLY];
  for(int i=0; i<NPOLY; i++)
    covar[i] = new double[NPOLY];

  Subs::llsqr(NDATA, yd, ed, NPOLY, func, coeff, covar);

  for(int j=0; j<NPOLY; j++)
    std::cout << j << " coeff " << coeff[j] << std::endl;
  double chisq = 0;
  for(int i=0; i<NDATA; i++){
    double sum = 0;
    for(int j=0; j<NPOLY; j++)
      sum += coeff[j]*func[i][j];
    fit[i] = sum;
    chisq += Subs::sqr(yd[i]-fit[i]);
  }
  std::cout << "first chi**2 = " << chisq << std::endl;
  cpgsci(2);
  pgline(x,fit);

  // SVD it
  Subs::svdcmp(a, w, v); 

  // Edit out feeble bits
  double wmax = w.max();
  int ntrim = 0;
  for(int i=0; i<w.size(); i++){
    if(w[i] < 1.e-14*wmax){
      w[i] = 0;
      ntrim++;
    }
  }

  std::cout << "trimmed " << ntrim << std::endl;

  Subs::Buffer1D<double> solution;
  Subs::svbksb(a, w, v, b, solution); 

  // Plot the safer fit
  
  chisq = 0;
  for(int j=0; j<NDATA; j++){
    double val = 1/Subs::sqr(e[j]), fac = (x[j]-(x[0]+x[NDATA-1])/2.)/((x[NDATA-1]-x[0])/2);
    int k;
    double sum = 0;
    for(k=0; k<NPOLY; k++){
      sum += solution[k]*func[j][k];
    }
    fit[j] = sum;
    chisq += Subs::sqr(yd[j]-fit[j]);
  }
  std::cout << "second chi**2 = " << chisq << std::endl;
  cpgsci(3);
  pgline(x,fit);

  delete[] ed;
  delete[] yd;
}
