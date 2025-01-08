#include <cmath>
#include <iostream>
#include "trm/subs.h"
#include "trm/buffer2d.h"

namespace Jacobi{

  inline void rotate(Subs::Buffer2D<double>& a, int i, int j, int k, int l, double s, double tau){
    double g, h;    
    g       = a[i][j];
    h       = a[k][l];
    a[i][j] = g-s*(h+g*tau);
    a[k][l] = h+s*(g-h*tau);
    return;
  }
};

/**
 * This function computes the eigenvalues and vectors of a real, symmetric matrix.
 * \param a n by n matrix. On output elements above the diagonal are destroyed. i.e elements like a[0][1] are changed
 * while a[0][0] and a[1][0] are unchanged.
 * \param d the eigenvalues
 * \param v n by n eigenvectors of a in column form. 
 * \param nrot the number of Jacobi rotations needed (returned).
 * \exception Subs::Subs_Error exceptions are thrown.
 */

void Subs::jacob(Buffer2D<double>& a, Buffer1D<double>& d, Buffer2D<double>& v, int &nrot){
  
  int n = a.nrow();
  if(!n) 
    throw Subs_Error("Subs::jacob(Buffer2D<double>&, Buffer1D<double>&, Buffer2D<double>&, int&): null matrix");
  if(n != a.ncol()) 
    throw Subs_Error("Subs::jacob(Buffer2D<double>&, Buffer1D<double>&, Buffer2D<double>&, int&): matrix not square");
  if(n != d.size()) 
    throw Subs_Error("Subs::jacob(Buffer2D<double>&, Buffer1D<double>&, Buffer2D<double>&, int&): matrix a and std::vector b not compatible");
  if(n != v.ncol() || n != v.nrow()) 
    throw Subs_Error("Subs::jacob(Buffer2D<double>&, Buffer1D<double>&, Buffer2D<double>&, int&): matrices a and v not compatible");

  int j, iq, ip, i;
  double tresh, theta, tau, t, sm, s, h, g, c;

  Buffer1D<double> b(n), z(n);

  // Initialise to identity matrix

  for(ip=0; ip<n; ip++){
    for(iq=0; iq<n; iq++) v[ip][iq] = 0.;
    v[ip][ip] = 1.;
  }

  // Initialise b and d to diagonal of a

  for(ip=0; ip<n; ip++){
    b[ip] = d[ip] = a[ip][ip];
    z[ip] = 0.;
  }
      
  nrot = 0;
  for(i=1; i<=50; i++){
    sm = 0.;

    // Sum off diagonal elements

    for(ip = 0; ip < n-1; ip++){ 
      for(iq=ip+1; iq<n; iq++) 
	sm += fabs(a[ip][iq]);
    }

    // Normal return

    if(sm == 0.) return;
    if(i < 4) 
      tresh = 0.2*sm/(n*n); 
    else
      tresh = 0.;
    for(ip=0; ip< n-1; ip++){
      for(iq=ip+1; iq<n; iq++){
	g=100.*fabs(a[ip][iq]);

	// After four sweeps skip rotation if off-diag elements is small

	if(i > 4 && (fabs(d[ip])+g) == fabs(d[ip]) &&
	   (fabs(d[iq])+g) == fabs(d[iq]))
	  a[ip][iq]= 0.;
	else if(fabs(a[ip][iq]) > tresh) {
	  h = d[iq] - d[ip];
	  if((fabs(h)+g) == fabs(h))
	    t = a[ip][iq]/h;
	  else {
	    theta = 0.5*h/a[ip][iq];
	    t = 1./(fabs(theta)+sqrt(1.+theta*theta));
	    if(theta < 0.) t = -t;
	  }
	  c = 1./sqrt(1.+t*t);
	  s = t*c;
	  tau= s/(1.+c);
	  h = t*a[ip][iq];
	  z[ip] -= h;
	  z[iq] += h;
	  d[ip] -= h;
	  d[iq] += h;
	  a[ip][iq] = 0.;
	  for(j=0;j<=ip-1;j++){
	    Jacobi::rotate(a,j,ip,j,iq,s,tau);
	  }
	  for(j=ip+1;j<=iq-1;j++){
	    Jacobi::rotate(a,ip,j,j,iq,s,tau);
	  }
	  for(j=iq+1;j<n;j++){
	    Jacobi::rotate(a,ip,j,iq,j,s,tau);
	  }
	  for(j=0;j<n;j++){
	    Jacobi::rotate(v,j,ip,j,iq,s,tau);
	  }
	  ++nrot;
	}
      }
    }
    for(ip=0; ip<n; ip++){
      b[ip] += z[ip];
      d[ip] = b[ip];
      z[ip] = 0.;
    }
  }
  throw Subs_Error("Subs::jacob(Buffer2D<double>&, Buffer1D<double>&, Buffer2D<double>&, int&): too many iterations in jacob");
}




	  



