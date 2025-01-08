#include <cmath>
#include <iostream>
#include <string>
#include "trm/subs.h"
#include "trm/buffer2d.h"

/**
 * The routine implements Gauss-Jordan elimination to solve linear equations of form
 * Ax = b where A is a square matrix, and b is a set of vectors in a 
 * 2D array.
 * \param a n by n matrix
 * \param b n by m matrix, i.e. m vectors of length n
 * \return On return a will be its inverse and b will contain the m solution vectors x.
 * \exception Throws Subs::Subs_Error exceptions if there are problems.
 */

void Subs::gaussj(Buffer2D<double>& a, Buffer2D<double>& b){

  size_t n = a.nrow();
  if(!n) 
    throw Subs_Error("void Subs::gaussj(Buffer2D<double>&, Buffer2D<double>&): matrix a is null");
  if(int(n) != a.ncol()) 
    throw Subs_Error("void Subs::gaussj(Buffer2D<double>&, Buffer2D<double>&): matrix a is not square");
  if(int(n) != b.nrow()) 
    throw Subs_Error("void Subs::gaussj(Buffer2D<double>&, Buffer2D<double>&): matrices a and b not compatible");
  size_t m = b.ncol();

  Buffer1D<size_t> indxc(n), indxr(n), ipiv(n);
  size_t i, icol = 0, irow = 0, j, k, l , ll;
  double big, dum, pivinv;
  
  for(j=0;j<n;j++) ipiv[j] = 0;
  for(i=0;i<n;i++){
    big = 0.;
    for(j=0;j<n;j++)
      if(ipiv[j] != 1)
	for(k=0;k<n;k++){
	  if(ipiv[k] == 0){
	    if(fabs(a[j][k]) >= big){
	      big = fabs(a[j][k]);
	      irow = j;
	      icol = k;
	    }
	  }else if (ipiv[k] > 1) 
	    throw Subs_Error("void Subs::gaussj(Buffer2D<double>&, Buffer2D<double>&): singular matrix 1");
	}
    ++(ipiv[icol]);

    if(irow != icol){
      for(l=0;l<n;l++) std::swap(a[irow][l],a[icol][l]);
      for(l=0;l<m;l++) std::swap(b[irow][l],b[icol][l]);
    }
    indxr[i] = irow;
    indxc[i] = icol;
    if(a[icol][icol] == 0.) 
      throw Subs_Error("void Subs::gaussj(Buffer2D<double>&, Buffer2D<double>&): singular matrix 2");

    pivinv = 1.0/a[icol][icol];
    a[icol][icol] = 1.0;
    for(l=0;l<n;l++) a[icol][l] *= pivinv;
    for(l=0;l<m;l++) b[icol][l] *= pivinv;

    for(ll=0;ll<n;ll++)
      if(ll != icol){
	dum = a[ll][icol];
	a[ll][icol] = 0.;
	for(l=0;l<n;l++) a[ll][l] -= a[icol][l]*dum;
	for(l=0;l<m;l++) b[ll][l] -= b[icol][l]*dum;
      }
  }

  // take care not to try to reach l = -1 given data type
  for(l=n;l>0;l--){
    if(indxr[l-1] != indxc[l-1])
      for(k=0;k<n;k++)
	std::swap(a[k][indxr[l-1]],a[k][indxc[l-1]]);
  }
}

/**
 * The routine implements Gauss-Jordan elimination to solve linear equations of form
 * Ax = b where A is a square matrix, and b is a set of vectors in a 2D array, but allows
 * one to override the usual dimensions in the sense of restricting them to a smaller value.
 * \param nover overrides the size of matrices
 * \param a m by m matrix, with m >= n
 * \param b m by p matrix, i.e. p vectors of length m >= n
 * \return On return a will be its inverse and b will contain the m solution vectors x.
 * \exception Throws Subs::Subs_Error exceptions if there are problems.
 */

void Subs::gaussj(int nover, Buffer2D<double>& a, Buffer2D<double>& b){

  if(nover > a.nrow() || nover > a.ncol())
    throw Subs_Error("void Subs::gaussj(int, Buffer2D<double>&, Buffer2D<double>&): nover > a.nrow() or a.ncol()");
  if(nover > b.nrow())
    throw Subs_Error("void Subs::gaussj(int, Buffer2D<double>&, Buffer2D<double>&): nover > b.nrow()");

  size_t m = b.ncol();

  Buffer1D<size_t> indxc(nover), indxr(nover), ipiv(nover);
  size_t i, icol = 0, irow = 0, j, k, l , ll;
  double big, dum, pivinv;
  size_t nnover = size_t(nover);
  for(j=0;j<nnover;j++) ipiv[j] = 0;
  for(i=0;i<nnover;i++){
    big = 0.;
    for(j=0;j<nnover;j++)
      if(ipiv[j] != 1)
	for(k=0;k<nnover;k++){
	  if(ipiv[k] == 0){
	    if(fabs(a[j][k]) >= big){
	      big = fabs(a[j][k]);
	      irow = j;
	      icol = k;
	    }
	  }else if (ipiv[k] > 1) 
	    throw Subs_Error("void Subs::gaussj(int, Buffer2D<double>&, Buffer2D<double>&): singular matrix 1");
	}
    ++(ipiv[icol]);

    if(irow != icol){
      for(l=0;l<nnover;l++) std::swap(a[irow][l],a[icol][l]);
      for(l=0;l<m;l++) std::swap(b[irow][l],b[icol][l]);
    }
    indxr[i] = irow;
    indxc[i] = icol;
    if(a[icol][icol] == 0.) 
      throw Subs_Error("void Subs::gaussj(int, Buffer2D<double>&, Buffer2D<double>&): singular matrix 2");

    pivinv = 1.0/a[icol][icol];
    a[icol][icol] = 1.0;
    for(l=0;l<nnover;l++) a[icol][l] *= pivinv;
    for(l=0;l<m;l++) b[icol][l] *= pivinv;

    for(ll=0;ll<nnover;ll++)
      if(ll != icol){
	dum = a[ll][icol];
	a[ll][icol] = 0.;
	for(l=0;l<nnover;l++) a[ll][l] -= a[icol][l]*dum;
	for(l=0;l<m;l++) b[ll][l] -= b[icol][l]*dum;
      }
  }

  // take care not to try to reach l = -1 given data type
  for(l=nnover;l>0;l--){
    if(indxr[l-1] != indxc[l-1])
      for(k=0;k<nnover;k++)
	std::swap(a[k][indxr[l-1]],a[k][indxc[l-1]]);
  }
}

/**
 * The routine implements Gauss-Jordan elimination to solve linear equations of form
 * Ax = b where A is a square matrix, and b is a set of vectors in a 
 * 2D array.
 * \param a n by n matrix
 * \param n dimension of matrix a
 * \param b n by m matrix, i.e. m vectors of length n
 * \param m second dimension of matrix b
 * \return On return a will be its inverse and b will contain the m solution vectors x.
 * \exception Throws Subs::Subs_Error exceptions if there are problems.
 */

void Subs::gaussj(double** a, int n, double** b, int m){

  if(n < 1) 
    throw Subs_Error("void Subs::gaussj(double**, int, double**, int): n < 1");

  if(m < 1) 
    throw Subs_Error("void Subs::gaussj(double**, int, double**, int): m < 1");

  Buffer1D<size_t> indxc(n), indxr(n), ipiv(n);
  size_t i, icol = 0, irow = 0, j, k, l , ll;
  double big, dum, pivinv;
  size_t nn = size_t(n), mm = size_t(m);
  
  for(j=0;j<nn;j++) ipiv[j] = 0;
  for(i=0;i<nn;i++){
    big = 0.;
    for(j=0;j<nn;j++)
      if(ipiv[j] != 1)
	for(k=0;k<nn;k++){
	  if(ipiv[k] == 0){
	    if(fabs(a[j][k]) >= big){
	      big = fabs(a[j][k]);
	      irow = j;
	      icol = k;
	    }
	  }else if (ipiv[k] > 1) 
	    throw Subs_Error("void Subs::gaussj(double**, int, double**, int): singular matrix 1");
	}
    ++(ipiv[icol]);

    if(irow != icol){
      for(l=0;l<nn;l++) std::swap(a[irow][l],a[icol][l]);
      for(l=0;l<mm;l++) std::swap(b[irow][l],b[icol][l]);
    }
    indxr[i] = irow;
    indxc[i] = icol;
    if(a[icol][icol] == 0.) 
      throw Subs_Error("void Subs::gaussj(double**, int, double**, int): singular matrix 2");

    pivinv = 1.0/a[icol][icol];
    a[icol][icol] = 1.0;
    for(l=0;l<nn;l++) a[icol][l] *= pivinv;
    for(l=0;l<mm;l++) b[icol][l] *= pivinv;

    for(ll=0;ll<nn;ll++)
      if(ll != icol){
	dum = a[ll][icol];
	a[ll][icol] = 0.;
	for(l=0;l<nn;l++) a[ll][l] -= a[icol][l]*dum;
	for(l=0;l<mm;l++) b[ll][l] -= b[icol][l]*dum;
      }
  }

  // take care not to try to reach l = -1 given data type
  for(l=nn;l>0;l--){
    if(indxr[l-1] != indxc[l-1])
      for(k=0;k<nn;k++)
	std::swap(a[k][indxr[l-1]],a[k][indxc[l-1]]);
  }
}
