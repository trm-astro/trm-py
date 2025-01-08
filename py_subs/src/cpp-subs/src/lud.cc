#include <cmath>
#include <iostream>
#include <string>
#include "trm/subs.h"
#include "trm/buffer2d.h"

/**
 * LU decomposition is a useful way to invert matrices and solve linear
 * equations. ludcmp is the main decomposition routine. The matrix \c a is modified
 * to the LU decomposition of a row-wise permutation of itself. \c indx
 * records the permutation and \c d is returned as +/- 1 depending upon
 * the odd/even nature of the row interchanges.
 * \param a    the square matrix to be decomposed
 * \param n    the dimension of the matrix
 * \param indx index array to record the permutations made
 * \param d    +/- 1 depending on sign of permutation
 * \exception Throws a Subs::Subs_Error if the matrix is singular
 */

void Subs::ludcmp(double **a, unsigned int n, unsigned *indx, double &d){
  unsigned int i, imax=0, j, k;
  double big, dum, sum, temp;
  double *vv;
  const double TINY = 1.e-20;

  vv = new double [n];
  d  = 1.;
  for(i=0;i<n;i++){
    big = 0.;
    for(j=0;j<n;j++) 
      if((temp=fabs(a[i][j])) > big) big = temp;
    if(big == 0.)
      throw Subs_Error("void Subs::ludcmp(double, unsigned int, unsigned*, double&): singular matrix");
    vv[i] = 1./big;
  }
  for(j=0;j<n;j++){
    for(i=0;i<j;i++){
      sum = a[i][j];
      for(k=0;k<i;k++) sum -= a[i][k]*a[k][j];
      a[i][j] = sum;
    }
    big = 0.;
    for(i=j;i<n;i++){
      sum = a[i][j];
      for(k=0;k<j;k++) sum -= a[i][k]*a[k][j];
      a[i][j] = sum;
      if( (dum=vv[i]*fabs(sum)) >= big){
        big  = dum;
        imax = i;
      }
    }
    if(j != imax){
      for(k=0;k<n;k++){
        dum = a[imax][k];
        a[imax][k] = a[j][k];
        a[j][k] = dum;
      }
      d = -d;
      vv[imax] = vv[j];
    }
    indx[j] = imax;
    if(a[j][j] == 0.) a[j][j] = TINY;
    
    if(j != n-1){
      dum = 1./a[j][j];
      for(i=j+1;i<n;i++) a[i][j] *= dum;
    }
  }
  delete vv;
}

/**
 * LU decomposition is a useful way to invert matrices and solve linear
 * equations. ludcmp is the main decomposition routine. The matrix \c a is modified
 * to the LU decomposition of a row-wise permutation of itself. \c indx
 * records the permutation and \c d is returned as +/- 1 depending upon
 * the odd/even nature of the row interchanges.
 * \param a    the square matrix to be decomposed
 * \param indx index array to record the permutations made
 * \param d    +/- 1 depending on sign of permutation
 * \exception Throws a Subs::Subs_Error if the matrix is singular or if there
 * are problems with the dimensions of the Arrays
 */

void Subs::ludcmp(Buffer2D<double>& a, Buffer1D<size_t>& indx, double &d){

  size_t n = a.nrow();
  if(!n) throw Subs_Error("void Subs::ludcmp(Buffer2D<double>&, Buffer1D<size_t>&, double&): null matrix");
  if(int(n) != a.ncol()) throw Subs_Error("void Subs::ludcmp(Buffer2D<double>&, Buffer1D<size_t>&, double&): matrix not square");
  
  size_t i, imax=0, j, k;
  double big, dum, sum, temp;
  const double TINY = 1.e-20;
  Buffer1D<double> vv(n);
  indx.resize(n);

  d  = 1.;
  for(i=0;i<n;i++){
    big = 0.;
    for(j=0;j<n;j++) 
      if((temp=fabs(a[i][j])) > big) big = temp;
    if(big == 0.)
      throw Subs_Error("void Subs::ludcmp(Buffer2D<double>&, Buffer1D<size_t>&, double&): singular matrix");

    vv[i] = 1./big;
  }
  for(j=0;j<n;j++){
    for(i=0;i<j;i++){
      sum = a[i][j];
      for(k=0;k<i;k++) sum -= a[i][k]*a[k][j];
      a[i][j] = sum;
    }
    big = 0.;
    for(i=j;i<n;i++){
      sum = a[i][j];
      for(k=0;k<j;k++) sum -= a[i][k]*a[k][j];
      a[i][j] = sum;
      if( (dum=vv[i]*fabs(sum)) >= big){
        big  = dum;
        imax = i;
      }
    }
    if(j != imax){
      for(k=0;k<n;k++){
        dum = a[imax][k];
        a[imax][k] = a[j][k];
        a[j][k] = dum;
      }
      d = -d;
      vv[imax] = vv[j];
    }
    indx[j] = imax;
    if(a[j][j] == 0.) a[j][j] = TINY;
    
    if(j != n-1){
      dum = 1./a[j][j];
      for(i=j+1;i<n;i++) a[i][j] *= dum;
    }
  }
}

 /**
  * ludksb performs back substitution given the results from
  * \c ludcmp to solve a set of linear equations.
  * \param a the matrix returned by ludcmp
  * \param n its dimension (n by n)
  * \param indx the index array returned by ldcmp
  * \param b the vector of the right-hand side, returned as the solution.
  */

void Subs::lubksb(double **a, unsigned int n, unsigned int *indx, double *b){
  int i, ii=-1, ip, j;
  double sum;
  
  for(i=0;i<int(n);i++){
    ip    = indx[i];
    sum   = b[ip];
    b[ip] = b[i];
    if(ii >= 0) 
      for(j=ii;j<i;j++) 
	sum -= a[i][j]*b[j];
    else if(sum) 
      ii = i;
    b[i] = sum;
  }
  for(i=int(n-1);i>=0;i--){
    sum = b[i];
    for(j=i+1;j<int(n);j++) sum -= a[i][j]*b[j];
    b[i] = sum/a[i][i];
  }
}

 /**
  * ludksb performs back substitution given the results from
  * \c ludcmp to solve a set of linear equations.
  * \param a the matrix returned by ludcmp
  * \param indx the index array returned by ldcmp
  * \param b the vector of the right-hand side, returned as the solution.
  */

void Subs::lubksb(const Buffer2D<double>& a, const Buffer1D<size_t>& indx, Buffer1D<double>& b){
  int i, ii=-1, ip, j;
  double sum;

  size_t n = a.nrow();
  if(!n) throw Subs_Error("Null matrix in ludksb");
  if(int(n) != a.ncol()) throw Subs_Error("ludksb: matrix not square");
  if(int(n) != indx.size()) throw Subs_Error("a and indx clash in ludksb");

  b.resize(n);
  
  for(i=0;i<int(n);i++){
    ip    = indx[i];
    sum   = b[ip];
    b[ip] = b[i];
    if(ii >= 0) for(j=ii;j<i;j++) sum -= a[i][j]*b[j];
    else if(sum) ii = i;
    b[i] = sum;
  }
  for(i=int(n-1);i>=0;i--){
    sum = b[i];
    for(j=i+1;j<int(n);j++) sum -= a[i][j]*b[j];
    b[i] = sum/a[i][i];
  }
}



