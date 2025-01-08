
#include <iostream.h>
#include "../base/subs.h"

void meml33(const int nd, const double amat[6][6], double val[6],
	      double vec[6][6]);

int main(){
  double  d[6][6],v[6],vec[6][6];
  d[0][0] = 1.;
  d[0][1] = d[1][0] = 67788.2;
  d[0][2] = d[2][0] = -0.024412;
  d[1][1] = 1.;
  d[1][2] = d[2][1] = 6628.52;  
  d[2][2] = 1.;
  meml33(3,d,v,vec);
  for(int i=0;i<3;i++)
    std::cout << v[i] << std::endl;
}
      

  void meml33(const int nd, const double amat[6][6], double val[6],
	      double vec[6][6]){

    //
    // purpose:
    //   diagonalisation of symmetric matrix
    //                                 k      k         k
    //   solves      amat(i,j) * vec(j)  = val  * vec(i)
    //   for eigenvalues val and eigenvectors vec, for k=1,2,...,nd
    //
    // parameters:
    //   argument type  i/o  dimension   description
    //    nd      i     i      -         number of dimensions
    //    amat    r*8   i     nd,nd      matrix (from 6*6 array)
    //    val     r*8     o    nd        eigenvalues
    //    vec     r*8     o   nd,nd      eigenvectors
    //
    // globals:
    //   variable  common block  i/o  description
    //    -
    //
    // external calls:
    //    -
    //
    // history:
    //   john skilling    8 nov 1985     initial release
    //
    //
    // notes:
    //   (1) eigenvalues val are returned in increasing order
    //   (2) eigenvectors vec are normalised to v.amat.v=1
    //   (3) input matrix amat is preserved
    //   (4) only the upper triangle j>i of amat(i,j) is read
    //
    //   (5) algorithm is to repeatedly square amat until the largest
    //       eigenvalue dominates, then to subtract off that eigenvector,
    //       and repeat until all the eigenvectors have been found.
    //       this nd**4 algorithm would be too slow for large matrices,
    //       but is robust for small.
    //

    double c, z[6][6], p[6][6], q[6][6], y[6], x[6];

    // eps is related to machine accuracy

    const double eps=3.0e-16;

    // copy amat to a positive definite matrix z (by adding c to eigenvals)
    
    double a=0.0, b=0.0, t=0.0;
    for(int i=0; i<nd; i++){
      a += amat[i][i];
      for(int j=0; j<nd; j++){
	c = amat[i][j];
	if(fabs(c)<1.0e-18) c=0.0;
	b += c*c;
      }
      t += 1.0;
    }
    c = max(b-a*a/t,0.0);
    c = a/t-sqrt(c)-eps*sqrt(b);
    std::cerr << "c = " << c << std::endl;
    for(int i=0; i<nd; i++){
      for(int j=0; j<nd; j++){
	z[i][j] = amat[i][j];
	z[j][i] = amat[i][j];
      }
      z[i][i] -= c;
    }

    // lmax = maximum number of inner loop iterates

    t = -log(eps)/(eps*log(2.0));
    t = log(t)/log(2.0);
    int lmax = int(t+2.0);

    int n = nd-1;
    
    std::cerr << "lmax= " << lmax << std::endl;

    // outer loop over n for successively smaller eigenvalues

    for(;;){
      t=0.0;
      for(int i=0; i<nd; i++){
	for(int j=0; j<nd; j++)
	  p[i][j] = z[i][j];
	t += p[i][i];
      }
      std::cerr << "t = " << t << std::endl;

      int l = 0;

      // inner loop over l for squaring p and setting t=trace

      do{
	//	cerr << "start "<< n << " " << t << endl;
	l++;

	t = 1.0/t;
	for(int i=0; i<nd; i++){
	  for(int j=0; j<nd; j++){
	    q[i][j] = t*p[i][j];
	    if(fabs(q[i][j]) < 1.0e-18) q[i][j] =0.0;
	  }
	}

	t=0.0;
	for(int i=0; i<nd; i++){
	  for(int j=0; j<nd; j++){
	    //	    cerr << "A " << p[i][j] << " " << q[i][j] << endl;
	    a=0.0;
	    for(int k=0; k<nd; k++)
	      a += q[i][k]*q[k][j];
	    if(fabs(a)<1.0e-18) a=0.0;
	    p[i][j] = a;
	    p[j][i] = a;
	  }
	  t += p[i][i];
	}
	std::cerr << "end t = " << t << std::endl;
      }while( t < 1.0-eps*10.0 && l <= lmax );

      // end inner loop when p is dyadic
      //
      // k = largest column = estimate of current largest eigenvector

      int k=0;
      a=0.0;
      for(int i=0; i<nd; i++){
	if(p[i][i] > a){
	  a = p[i][i];
	  k = i;
	}
      }

      // p(p(largest column)) = better estimate x of eigenvector

      for(int i=0; i<nd; i++){
	a=0.0;
	for(int j=0; j<nd; j++)
	  a += p[i][j]*p[j][k];
	
	if(fabs(a) < 1.0e-18) a=0.0;
	y[i]=a;
      }

      t=0.0;
      for(int i=0; i<nd; i++){
	a=0.0;
	for(int j=0; j<nd; j++)
	  a += p[i][j]*y[j];
	if(fabs(a)<1.0e-18) a=0.0;
	t += a*a;
	x[i]=a;
      }

      // repeat..  orthogonalise x to previous eigenvectors

      k=nd-1;
      if(k != n){
	do{
	  a=0.0;
	  for(int i=0; i<nd; i++)
            a += x[i]*vec[i][k];
	  for(int i=0; i<nd; i++){
            b=x[i]-a;
            if(fabs(b)<1.0e-18) b=0.0;
            x[i]=b;
	  }
          k--;
	}while(k>n);
      }

      // ..while
      // normalise eigenvector x

      a=0.0;
      for(int i=0; i<nd; i++)
	a += x[i]*x[i];

      a =1.0/sqrt(a);

      // copy eigenvector x into output array vec

      for(int i=0; i<nd; i++){
	x[i] *= a;
	vec[i][n] = x[i];
      }

      // set eigenvalue val directly from the eigenvector

      a=0.0;
      for(int i=0; i<nd; i++){
	for(int j=0; j<nd; j++){
	  b = x[i]*x[j];
	  if(fabs(b)<1.0e-18) b=0.0;
	  a += z[i][j]*b;
	}
      }
      val[n] = a+c;
      //      cerr << "... n = " << n << ", eval = " << val[n] << endl;

      // finish ?   (if full set of eigenvectors has been found)

      if(--n < 0) return;

      // otherwise, remove dyadic  x.xtranspose  from matrix z

      for(int i=0; i<nd; i++){
	for(int j=i; j<nd; j++){
	  b =x[i]*x[j];
	  if(fabs(b)<1.0e-18) b=0.0;
	  z[i][j] -= a*b;
	  if(fabs(z[i][j])<1.0e-18) z[i][j] =0.0;
	  z[j][i] = z[i][j];
	  std::cerr << "z " << z[i][j] << std::endl;
	}
      }
    }
  }
