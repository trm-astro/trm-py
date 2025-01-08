#include <string>
#include "trm/subs.h"
#include "trm/buffer2d.h"

/** svdfit uses singular value decomposition to solve linear
 * least squares problems. This is much safer than the "normal" equation
 * method as it can cope with degeneracy. The arguments 
 * u, v and w can be used to obtain the covariance matrix. u will only
 * be ndat by n on output where ndat is the number of valid data points.
 *
 * \param data 1D data array. This includes the x values for convenience
 * of the calling routine rather than necessity. Points with negative
 * errors are ignored. ndata points.
 * \param a the nc fit coefficients
 * \param vect ndata by nc array of function values
 * \param u ndata by nc work array (it will be resized by the routine
 * if there are negative error bars)
 * \param v nc by nc work array (will be resized)
 * \param w nc element vector (will be resized)
 */

double Subs::svdfit(const Buffer1D<rv>& data, Buffer1D<float>& a, const Buffer2D<float>& vect, Buffer2D<float>& u,
		    Buffer2D<float>& v, Buffer1D<float>& w){
    
    if(a.size() != vect.get_nx())
	throw Subs_Error("svdfit[float]: number of coefficients = " + Subs::str(a.size()) + 
			 " in parameter vector does not match number in function array = " + Subs::str(vect.get_nx()));
    if(data.size() != vect.get_ny())
	throw Subs_Error("svdfit[float]: number of data = " + Subs::str(data.size()) + 
			 " does not match number in function array = " + Subs::str(vect.get_ny()));

    size_t ndata = data.size();
    size_t nc    = a.size();
    size_t ndat = 0;
    size_t i, j, k;
    for(j=0; j<ndata; j++)
	if(data[j].z>0.) ndat++;
    
    u.resize(ndat,nc);
    v.resize(nc,nc);
    w.resize(nc);

    const float TOL = 1.e-5;
    float tmp;
    
    Buffer1D<float> b(ndat);
    for(i=0,k=0; i<ndata; i++){
	if(data[i].z>0.){
	    tmp = 1./data[i].z;
	    for(j=0; j<nc; j++)
		u[k][j] = tmp*vect[i][j];
	    b[k] = tmp*data[i].y;
	    k++;
	}
    }
    
    svdcmp(u, w, v);
    
    // Edit singular values
    float wmax = 0.;
    for(i=0; i<nc; i++) 
	if(w[i] > wmax) wmax = w[i];
    
    float thresh = TOL*wmax;
    for(i=0; i<nc; i++)
	if(w[i] < thresh) w[i] = 0.;
    
    // Carry on
    svbksb(u, w, v, b, a);
    
    double sum, chisq = 0.;
    
    for(i=0,k=0; i<ndata; i++){
	if(data[i].z>0.){
	    for(j=0, sum=0.; j<nc; j++)
		sum += a[j]*vect[i][j];
	    chisq += sqr((data[i].y-sum)/data[i].z);
	}
    }
    return chisq;
}

/** svdfit uses singular value decomposition to solve linear
 * least squares problems. This is much safer than the "normal" equation
 * method as it can cope with degeneracy. The arguments 
 * u, v and w can be used to obtain the covariance matrix. u will only
 * be ndat by n on output where ndat is the number of valid data points.
 *
 * \param data 1D data array. This includes the x values for convenience
 * of the calling routine rather than necessity. Points with negative
 * errors are ignored. ndata points.
 * \param a the nc fit coefficients
 * \param vect ndata by nc array of function values
 * \param u ndata by nc work array (it will be resized by the routine
 * if there are negative error bars)
 * \param v nc by nc work array
 * \param w nc element vector
 */

double Subs::svdfit(const Buffer1D<ddat>& data, Buffer1D<double>& a, const Buffer2D<double>& vect, 
		    Buffer2D<double>& u, Buffer2D<double>& v, Buffer1D<double>& w){
    
    if(a.size() != vect.get_nx())
	throw Subs_Error("svdfit[double]: number of coefficients = " + Subs::str(a.size()) + 
			 " in parameter vector does not match number in function array = " + Subs::str(vect.get_nx()));
    if(data.size() != vect.get_ny())
	throw Subs_Error("svdfit[double]: number of data = " + Subs::str(data.size()) + 
			 " does not match number in function array = " + Subs::str(vect.get_ny()));

    size_t ndata = data.size();
    size_t nc    = a.size();
    size_t ndat = 0;
    size_t i, j, k;
    for(j=0; j<ndata; j++)
	if(data[j].z>0.) ndat++;
    
    u.resize(ndat,nc);
    v.resize(nc,nc);
    w.resize(nc);
    
    const double TOL = 1.e-5;
    double tmp;
    
    Buffer1D<double> b(ndat);
    for(i=0,k=0; i<ndata; i++){
	if(data[i].z>0.){
	    tmp = 1./data[i].z;
	    for(j=0; j<nc; j++)
		u[k][j] = tmp*vect[i][j];
	    b[k] = tmp*data[i].y;
	    k++;
	}
    }
    
    svdcmp(u, w, v);
    
    // Edit singular values
    double wmax = 0.;
    for(i=0; i<nc; i++) 
	if(w[i] > wmax) wmax = w[i];
    
    double thresh = TOL*wmax;
    for(i=0; i<nc; i++)
	if(w[i] < thresh) w[i] = 0.;
    
    // Carry on
    svbksb(u, w, v, b, a);
    
    double sum, chisq = 0.;
    
    for(i=0,k=0; i<ndata; i++){
	if(data[i].z>0.){
	    for(j=0, sum=0.; j<nc; j++)
		sum += a[j]*vect[i][j];
	    chisq += sqr((data[i].y-sum)/data[i].z);
	}
    }
    return chisq;
}

/** Fits sinusoid using singular value decomposition.
 * \param data 1D data array. This includes the x values for convenience
 * of the calling routine rather than necessity. Points with negative
 * errors are ignored. ndata points.
 * \param a The nc fit coefficients
 * \param cosine ndata array of cosine values
 * \param sine ndata array of sine values
 * \param u ndata by nc work array (it will be resized by the routine
 * if there are negative error bars)
 * \param v nc by nc work array
 * \param w nc element vector
 */

double Subs::svdfit(const Buffer1D<rv>& data, Buffer1D<float>& a,
		    const Buffer1D<double>& cosine, const Buffer1D<double>& sine, 
		    Buffer2D<float>& u, Buffer2D<float>& v, Buffer1D<float>& w){
    
    size_t ndata = data.size();
    size_t nc    = a.size();
    if(nc != 3) throw Subs_Error("Expected 3 coefficients in sinusoid svdfit");
    size_t ndat = 0;
    size_t i, k;
    for(i=0; i<ndata; i++)
	if(data[i].z>0.) ndat++;
    
    u.resize(ndat,nc);
    v.resize(nc,nc);
    w.resize(nc);
    
    const float TOL = 1.e-5;
    float tmp;
    
    Buffer1D<float> b(ndat);
    for(i=0,k=0; i<ndata; i++){
	if(data[i].z>0.){
	    tmp     = 1./data[i].z;
	    u[k][0] = tmp;
	    u[k][1] = tmp*cosine[i];
	    u[k][2] = tmp*sine[i];
	    b[k]    = tmp*data[i].y;
	    k++;
	}
    }
    
    svdcmp(u, w, v);
    
    // Edit singular values
    float wmax = 0.;
    for(i=0; i<nc; i++) 
	if(w[i] > wmax) wmax = w[i];
    
    float thresh = TOL*wmax;
    for(i=0; i<nc; i++)
	if(w[i] < thresh) w[i] = 0.;
    
    // Carry on
    svbksb(u, w, v, b, a);
    
    double sum, chisq = 0.;
    
    for(i=0; i<ndata; i++){
	if(data[i].z>0.){
	    sum    = a[0] + a[1]*cosine[i] + a[2]*sine[i];
	    chisq += sqr((data[i].y-sum)/data[i].z);
	}
    }
    return chisq;
}







