#include "trm/memsys.h"
#include "doppler.h"



/* op computes the image to data transform that is the heart of Doppler
 * tomography. i.e it projects the image(s) to compute data corresponding to
 * each spectrum.
 *
 * Arguments::
 *
 *  image  : multiple images, stored in a 1D array for compatibility with
 *           mem routines (input)
 *  nxyz   : dimensions of each image, one/image, (input)
 *  vxy    : Vx-Vy pixel size(s), one/image, (input)
 *  vz     : Vz pixel size(s), one/image, (input)
 *  wavel  : central wavelengths for each image, vector/image, (input)
 *  gamma  : systemic velocities for each wavelength, vector/image, (input)
 *  scale  : scaling factors for each wavelength, vector/image, (input)
 *  itype  : types for each image, (input)
 *  tzero  : ephemeris zero point, same units as the times (input)
 *  period : ephemeris period, same units as the times (input)
 *  quad   : quadratic term of ephemeris, same units as the times (input)
 *  vfine  : km/s to use for internal fine array (input)
 *  sfac   : global scaling factor
 *
 *  data   : computed data, again as a 1D array for mem routines (output)
 *  wave   : wavelength array, matching the data array, (input)
 *  nwave  : number of wavelengths, one/dataset [effectively an X dimension],
 *           (input)
 *  nspec  : number of spectra, one/dataset [effectively a Y dimension], (input)
 *  time   : central times for each spectrum (input)
 *  expose : exposure lengths for each spectrum (input)
 *  nsub   : sub-division factors for each spectrum (input)
 *  fwhm   : FWHM resolutions, kms/s, one/dataset, (input)
 *
 * Note the arguments 'image' through to 'sfac' are associated with the
 * Doppler map file, while the rest are associated with the data. 'data' is
 * the only output argument.
 */

void op(const float* image, const std::vector<Nxyz>& nxyz,
    const std::vector<double>& vxy, const std::vector<double>& vz,
    const std::vector<std::vector<double> >& wavel,
    const std::vector<std::vector<float> >& gamma,
    const std::vector<std::vector<float> >& scale,
    const std::vector<Itype>& itype,
    double tzero, double period, double quad, double vfine, double sfac,
    float* data, const double* wave,
    const std::vector<size_t>& nwave, const std::vector<size_t>& nspec,
    const std::vector<std::vector<double> >& time,
    const std::vector<std::vector<float> >& expose,
    const std::vector<std::vector<int> >& nsub,
    const std::vector<double>& fwhm){

// Each image is first projected onto a finely-spaced array which is
// blurred once the projection is finished and then added into the output
// data array.

// We first need to know how many pixels are needed for the blurring
// function. Loop over the data sets to compute the maximum. Also take
// chance to compute number of data pixels.
size_t nblurr = 0, ndpix = 0;
double psigma;
for(size_t nd=0; nd<nspec.size(); nd++){
    // sigma of blurring in terms of fine array pixels
    psigma = fwhm[nd]/vfine/EFAC;

    // total blurr width will be 2*nblurr-1
    nblurr = std::max(nblurr,size_t(round(6.*psigma)+1));
    ndpix += nwave[nd]*nspec[nd];
}

// Now calculate vmax = 2*(maximum projected velocity) of any pixel
// relative to the (possibly gamma shifted) line centre
size_t ny, nx, nz, nimage = nxyz.size();
double vmax = 0.;
for(size_t nim=0; nim<nimage; nim++){
    nx = nxyz[nim].nx;
    ny = nxyz[nim].ny;
    nz = nxyz[nim].nz;
    if(nz > 1){
        vmax = std::max(std::sqrt(std::pow(vxy[nim],2)
                                  *(std::pow(nx,2)+std::pow(ny,2))+
                                  std::pow(vz[nim]*nz,2)), vmax);
    }else{
        vmax = std::max(std::sqrt(std::pow(vxy[nim],2)
                                  *(std::pow(nx,2)+std::pow(ny,2))), vmax);
    }
}

// NFINE is the number of pixels needed for the fine array at its maximum
// which should be enough in all cases even if its too many on occasion.
// We add in extra pixels beyond the maximum velocity in order to allow
// for blurring.
const int nfine = std::max(5,int(std::ceil(vmax/vfine)+2*nblurr));

// adjust vmax to be the maximum velocity of the fine array measured
// from zero.
vmax = nfine*vfine/2.;

// The blurring is carried out with FFTs requiring zero padding.  Thus the
// actual number of pixels to grab is the nearest power of 2 larger than
// nfine+nblurr Call this NFINE.
const size_t NFINE =
    size_t(std::pow(2,int(std::ceil(std::log(nfine+nblurr)/std::log(2.)))));

// grab memory for blurring array and associated FFTs
double *blurr = fftw_alloc_real(NFINE);

// this for the FFT of the blurring array. This could in fact be done in
// place, but the amount of memory required here is small, so it is better
// to have an explicit dedicated array
const size_t NFFT = NFINE/2+1;
complex_fft *bfft = (complex_fft *)fftw_malloc(sizeof(complex_fft)*NFFT);

// this for the FFT of the fine pixel array
complex_fft *fpfft = (complex_fft *)fftw_malloc(sizeof(complex_fft)*NFFT);

// Grab space for fine arrays
double *fine = fftw_alloc_real(NFINE);

// create plans for blurr, fine pixel and final inverse FFTs
// must be here not inside parallelised loop because they
// are not thread-safe
fftw_plan pblurr = fftw_plan_dft_r2c_1d(NFINE, blurr, RECAST(bfft), FFTW_ESTIMATE);
fftw_plan pforw  = fftw_plan_dft_r2c_1d(NFINE, fine, RECAST(fpfft), FFTW_ESTIMATE);
fftw_plan pback  = fftw_plan_dft_c2r_1d(NFINE,  RECAST(fpfft), fine, FFTW_ESTIMATE);

// various local variables
size_t m, n;
double norm, prf, v1, v2;
double w1, w2, wv;
float gm;

// image velocity steps xy and z
double vxyi, vzi = 0.;

// temporary image & data pointers
const float *iptr;
const double *wptr;

// large series of nested loops coming up.
// First zero the data array
memset(data, 0, ndpix*sizeof(float));

// Loop over each data set. doff is a pointer
// to get to the start of the dataset
for(size_t nd=0; nd<nspec.size(); nd++){

    // compute blurring array (specific to the dataset so can't be done
    // earlier). End by FFT-ing it in order to perform the convolution
    // later.
    norm = blurr[0] = 1.;
    psigma = fwhm[nd]/vfine/EFAC;

    for(m=1, n=NFINE-1; m<nblurr; m++, n--){
        prf = std::exp(-std::pow(m/psigma,2)/2.);
        blurr[m] = blurr[n] = prf;
        norm += 2.*prf;
    }

    // zero the centre part
    memset(blurr+nblurr, 0, (NFINE-2*nblurr+1)*sizeof(double));

    // Next line gets overall FFT/iFFT normalisation right
    // cutting out the need to normalise by NFINE later
    norm *= NFINE;

    // normalise
    blurr[0] /= norm;
    for(m=1, n=NFINE-1; m<nblurr; m++, n--){
        blurr[m] /= norm;
        blurr[n] /= norm;
    }

    // take FFT
    fftw_execute(pblurr);

    // Loop over each image. iptr updated at end of loop
    iptr = image;
    for(size_t ni=0; ni<nxyz.size(); ni++){

        // Calculate whether we can skip this particular
        // dataset / image combination because we can save
        // much fuss now if we can
        bool skip = true;
        for(size_t nw=0; nw<wavel[ni].size(); nw++){
            wv = wavel[ni][nw];
            gm = gamma[ni][nw];

            wptr = wave;
            for(size_t ns=0; ns<nspec[nd]; ns++, wptr+=nwave[nd]){
                w1 = wptr[0];
                w2 = wptr[nwave[nd]-1];
                v1 = CKMS*(w1/wv-1.)-gm;
                v2 = CKMS*(w2/wv-1.)-gm;

                // test for overlap. If yes, break the loop because we
                // can't skip
                if(v1 < vmax && v2 > -vmax){
                    skip = false;
                    break;
                }
            }
            if(!skip) break;
        }
        if(skip) continue;

        // extract image dimensions and velocity scales
        nx   = nxyz[ni].nx;
        ny   = nxyz[ni].ny;
        nz   = nxyz[ni].nz;
        vxyi = vxy[ni];
        if(nz > 1) vzi = vz[ni];

        // loop through each spectrum of the data set

        // Next loop contains the main effort, which should
        // be pretty much the same per spectrum, hence we
        // parallelise it
#pragma omp parallel for
        for(int ns=0; ns<int(nspec[nd]); ns++){

            // declare variables here so they are unique to each thread
            int nf, ifp1, ifp2;
            size_t ix, iy, iz, k, m;
            double cosp, sinp, phase, tsub, corr, deriv, sum;
            double pxoff, pyoff, pzoff, pxstep, pystep, pzstep;
            double weight, itfac, wv, sc, v1, v2, fp1, fp2;
            float gm;
            const float *iiptr;

            // unique pointers for each thread
            float *dptr = data + nwave[nd]*ns;
            const double *wptr = wave + nwave[nd]*ns;

            // this for the FFT of the fine pixel array
            complex_fft *fpffti = (complex_fft *)fftw_malloc(sizeof(complex_fft)*NFFT);

            // Grab space for fine arrays
            double *fine = fftw_alloc_real(NFINE);
            double *tfine = fftw_alloc_real(nfine);

            // Zero the fine array
            memset(fine, 0, NFINE*sizeof(double));

            // Loop over sub-spectra to simulate finite exposures
            int ntdiv = nsub[nd][ns];
            for(int nt=0; nt<ntdiv; nt++){

                // Zero the fine sub array
                memset(tfine, 0, nfine*sizeof(double));

                // Compute phase over uniformly spaced set from start to
                // end of exposure. Times are assumed to be mid-exposure
                tsub = time[nd][ns]+expose[nd][ns]*(float(nt)-float(ntdiv-1)/2.)/
                    std::max(ntdiv-1,1);

                // Linear estimate of phase to start
                phase = (tsub-tzero)/period;

                // Two Newton-Raphson corrections for quadratic term
                for(int nc=0; nc<2; nc++){
                    corr  = tzero+period*phase+quad*std::pow(phase,2) - tsub;
                    deriv = period+2.*quad*phase;
                    phase -= corr/deriv;
                }

                cosp = cos(2.*M_PI*phase);
                sinp = sin(2.*M_PI*phase);

                // Image type factor
                switch(itype[ni])
                {
                case PUNIT:
                    itfac = +1.;
                    break;
                case NUNIT:
                    itfac = -1.;
                    break;
                case PSINE:
                    itfac = +sinp;
                    break;
                case NSINE:
                    itfac = -sinp;
                    break;
                case PCOSINE:
                    itfac = +cosp;
                    break;
                case NCOSINE:
                    itfac = -cosp;
                    break;
                case PSINE2:
                    itfac = +sin(4.*M_PI*phase);
                    break;
                case NSINE2:
                    itfac = -sin(4.*M_PI*phase);
                    break;
                case PCOSINE2:
                    itfac = +cos(4.*M_PI*phase);
                    break;
                case NCOSINE2:
                    itfac = -cos(4.*M_PI*phase);
                    break;
                default:
                    std::cerr << "unrecognised image type; programming error in op" << std::endl;
                    itfac = 1.;
                }

                // Inner loops are coming, so time to pre-compute some
                // stuff for speed.

                // If nz > 1 (3D tomog), assume equal spaced placed in Vz
                // space spaced by vz symmetrically around 0 (there is in
                // addition an overall gamma that is added later). We work
                // in the fine array pixel space (denoted by p) to avoid
                // too many divisions by vfine in the inner loops. A pixel
                // with zero projected velocity should fall on the centre
                // of the fine array at (NFINE-1)/2. Overall systemic
                // velocity shifts are added later when the fine array is
                // added into the data array.

                pxstep = -vxyi/vfine*cosp;
                pystep = +vxyi/vfine*sinp;
                if(nz > 1){
                    pzoff  = double(nfine-1)/2. + vzi*double(nz+1)/2./vfine;
                    pzstep = -vzi/vfine;
                }else{
                    pzoff  = double(nfine-1)/2.;
                    pzstep = 0.;
                }

                // Project onto the fine array. These are the innermost
                // loops which need to be as fast as possible. Each loop
                // initialises a pixel index and a pixel offset that are
                // both added to as the loop progresses. All the
                // projection has to do is calculate where to add into the
                // tfine array, check that it lies within range and add the
                // value in if it is. Per cycle of the innermost loop there
                // are 4 additions (2 integer, 2 double), 1 rounding,
                // 1 conversion to an integer and 2 comparisons.

                iiptr = iptr;
                for(iz=0; iz<nz; iz++){
                    pzoff += pzstep;
                    pyoff  = pzoff - pystep*double(ny+1)/2.;
                    for(iy=0; iy<ny; iy++){
                        pyoff += pystep;
                        pxoff  = pyoff - pxstep*double(nx+1)/2.;
                        for(ix=0; ix<nx; ix++, iiptr++){
                            pxoff += pxstep;
                            nf = int(round(pxoff));
                            if(nf >= 0 && nf < nfine) tfine[nf] += *iiptr;
                        }
                    }
                }

                // Now add in with correct weight (trapezoidal) to the
                // fine buffer The vxyi squared factor is to give a
                // similar intensity regardless of the pixel
                // size. i.e. the pixel intensities should be thought of
                // as being per (km/s)**2. The normalisation also ensures
                // relative independence with respect to the value of
                // ntdiv. Also take opportunity to fold in image type
                // factor
                if(ntdiv > 1 && (nt == 0 || nt == ntdiv - 1)){
                    weight = itfac*std::pow(vxyi,2)/(2*(ntdiv-1));
                }else{
                    weight = itfac*std::pow(vxyi,2)/std::max(1,ntdiv-1);
                }
                for(nf=0; nf<nfine; nf++) fine[nf] += weight*tfine[nf];
            }

            // At this point 'fine' contains the projection of the current
            // image for the current spectrum. We now applying the blurring.

            // Take FFT of fine array
            fftw_execute_dft_r2c(pforw, fine, RECAST(fpffti));

            // multiply the FFT by the FFT of the blurring (in effect a
            // convolution)
            for(k=0; k<NFFT; k++) fpffti[k] *= bfft[k];

            // Take the inverse FFT
            fftw_execute_dft_c2r(pback, RECAST(fpffti), fine);

            // We now need to add the blurred array into the spectrum once
            // for each wavelength associated with the current image. Do
            // this by calculating the start and end velocities of each
            // data pixel relative to the projection. The first and last
            // pixels of the data array are ignored because we don't have
            // the surrounding pixels needed to compute their extent.

            // loop over each line associated with this image
            for(k=0; k<wavel[ni].size(); k++){
                wv = wavel[ni][k];
                gm = gamma[ni][k];
                sc = sfac*scale[ni][k];

                // left-hand side of first pixel
                v1 = CKMS*((wptr[0]+wptr[1])/2./wv)-gm;

                // loop over every pixel in the spectrum, except
                // first and last.
                for(m=1; m<nwave[nd]-1; m++){
                    // velocity of right-hand side of pixel
                    v2 = CKMS*((wptr[m]+wptr[m+1])/2./wv-1.)-gm;

                    if(v1 < vmax && v2 > -vmax){

                        // fp1, fp2 -- start and end limits of data pixel
                        // in fine array pixels
                        fp1  = v1/vfine + double(nfine-1)/2.;
                        fp2  = v2/vfine + double(nfine-1)/2.;

                        // ifp1, ifp2 -- fine pixel range fully inside the
                        // data in traditional C form (i.e. ifp1<= <ifp2)
                        ifp1 = std::min(nfine,std::max(0,int(std::ceil(fp1+0.5))));
                        ifp2 = std::min(nfine,std::max(0,int(std::floor(fp2+0.5))));

                        // add full pixels
                        sum = 0.;
                        for(nf=ifp1; nf<ifp2; nf++) sum += fine[nf];

                        // add partial pixels
                        if(ifp1 > 0) sum += (ifp1-0.5-fp1)*fine[ifp1-1];
                        if(ifp2 < nfine) sum += (fp2-ifp2+0.5)*fine[ifp2];

                        // finally add it in with scaling dividing by
                        // Delta v / lambda representing the frequency
                        // width
                        dptr[m] += sc*sum/(std::abs(v2-v1)/wptr[m]);
                    }

                    // move on a click
                    v1 = v2;
                }
            }

            // clear memory.
            fftw_free(fine);
            fftw_free(tfine);
            fftw_free(fpffti);

        } // end of parallel section

        // advance the image pointer
        iptr += nz*ny*nx;
    }

    // advance the data and wavelength pointers
    data += nspec[nd]*nwave[nd];
    wave += nspec[nd]*nwave[nd];
}

// cleanup in reverse order of allocation. Have to take
// care according to how memeory was allocated.
fftw_destroy_plan(pback);
fftw_destroy_plan(pforw);
fftw_destroy_plan(pblurr);
fftw_free(fine);
fftw_free(fpfft);
fftw_free(bfft);
fftw_free(blurr);
}

/* 'tr' is the counterpart of op calculating its transpose.
*
* Arguments::
*
*  image  : multiple images, stored in a 1D array for compatibility with
*           mem routines (output)
*  nxyz   : dimensions of each image, one/image, (input)
*  vxy    : Vx-Vy pixel size(s), one/image, (input)
*  vz     : Vz pixel size(s), one/image, (input)
*  wavel  : central wavelengths for each image, vector/image, (input)
*  gamma  : systemic velocities for each wavelength, vector/image, (input)
*  scale  : scaling factors for each wavelength, vector/image, (input)
*  itype  : types for each image, (input)
*  tzero  : ephemeris zero point, same units as the times (input)
*  period : ephemeris period, same units as the times (input)
*  quad   : quadratic term of ephemeris
*  vfine  : km/s to use for internal fine array (input)
*  sfac   : global scaling factor to keep numbers within nice range
*
*  data   : computed data, again as a 1D array for mem routines (input)
*  wave   : wavelength array, matching the data array, (input)
*  nwave  : number of wavelengths, one/dataset [effectively an X dimension],
*           (input)
*  nspec  : number of spectra, one/dataset [effectively a Y dimension], (input)
*  time   : central times for each spectrum (input)
*  expose : exposure lengths for each spectrum (input)
*  nsub   : sub-division factors for each spectrum (input)
*  fwhm   : FWHM resolutions, kms/s, one/dataset, (input)
*
* Note the arguments 'image' through to 'vfine' are associated with the
* Doppler map file, while the rest are associated with the data. 'data' is
* the only output argument.
*/

void tr(float* image, const std::vector<Nxyz>& nxyz,
    const std::vector<double>& vxy, const std::vector<double>& vz,
    const std::vector<std::vector<double> >& wavel,
    const std::vector<std::vector<float> >& gamma,
    const std::vector<std::vector<float> >& scale,
    const std::vector<Itype>& itype,
    double tzero, double period, double quad, double vfine, double sfac,
    const float* data, const double* wave,
    const std::vector<size_t>& nwave, const std::vector<size_t>& nspec,
    const std::vector<std::vector<double> >& time,
    const std::vector<std::vector<float> >& expose,
    const std::vector<std::vector<int> >& nsub,
    const std::vector<double>& fwhm){

// See op for what is going on. This routine computes a transposed version
// which is hard to explain except by saying "here we carry out the
// transpose of what happens in op". The actual code ends up looking very
// similar, bar a reversal of ordering.

// Now we need to know how many pixels at most are needed for the blurring
// function. Loop over the data sets. Note that cf op, we don't compute
// number of data pixels
size_t nblurr = 0;
double psigma;
for(size_t nd=0; nd<nspec.size(); nd++){
    // sigma of blurring in terms of fine array pixels
    psigma = fwhm[nd]/vfine/EFAC;

    // total blurr width will be 2*nblurr-1
    nblurr = std::max(nblurr,size_t(round(6.*psigma)+1));
}

// cf op: we compute number of *image* pixels
size_t nimage = nxyz.size(), nx, ny, nz, nipix=0;
double vmax = 0.;
for(size_t nim=0; nim<nimage; nim++){
    nx = nxyz[nim].nx;
    ny = nxyz[nim].ny;
    nz = nxyz[nim].nz;
    if(nz > 1){
        vmax = std::max(std::sqrt(std::pow(vxy[nim],2)
                                  *(std::pow(nx,2)+std::pow(ny,2))+
                                  std::pow(vz[nim]*nz,2)), vmax);
    }else{
        vmax = std::max(std::sqrt(std::pow(vxy[nim],2)
                                  *(std::pow(nx,2)+std::pow(ny,2))), vmax);
    }
    nipix += nx*ny*nz;
}

// NFINE is the number of pixels needed for the fine array.
const int nfine = std::max(5,int(std::ceil(vmax/vfine)+2*nblurr));

// adjust vmax to be the maximum velocity of the fine array (from zero not
// full extent)
vmax = nfine*vfine/2.;

// The blurring is carried out with FFTs requiring zero padding.  Thus the
// actual number of pixels to grab is the nearest power of 2 larger than
// nfine+nblurr. Call this NFINE.
const size_t NFINE =
    size_t(std::pow(2,int(std::ceil(std::log(nfine+nblurr)/std::log(2.)))));

// grab memory for blurring array and associated FFTs
double *blurr = fftw_alloc_real(NFINE);

// this for the FFT of the blurring array. This could in fact be done in
// place, but the amount of memory required here is small, so it is better
// to have an explicit dedicated array
const size_t NFFT = NFINE/2+1;
complex_fft *bfft = (complex_fft *)fftw_malloc(sizeof(complex_fft)*NFFT);

// this for the FFT of the fine pixel array
complex_fft *fpfft = (complex_fft *)fftw_malloc(sizeof(complex_fft)*NFFT);

// Grab space for fine arrays
double *fine  = fftw_alloc_real(NFINE);
double *tfine = fftw_alloc_real(nfine);

// create plans for the blurring, fine pixel and final inverse FFTs
// (must be called outside multi-threaded parts)
fftw_plan pblurr = fftw_plan_dft_r2c_1d(NFINE, blurr, RECAST(bfft), FFTW_ESTIMATE);
fftw_plan pforw  = fftw_plan_dft_r2c_1d(NFINE, fine, RECAST(fpfft), FFTW_ESTIMATE);
fftw_plan pback  = fftw_plan_dft_c2r_1d(NFINE, RECAST(fpfft), fine, FFTW_ESTIMATE);

// various local variables
size_t m, n, k;
double norm, prf;
double w1, w2, wv, cosp, sinp, phase, tsub, corr, deriv;
double add,pxstep, pystep, pzstep, weight, itfac, sc;
double v1, v2, fp1, fp2;
float gm;
int ifp1, ifp2;

// image velocity steps xy and z
double vxyi, vzi = 0.;

// temporary image & data pointers
float *iptr;
const float *dptr;
const double *wptr;

// large series of nested loops coming up.
// First zero the image array [cf op]
memset(image, 0, nipix*sizeof(float));

// Loop over each data set
for(size_t nd=0; nd<nspec.size(); nd++){

    // compute blurring array (specific to the dataset so can't be done
    // earlier). End by FFT-ing it in order to perform the convolution
    // later.
    norm = blurr[0] = 1.;
    psigma = fwhm[nd]/vfine/EFAC;
    for(m=1, n=NFINE-1; m<nblurr; m++, n--){
        prf = std::exp(-std::pow(m/psigma,2)/2.);
        blurr[m] = blurr[n] = prf;
        norm += 2.*prf;
    }

    // zero the centre part
    memset(blurr+nblurr, 0, (NFINE-2*nblurr+1)*sizeof(double));

    // Next line gets overall FFT/iFFT normalisation right cutting out any
    // need to normalise by NFINE later
    norm *= NFINE;

    // normalise
    blurr[0] /= norm;
    for(m=1, n=NFINE-1; m<nblurr; m++, n--){
        blurr[m] /= norm;
        blurr[n] /= norm;
    }

    // take FFT
    fftw_execute(pblurr);

    // Loop over each image. iptr updated at end of loop
    iptr = image;
    for(size_t ni=0; ni<nxyz.size(); ni++){

        // Calculate whether we can skip this particular dataset / image
        // combination because we can save much fuss now if we can
        bool skip = true;
        for(size_t nw=0; nw<wavel[ni].size(); nw++){
            wv = wavel[ni][nw];
            gm = gamma[ni][nw];

            wptr = wave;
            for(size_t ns=0; ns<nspec[nd]; ns++, wptr+=nwave[nd]){
                w1 = wptr[0];
                w2 = wptr[nwave[nd]-1];
                v1 = CKMS*(w1/wv-1.)-gm;
                v2 = CKMS*(w2/wv-1.)-gm;

                // test for overlap. If yes, break the loop because we
                // can't skip
                if(v1 < vmax && v2 > -vmax){
                    skip = false;
                    break;
                }
            }
            if(!skip) break;
        }
        if(skip) continue;

        // extract image dimensions and velocity scales
        nx   = nxyz[ni].nx;
        ny   = nxyz[ni].ny;
        nz   = nxyz[ni].nz;
        vxyi = vxy[ni];
        if(nz > 1) vzi = vz[ni];

        // loop through each spectrum of the data set
        for(int ns=0; ns<int(nspec[nd]); ns++){
            dptr = data + nwave[nd]*ns;
            wptr = wave + nwave[nd]*ns;

            // Zero the fine array
            memset(fine, 0, NFINE*sizeof(double));

            // [cf op]. This is the final step in op, here it comes the
            // first

            // loop over each line associated with this image
            for(k=0; k<wavel[ni].size(); k++){
                wv = wavel[ni][k];
                gm = gamma[ni][k];
                sc = sfac*scale[ni][k];

                // left-hand side of first pixel
                v1 = CKMS*((wptr[0]+wptr[1])/2./wv)-gm;

                // loop over every pixel in the spectrum.
                for(m=1; m<nwave[nd]-1; m++){
                    // velocity of right-hand side of pixel
                    v2 = CKMS*((wptr[m]+wptr[m+1])/2./wv-1.)-gm;

                    if(v1 < vmax && v2 > -vmax){

                        // fp1, fp2 -- start and end limits of data pixel
                        // in fine array pixels
                        fp1  = v1/vfine + double(nfine-1)/2.;
                        fp2  = v2/vfine + double(nfine-1)/2.;

                        // ifp1, ifp2 -- fine pixel range fully inside the
                        // data in traditional C form (i.e. ifp1<= <ifp2)
                        ifp1 = std::min(nfine,std::max(0,int(std::ceil(fp1+0.5))));
                        ifp2 = std::min(nfine,std::max(0,int(std::floor(fp2+0.5))));

                        // [cf op]
                        add = sc*dptr[m]/(std::abs(v2-v1)/wptr[m]);
                        for(int nf=ifp1; nf<ifp2; nf++) fine[nf] += add;

                        // add partial pixels
                        if(ifp1 > 0) fine[ifp1-1] += (ifp1-0.5-fp1)*add;
                        if(ifp2 < nfine) fine[ifp2] += (fp2-ifp2+0.5)*add;

                    }

                    // move on a click
                    v1 = v2;
                }
            }

            // next section, blurring of fine array, stays same cf op

            // Take FFT of fine array
            fftw_execute_dft_r2c(pforw, fine, RECAST(fpfft));

            // multiply the FFT by the FFT of the blurring (in effect a
            // convolution)
            for(k=0; k<NFFT; k++) fpfft[k] *= bfft[k];

            // Take the inverse FFT
            fftw_execute_dft_c2r(pback, RECAST(fpfft), fine);

            // Loop over sub-spectra to simulate finite exposures
            int ntdiv = nsub[nd][ns];
            for(int nt=0; nt<ntdiv; nt++){

                // Compute phase over uniformly spaced set from start to
                // end of exposure. Times are assumed to be mid-exposure
                tsub = time[nd][ns]+expose[nd][ns]*(float(nt)-float(ntdiv-1)/2.)/
                    std::max(ntdiv-1,1);

                // Linear estimate of phase to start
                phase = (tsub-tzero)/period;

                // Two Newton-Raphson corrections for quadratic term
                for(int nc=0; nc<2; nc++){
                    corr   = tzero+period*phase+quad*std::pow(phase,2) - tsub;
                    deriv  = period+2.*quad*phase;
                    phase -= corr/deriv;
                }

                cosp  = cos(2.*M_PI*phase);
                sinp  = sin(2.*M_PI*phase);

                // Image type factor
                switch(itype[ni])
                {
                case PUNIT:
                    itfac = +1.;
                    break;
                case NUNIT:
                    itfac = -1.;
                    break;
                case PSINE:
                    itfac = +sinp;
                    break;
                case NSINE:
                    itfac = -sinp;
                    break;
                case PCOSINE:
                    itfac = +cosp;
                    break;
                case NCOSINE:
                    itfac = -cosp;
                    break;
                case PSINE2:
                    itfac = +sin(4.*M_PI*phase);
                    break;
                case NSINE2:
                    itfac = -sin(4.*M_PI*phase);
                    break;
                case PCOSINE2:
                    itfac = +cos(4.*M_PI*phase);
                    break;
                case NCOSINE2:
                    itfac = -cos(4.*M_PI*phase);
                    break;
                default:
                    std::cerr << "unrecognised image type; programming error in tr" << std::endl;
                    itfac = 1.;
                }

                // Inner loops are coming, so time to pre-compute some
                // stuff for speed.

                // If nz > 1 (3D tomog), assume equal spaced placed in Vz
                // space spaced by vz symmetrically around 0 (there is in
                // addition an overall gamma that is added later). We work
                // in the fine array pixel space (denoted by p) to avoid
                // too many divisions by vfine in the inner loops. A pixel
                // with zero projected velocity should fall on the centre
                // of the fine array at (NFINE-1)/2. Overall systemic
                // velocity shifts are added later when the fine array is
                // added into the data array.

                pxstep = -vxyi/vfine*cosp;
                pystep =  vxyi/vfine*sinp;
                if(nz > 1){
                    pzstep = -vzi/vfine;
                }else{
                    pzstep = 0.;
                }

                // Project onto the fine array. These are the innermost
                // loops which need to be as fast as possible. Each loop
                // initialises a pixel index and a pixel offset that are
                // both added to as the loop progresses. All the
                // projection has to do is calculate where to add into the
                // tfine array, check that it lies within range and add
                // the value in if it is. Per cycle of the innermost loop
                // there are 4 additions (2 integer, 2 double), 1
                // rounding, 1 conversion to an integer and 2 comparisons.

                // [cf op]
                if(ntdiv > 1 && (nt == 0 || nt == ntdiv - 1)){
                    weight = itfac*std::pow(vxyi,2)/(2*(ntdiv-1));
                }else{
                    weight = itfac*std::pow(vxyi,2)/std::max(1,ntdiv-1);
                }
                for(int nf=0; nf<nfine; nf++) tfine[nf] = weight*fine[nf];

                if(nz < 4){
                    // Small nz, parallelize the y loop
                    double pzoff  = double(nfine-1)/2. - pzstep*double(nz+1)/2.;
                    for(size_t iz=0; iz<nz; iz++){
                        pzoff += pzstep;
                        double pyoff0  = pzoff - pystep*double(ny-1)/2.;
                        float *iiptr = iptr + ny*nx*iz;

#pragma omp parallel for
                        for(int iy=0; iy<int(ny); iy++){
                            int nf;
                            double pyoff = pyoff0 + pystep*iy;
                            double pxoff = pyoff  - pxstep*double(nx+1)/2.;
                            float *iiiptr = iiptr + nx*iy;
                            for(size_t ix=0; ix<nx; ix++, iiiptr++){
                                pxoff += pxstep;
                                nf = int(round(pxoff));
                                if(nf >= 0 && nf < nfine) *iiiptr += tfine[nf];
                            }
                        } // end of parallel section
                    }

                }else{

                    // Large nz, parallelize the z loop
                    double pzoff0  = double(nfine-1)/2. - \
                        pzstep*double(nz-1)/2.;
#pragma omp parallel for
                    for(int iz=0; iz<int(nz); iz++){
                        int nf;
                        double pxoff;
                        double pzoff = pzoff0 + pzstep*iz;
                        double pyoff = pzoff - pystep*double(ny+1)/2.;
                        float *iiptr = iptr + ny*nx*iz;
                        for(size_t iy=0; iy<ny; iy++){
                            pyoff += pystep;
                            pxoff  = pyoff  - pxstep*double(nx+1)/2.;
                            float *iiiptr = iiptr + nx*iy;
                            for(size_t ix=0; ix<nx; ix++, iiiptr++){
                                pxoff += pxstep;
                                nf = int(round(pxoff));
                                if(nf >= 0 && nf < nfine) *iiiptr += tfine[nf];
                            }
                        }
                    } // end of parallel section
                }
            }
        }

        // advance the image pointer
        iptr += nz*ny*nx;
    }

    // advance the data and wavelength pointers
    data += nspec[nd]*nwave[nd];
    wave += nspec[nd]*nwave[nd];
}

// cleanup in reverse order of allocation
fftw_destroy_plan(pback);
fftw_destroy_plan(pforw);
fftw_destroy_plan(pblurr);
fftw_free(tfine);
fftw_free(fine);
fftw_free(fpfft);
fftw_free(bfft);
fftw_free(blurr);
}

/* gaussdef computes a gaussian default image by blurring in all three
 * directions one after the other. The blurring is carried out using FFTs.
 * This routine applies the blurring to one image
 *
 * input  : input image (input)
 * nxyz   : dimensions. Only axes with dimensions > 1 are blurred. (input)
 * fwhmx  : FWHM blurr in X (pixels), <=0 to ignore. (input)
 * fwhmy  : FWHM blurr in Y (pixels), <=0 to ignore. (input)
 * fwhmz  : FWHM blurr in Z (pixels), <=0 to ignore. (input)
 * output : output (blurred) image (output)
 */

void gaussdef(const float *input, const Nxyz& nxyz, double fwhmx,
    double fwhmy, double fwhmz, float *output){

    // copy input to output
    memcpy(output, input, nxyz.ntot()*sizeof(float));

    // some repeatedly used variables
    double norm, sigma, prf;
    size_t ix, iy, iz, nstep, k, m, n;
    size_t nadd, ntot, NTOT, NFFT;
    float *iptr, *iiptr, *optr, *ooptr;
    double *array, *blurr;
    complex_fft *bfft, *afft;
    fftw_plan pblurr, pforw, pback;

    // Blurr in X
    if(nxyz.nx > 1 && fwhmx >= 0.){

    // Work out buffer size needed for FFTs
    nadd = 2*size_t(3.*fwhmx+1.);
    ntot = nxyz.nx + nadd;
    NTOT = size_t(std::pow(2,int(std::ceil(std::log(ntot)/std::log(2.)))));
    NFFT = NTOT/2 + 1;

    // input workspace
    array = fftw_alloc_real(NTOT);

    // blurring array
    blurr = fftw_alloc_real(NTOT);

    // FFT of blurring array
    bfft = (complex_fft *)fftw_malloc(sizeof(complex_fft)*NFFT);

    // FFT of the array to be blurred
    afft = (complex_fft *)fftw_malloc(sizeof(complex_fft)*NFFT);

    // create FFT plans for the blurr, work and final inverse
    pblurr = fftw_plan_dft_r2c_1d(NTOT, blurr, RECAST(bfft), FFTW_ESTIMATE);
    pforw = fftw_plan_dft_r2c_1d(NTOT, array, RECAST(afft), FFTW_ESTIMATE);
    pback = fftw_plan_dft_c2r_1d(NTOT, RECAST(afft), array, FFTW_ESTIMATE);

    // all memory allocated. Get on with it.
    norm = blurr[0] = 1.;
    sigma = fwhmx/EFAC;
    for(m=1, n=NTOT-1; m<nadd/2; m++, n--){
    prf = std::exp(-std::pow(m/sigma,2)/2.);
    blurr[m] = blurr[n] = prf;
    norm += 2.*prf;
    }

    // zero the middle part
    memset(blurr+nadd/2, 0, (NTOT-nadd+1)*sizeof(double));

    // next line gets the FFT/inverse-FFT pair normalisation right.
    norm *= NTOT;

    // normalise
    for(m=0;m<NTOT;m++)
    blurr[m] /= norm;

    // FFT the blurring array
    fftw_execute(pblurr);

    iptr = output;
    optr = output;
    nstep = nxyz.nx;

    for(iz=0; iz<nxyz.nz; iz++){
    for(iy=0; iy<nxyz.ny; iy++, iptr+=nstep, optr+=nstep){

        // transfer data to double work array
        for(ix=0; ix<nxyz.nx; ix++)
            array[ix] = double(iptr[ix]);

        // zeropad
        memset(array+nxyz.nx, 0, (NTOT-nxyz.nx)*sizeof(double));

        // FFT
        fftw_execute(pforw);

        // multiply by the FFT of the blurr
        for(k=0; k<NFFT; k++) afft[k] *= bfft[k];

        // inverse FFT
        fftw_execute(pback);

        // transfer result to output
        for(ix=0; ix<nxyz.nx; ix++)
            optr[ix] = float(array[ix]);
    }
    }

    // Recover memory
    fftw_destroy_plan(pback);
    fftw_destroy_plan(pforw);
    fftw_destroy_plan(pblurr);
    fftw_free(afft);
    fftw_free(bfft);
    fftw_free(blurr);
    fftw_free(array);
    }

    // Blurr in Y
    if(nxyz.ny > 1 && fwhmy >= 0.){

    // Work out buffer size needed for FFTs
    nadd = 2*int(3.*fwhmy+1.);
    ntot = nxyz.ny + nadd;
    NTOT = size_t(std::pow(2,int(std::ceil(std::log(ntot)/std::log(2.)))));
    NFFT = NTOT/2 + 1;

    // input workspace
    array = fftw_alloc_real(NTOT);

    // blurring array
    blurr = fftw_alloc_real(NTOT);

    // FFT of blurring array
    bfft = (complex_fft *)fftw_malloc(sizeof(complex_fft)*NFFT);

    // FFT of the array to be blurred
    afft = (complex_fft *)fftw_malloc(sizeof(complex_fft)*NFFT);

    // create FFT plans for the blurr, work and final inverse
    pblurr = fftw_plan_dft_r2c_1d(NTOT, blurr, RECAST(bfft), FFTW_ESTIMATE);
    pforw  = fftw_plan_dft_r2c_1d(NTOT, array, RECAST(afft), FFTW_ESTIMATE);
    pback = fftw_plan_dft_c2r_1d(NTOT, RECAST(afft), array, FFTW_ESTIMATE);

    // all memory allocated, get on with it
    norm  = blurr[0] = 1.;
    sigma = fwhmy/EFAC;
    for(m=1, n=NTOT-1; m<nadd/2; m++, n--){
    prf = std::exp(-std::pow(m/sigma,2)/2.);
    blurr[m] = blurr[n] = prf;
    norm += 2.*prf;
    }

    // zero the middle part
    memset(blurr+nadd/2, 0, (NTOT-nadd+1)*sizeof(double));

    // next line gets the FFT/inverse-FFT pair normalisation right.
    norm *= NTOT;

    // normalise
    for(m=0;m<NTOT;m++)
    blurr[m] /= norm;

    // FFT the blurring array
    fftw_execute(pblurr);
    nstep = nxyz.nx;
    for(iz=0; iz<nxyz.nz; iz++){
    iiptr = output + nxyz.nx*nxyz.ny*iz;
    ooptr = output + nxyz.nx*nxyz.ny*iz;

    for(ix=0; ix<nxyz.nx; ix++, iiptr++, ooptr++){

        // transfer data to double work array
        iptr = iiptr;
        for(iy=0; iy<nxyz.ny; iy++, iptr+=nstep)
            array[iy] = double(*iptr);

        // zeropad
        memset(array+nxyz.ny, 0, (NTOT-nxyz.ny)*sizeof(double));

        // FFT
        fftw_execute(pforw);

        // multiply by the FFT of the blurr
        for(k=0; k<NFFT; k++) afft[k] *= bfft[k];

        // inverse FFT
        fftw_execute(pback);

        // transfer result to output
        optr = ooptr;
        for(iy=0; iy<nxyz.ny; iy++, optr+=nstep)
            *optr = float(array[iy]);
    }
    }

    // Recover memory
    fftw_destroy_plan(pback);
    fftw_destroy_plan(pforw);
    fftw_destroy_plan(pblurr);
    fftw_free(afft);
    fftw_free(bfft);
    fftw_free(blurr);
    fftw_free(array);
}

// Blurr in Z
if(nxyz.nz > 1 && fwhmz >= 0.){

// Work out buffer size needed for FFTs
nadd = 2*int(3.*fwhmz+1.);
ntot = nxyz.nz + nadd;
NTOT = size_t(std::pow(2,int(std::ceil(std::log(ntot)/std::log(2.)))));
NFFT = NTOT/2 + 1;

// input workspace
array = fftw_alloc_real(NTOT);

// blurring array
blurr = fftw_alloc_real(NTOT);

// FFT of blurring array
bfft = (complex_fft *)fftw_malloc(sizeof(complex_fft)*NFFT);

// FFT of the array to be blurred
afft = (complex_fft *)fftw_malloc(sizeof(complex_fft)*NFFT);

// create FFT plans for the blurr, work and final inverse
pblurr = fftw_plan_dft_r2c_1d(NTOT, blurr, RECAST(bfft), FFTW_ESTIMATE);
pforw  = fftw_plan_dft_r2c_1d(NTOT, array, RECAST(afft), FFTW_ESTIMATE);
pback = fftw_plan_dft_c2r_1d(NTOT, RECAST(afft), array, FFTW_ESTIMATE);

// all memory allocated, get on with it
norm  = blurr[0] = 1.;
sigma = fwhmz/EFAC;
for(m=1, n=NTOT-1; m<nadd/2; m++, n--){
prf = std::exp(-std::pow(m/sigma,2)/2.);
blurr[m] = blurr[n] = prf;
norm += 2.*prf;
}

// zero the middle part
memset(blurr+nadd/2, 0, (NTOT-nadd+1)*sizeof(double));

// next line gets the FFT/inverse-FFT pair normalisation right.
norm *= NTOT;

// normalise
for(m=0;m<NTOT;m++)
blurr[m] /= norm;

// FFT the blurring array
fftw_execute(pblurr);

nstep = nxyz.nx*nxyz.ny;
iiptr = ooptr = output;
for(iy=0; iy<nxyz.ny; iy++){
for(ix=0; ix<nxyz.nx; ix++, iiptr++, ooptr++){

    // transfer data to double work array
    iptr = iiptr;
    for(iz=0; iz<nxyz.nz; iz++, iptr+=nstep)
        array[iz] = double(*iptr);

    // zeropad
    memset(array+nxyz.nz, 0, (NTOT-nxyz.nz)*sizeof(double));

    // FFT
    fftw_execute(pforw);

    // multiply by the FFT of the blurr
    for(k=0; k<NFFT; k++) afft[k] *= bfft[k];

    // inverse FFT
    fftw_execute(pback);

    // transfer result to output
    optr = ooptr;
    for(iz=0; iz<nxyz.nz; iz++, optr+=nstep)
        *optr = float(array[iz]);
}
}

// Get memory back
fftw_destroy_plan(pback);
fftw_destroy_plan(pforw);
fftw_destroy_plan(pblurr);
fftw_free(afft);
fftw_free(bfft);
fftw_free(blurr);
fftw_free(array);
}
}

// Wrap globals to get through to opus
// and tropus inside a namespace

namespace Dopp {
    // For Map objects
    std::vector<Nxyz> nxyz;
    std::vector<DefOpt> def;
    std::vector<Itype> itype;
    std::vector<double> vxy, vz, bias, fwhmxy, fwhmz, squeeze, sqfwhm;
    std::vector<std::vector<double> > wavel;
    std::vector<std::vector<float> > gamma, scale;
    double tzero, period, quad, vfine, sfac;

    // For Data objects
    std::vector<size_t> nwave, nspec;
    std::vector<std::vector<double> > time;
    std::vector<std::vector<float> > expose;
    std::vector<std::vector<int> > nsub;
    std::vector<double> fwhm;
    double *wave;
}

// Implementation of the memsys routines opus and tropus
/* The opus routine needed for memsys
 */
namespace Mem {
    void opus(const int j, const int k){

        std::cerr << "    OPUS " << j+1 << " ---> " << k+1 << std::endl;

        op(Mem::Gbl::st+Mem::Gbl::kb[j], Dopp::nxyz, Dopp::vxy, Dopp::vz,
        Dopp::wavel, Dopp::gamma, Dopp::scale, Dopp::itype,
        Dopp::tzero, Dopp::period, Dopp::quad, Dopp::vfine, Dopp::sfac,
        Mem::Gbl::st+Mem::Gbl::kb[k], Dopp::wave, Dopp::nwave, Dopp::nspec,
        Dopp::time, Dopp::expose, Dopp::nsub, Dopp::fwhm);
    }

    /* The tropus routine needed for memsys
    */
    void tropus(const int k, const int j){

        std::cerr << "  TROPUS " << j+1 << " <--- " << k+1 << std::endl;

        tr(Mem::Gbl::st+Mem::Gbl::kb[j], Dopp::nxyz, Dopp::vxy, Dopp::vz,
        Dopp::wavel, Dopp::gamma, Dopp::scale, Dopp::itype,
        Dopp::tzero, Dopp::period, Dopp::quad, Dopp::vfine, Dopp::sfac,
        Mem::Gbl::st+Mem::Gbl::kb[k], Dopp::wave, Dopp::nwave, Dopp::nspec,
        Dopp::time, Dopp::expose, Dopp::nsub, Dopp::fwhm);
    }
}
// memsys opus and tropus routines implemented


/* npix_data -- returns with the number of pixels needed to allocate memory
 * for the images in a Doppler Map.
 *
 * Arguments::
 *
 *  Map     :  the Doppler image (input)
 *  npix    :  the total number of pixels (output)
 *
 * Returns true/false according to whether it has succeeded. If false
 * it sets an exception so you don't need to set another.
 */

bool npix_map(const py::object& Map, size_t& npix) {
    // Check if Map is None
    if (Map.is_none()) {
        throw std::invalid_argument("doppler.npix_map: Map is NULL");
    }

    npix = 0;

    try {
        // Get the "data" attribute from the Map object
        py::object data = Map.attr("data");

        // Ensure "data" is a sequence
        if (!py::isinstance<py::sequence>(data)) {
            throw std::invalid_argument("doppler.npix_map: Map.data is not a sequence");
        }

        // Iterate over the sequence of Images in "data"
        for (const py::handle& image : data) {
            // Get the "data" attribute of the Image
            py::object idata = image.attr("data");

            // Ensure "idata" is a numpy array
            if (!py::isinstance<py::array>(idata)) {
                throw std::invalid_argument("doppler.npix_map: Image.data is not a numpy array");
            }

            // Check the number of dimensions (2D or 3D)
            py::array idata_array = py::cast<py::array>(idata);
            if (idata_array.ndim() != 2 && idata_array.ndim() != 3) {
                throw std::invalid_argument("doppler.npix_map: Image.data has invalid dimensions (must be 2D or 3D)");
            }

            // Add the total number of elements in the array to npix
            npix += idata_array.size();
        }

        return true; // Success
    } catch (const std::exception& e) {
        throw std::runtime_error(std::string("doppler.npix_map: ") + e.what());
    }
}


bool read_map(py::object& map,
              float* images,
              std::vector<Nxyz>& nxyz,
              std::vector<double>& vxy,
              std::vector<double>& vz,
              std::vector<std::vector<double>>& wave,
              std::vector<std::vector<float>>& gamma,
              std::vector<std::vector<float>>& scale,
              std::vector<Itype>& itype,
              std::vector<DefOpt>& def,
              std::vector<double>& bias,
              std::vector<double>& fwhmxy,
              std::vector<double>& fwhmz,
              std::vector<double>& squeeze,
              std::vector<double>& sqfwhm,
              double& tzero,
              double& period,
              double& quad,
              double& vfine,
              double& sfac) {
    try {
        // Check if the map is None
        if (map.is_none()) {
            throw std::invalid_argument("doppler.read_map: map is NULL");
        }

        // Clear output vectors
        nxyz.clear();
        vxy.clear();
        vz.clear();
        wave.clear();
        gamma.clear();
        scale.clear();
        itype.clear();
        def.clear();
        bias.clear();
        fwhmxy.clear();
        fwhmz.clear();
        squeeze.clear();
        sqfwhm.clear();

        // Extract scalar attributes
        tzero = map.attr("tzero").cast<double>();
        period = map.attr("period").cast<double>();
        quad = map.attr("quad").cast<double>();
        vfine = map.attr("vfine").cast<double>();
        sfac = map.attr("sfac").cast<double>();

        // Extract the "data" attribute (list of images)
        py::list data = map.attr("data");

        // Iterate over each image in the "data" list
        for (const py::handle& image : data) {
            // Extract attributes of the image
            py::array_t<float> idata = image.attr("data").cast<py::array_t<float>>();
            py::array_t<double> iwave = image.attr("wave").cast<py::array_t<double>>();
            py::array_t<float> igamma = image.attr("gamma").cast<py::array_t<float>>();
            py::array_t<float> iscale = image.attr("scale").cast<py::array_t<float>>();
            double ivxy = image.attr("vxy").cast<double>();
            double ivz = image.attr("vz").cast<double>();
            int iitype = image.attr("itype").cast<int>();
            py::object idef = image.attr("default");

            // Extract the "default" attributes
            int doption = idef.attr("option").cast<int>();
            double dbias = idef.attr("bias").cast<double>();
            double dfwhmxy = idef.attr("fwhmxy").cast<double>();
            double dfwhmz = idef.attr("fwhmz").cast<double>();
            double dsqueeze = idef.attr("squeeze").cast<double>();
            double dsqfwhm = idef.attr("sqfwhm").cast<double>();

            // Validate the dimensions of idata (2D or 3D)
            if (idata.ndim() != 2 && idata.ndim() != 3) {
                throw std::invalid_argument("doppler.read_map: Image.data must be 2D or 3D");
            }

            // Add the image data to the output buffer
            std::memcpy(images, idata.data(), idata.nbytes());
            images += idata.size();

            // Store the image dimensions
            if (idata.ndim() == 2) {
                nxyz.emplace_back(idata.shape(1), idata.shape(0));
            } else {
                nxyz.emplace_back(idata.shape(2), idata.shape(1), idata.shape(0));
            }

            // Store wavelengths
            wave.emplace_back(iwave.data(), iwave.data() + iwave.size());

            // Store systemic velocities
            if (iwave.size() != igamma.size()) {
                throw std::invalid_argument("doppler.read_map: wave and gamma sizes do not match");
            }
            gamma.emplace_back(igamma.data(), igamma.data() + igamma.size());

            // Store scale factors
            if (iwave.size() != iscale.size()) {
                throw std::invalid_argument("doppler.read_map: wave and scale sizes do not match");
            }
            scale.emplace_back(iscale.data(), iscale.data() + iscale.size());

            // Store pixel sizes and spacings
            vxy.push_back(ivxy);
            if (idata.ndim() == 3) {
                vz.push_back(ivz);
            }

            // Store image type
            switch (iitype) {
                case 1: itype.push_back(PUNIT); break;
                case 2: itype.push_back(NUNIT); break;
                case 3: itype.push_back(PSINE); break;
                case 4: itype.push_back(NSINE); break;
                case 5: itype.push_back(PCOSINE); break;
                case 6: itype.push_back(NCOSINE); break;
                case 7: itype.push_back(PSINE2); break;
                case 8: itype.push_back(NSINE2); break;
                case 9: itype.push_back(PCOSINE2); break;
                case 10: itype.push_back(NCOSINE2); break;
                default:
                    throw std::invalid_argument("doppler.read_map: Invalid itype value");
            }

            // Store default options
            switch (doption) {
                case 1: def.push_back(UNIFORM); break;
                case 2: def.push_back(GAUSS2D); break;
                case 3: def.push_back(GAUSS3D); break;
                default:
                    throw std::invalid_argument("doppler.read_map: Invalid default option");
            }

            // Store default parameters
            bias.push_back(dbias);
            fwhmxy.push_back(dfwhmxy);
            if (doption == GAUSS3D) {
                fwhmz.push_back(dfwhmz);
                squeeze.push_back(dsqueeze);
                sqfwhm.push_back(dsqfwhm);
            }
        }

        return true; // Success
    } catch (const std::exception& e) {
        throw std::runtime_error(std::string("doppler.read_map: ") + e.what());
    }
}

/* update_map -- copies image data into a Map object for output. The idea is
 * you have created an image inside a C-array that was originally defined by
 * the structure of a Map object using read_map. You now copy it over into the
 * Map object. The very limited nature of the modification here is because
 * most manipulations should be carried out in Python.
 *
 * Arguments::
 *
 *  images : pointer to the images to load into the Map object. Contiguous array
 *         of floats. It must match the total size of the image arrays in map
 *         (input)
 *
 *  map    : the Doppler map (output)
 *
 *  Returns True/False according to success. If False, the outputs above
 *  will not be correctly assigned. If False, a Python exception is raised
 *  and you should return NULL from the calling routine.
 */

bool
update_map(float const* images, py::object& Map){
    if (Map.is_none()) {
        throw std::invalid_argument("doppler.update_map: map is NULL");
    }
    try{
        py::object data = Map.attr("data");
        // Ensure "data" is a sequence
        if (!py::isinstance<py::sequence>(data)) {
            throw std::invalid_argument("doppler.npix_map: Map.data is not a sequence");
        }
        for (const py::handle& image : data){
            py::object idata = image.attr("data");
            
            // Ensure "idata" is a NumPy array
            py::array_t<float> array = py::cast<py::array_t<float>>(idata);

            // Check the number of dimensions (2D or 3D)
            if (array.ndim() != 2 && array.ndim() != 3) {
                throw std::invalid_argument("doppler.update_map: Image.data has invalid dimensions (must be 2D or 3D)");
            }

            // Copy the data from the C++ array into the NumPy array
            std::memcpy(array.mutable_data(), images, array.nbytes());

            // Advance the pointer in the C++ array
            images += array.size();
        }
        return true; // Success
    } catch (const std::exception& e) {
        throw std::runtime_error(std::string("doppler.update_map: ") + e.what());
    }
}


bool npix_data(py::object& Data, size_t& npix) {
    // Check if Data is None
    if (Data.is_none()) {
        throw std::invalid_argument("doppler.npix_data: Data is NULL");
    }

    npix = 0;

    try {
        // Get the "data" attribute from the Data object
        py::object data = Data.attr("data");

        // Ensure "data" is a sequence
        if (!py::isinstance<py::sequence>(data)) {
            throw std::invalid_argument("doppler.npix_data: Data.data is not a sequence");
        }

        // Iterate over the sequence of Spectra objects in "data"
        for (const py::handle& spectra : data) {
            // Get the "flux" attribute of the Spectra
            py::object sflux = spectra.attr("flux");

            // Ensure "sflux" is a NumPy array
            py::array_t<float> flux_array = py::cast<py::array_t<float>>(sflux);

            // Check the number of dimensions (must be 2D)
            if (flux_array.ndim() != 2) {
                throw std::invalid_argument("doppler.npix_data: Spectra.flux must be a 2D array");
            }

            // Add the total number of elements in the array to npix
            npix += flux_array.size();
        }

        return true; // Success
    } catch (const std::exception& e) {
        throw std::runtime_error(std::string("doppler.npix_data: ") + e.what());
    }
}

bool read_data(py::object& Data, float* flux, float* ferr, double* wave,
               std::vector<size_t>& nwave, std::vector<size_t>& nspec,
               std::vector<std::vector<double>>& time,
               std::vector<std::vector<float>>& expose,
               std::vector<std::vector<int>>& nsub,
               std::vector<double>& fwhm) {
    // Check if Data is None
    if (Data.is_none()) {
        throw std::invalid_argument("doppler.read_data: Data is NULL");
    }

    // Temporary pointers for flux, ferr, and wave
    float* tflux = flux;
    float* tferr = ferr;
    double* twave = wave;

    // Clear output vectors
    nwave.clear();
    nspec.clear();
    time.clear();
    expose.clear();
    nsub.clear();
    fwhm.clear();

    try {
        // Get the "data" attribute from the Data object
        py::object data = Data.attr("data");

        // Ensure "data" is a sequence
        if (!py::isinstance<py::sequence>(data)) {
            throw std::invalid_argument("doppler.read_data: Data.data is not a sequence");
        }

        // Iterate over the sequence of Spectra objects in "data"
        for (const py::handle& spectra : data) {
            // Get attributes of the Spectra
            py::array_t<float> sflux = py::cast<py::array_t<float>>(spectra.attr("flux"));
            py::array_t<float> sferr = py::cast<py::array_t<float>>(spectra.attr("ferr"));
            py::array_t<double> swave = py::cast<py::array_t<double>>(spectra.attr("wave"));
            py::array_t<double> stime = py::cast<py::array_t<double>>(spectra.attr("time"));
            py::array_t<float> sexpose = py::cast<py::array_t<float>>(spectra.attr("expose"));
            py::array_t<int> snsub = py::cast<py::array_t<int>>(spectra.attr("nsub"));
            double sfwhm = py::cast<double>(spectra.attr("fwhm"));

            // Check dimensions of flux, ferr, and wave
            if (sflux.ndim() != 2 || sferr.ndim() != 2 || swave.ndim() != 2) {
                throw std::invalid_argument("doppler.read_data: flux, ferr, or wave must be 2D arrays");
            }

            // Ensure dimensions match
            if (sflux.shape(0) != sferr.shape(0) || sflux.shape(1) != sferr.shape(1) ||
                sflux.shape(0) != swave.shape(0) || sflux.shape(1) != swave.shape(1)) {
                throw std::invalid_argument("doppler.read_data: flux, ferr, and wave dimensions do not match");
            }

            // Copy flux data
            std::memcpy(tflux, sflux.data(), sflux.nbytes());
            tflux += sflux.size();

            // Copy flux error data
            std::memcpy(tferr, sferr.data(), sferr.nbytes());
            tferr += sferr.size();

            // Copy wavelength data
            std::memcpy(twave, swave.data(), swave.nbytes());
            twave += swave.size();

            // Store number of spectra and wavelengths
            nspec.push_back(sflux.shape(0));
            nwave.push_back(sflux.shape(1));

            // Copy time data
            std::vector<double> times(stime.data(), stime.data() + stime.size());
            time.push_back(times);

            // Copy exposure data
            std::vector<float> exposes(sexpose.data(), sexpose.data() + sexpose.size());
            expose.push_back(exposes);

            // Copy nsub data
            std::vector<int> nsubs(snsub.data(), snsub.data() + snsub.size());
            nsub.push_back(nsubs);

            // Store FWHM
            fwhm.push_back(sfwhm);
        }

        return true; // Success
    } catch (const std::exception& e) {
        throw std::runtime_error(std::string("doppler.read_data: ") + e.what());
    }
}


bool update_data(const float* flux, py::object& Data) {
    // Check if Data is None
    if (Data.is_none()) {
        throw std::invalid_argument("doppler.update_data: Data is NULL");
    }

    try {
        // Get the "data" attribute from the Data object
        py::object data = Data.attr("data");

        // Ensure "data" is a sequence
        if (!py::isinstance<py::sequence>(data)) {
            throw std::invalid_argument("doppler.update_data: Data.data is not a sequence");
        }

        // Iterate over the sequence of Spectra objects in "data"
        for (const py::handle& spectra : data) {
            // Get the "flux" attribute of the Spectra
            py::object sflux = spectra.attr("flux");

            // Ensure "sflux" is a NumPy array
            py::array_t<float> flux_array = py::cast<py::array_t<float>>(sflux);

            // Check the number of dimensions (must be 2D)
            if (flux_array.ndim() != 2) {
                throw std::invalid_argument("doppler.update_data: Spectra.flux must be a 2D array");
            }

            // Copy the data from the C++ array into the NumPy array
            std::memcpy(flux_array.mutable_data(), flux, flux_array.nbytes());

            // Advance the pointer in the C++ array
            flux += flux_array.size();
        }

        return true; // Success
    } catch (const std::exception& e) {
        throw std::runtime_error(std::string("doppler.update_data: ") + e.what());
    }
}

py::none doppler_comdat(py::object& map, py::object& data) {
    // Check if map or data is None
    if (map.is_none()) {
        throw std::invalid_argument("doppler.comdat: map is NULL");
    }
    if (data.is_none()) {
        throw std::invalid_argument("doppler.comdat: data is NULL");
    }

    // Read sizes
    size_t nipix, ndpix;
    if (!npix_map(map, nipix)) {
        throw std::runtime_error("doppler.comdat: Failed to calculate number of pixels in map");
    }
    if (!npix_data(data, ndpix)) {
        throw std::runtime_error("doppler.comdat: Failed to calculate number of pixels in data");
    }

    // Allocate memory for Map data
    std::vector<float> image(nipix);
    std::vector<Nxyz> nxyz;
    std::vector<DefOpt> def;
    std::vector<double> vxy, vz, bias, fwhmxy, fwhmz, squeeze, sqfwhm;
    std::vector<std::vector<double>> wavel;
    std::vector<std::vector<float>> gamma, scale;
    std::vector<Itype> itype;
    double tzero, period, quad, vfine, sfac;

    // Read the Map
    if (!read_map(map, image.data(), nxyz, vxy, vz, wavel, gamma, scale, itype,
                  def, bias, fwhmxy, fwhmz, squeeze, sqfwhm,
                  tzero, period, quad, vfine, sfac)) {
        throw std::runtime_error("doppler.comdat: Failed to read map");
    }

    // Allocate memory for Data data
    std::vector<size_t> nwave, nspec;
    std::vector<std::vector<double>> time;
    std::vector<std::vector<float>> expose;
    std::vector<std::vector<int>> nsub;
    std::vector<double> fwhm;
    std::vector<float> flux(ndpix);
    std::vector<float> ferr(ndpix);
    std::vector<double> wave(ndpix);

    // Read the Data
    if (!read_data(data, flux.data(), ferr.data(), wave.data(), nwave,
                   nspec, time, expose, nsub, fwhm)) {
        throw std::runtime_error("doppler.comdat: Failed to read data");
    }

    // Release the GIL for computation
    py::gil_scoped_release release;

    // Perform the image-to-data transformation
    op(image.data(), nxyz, vxy, vz, wavel, gamma, scale, itype, tzero, period, quad,
       vfine, sfac, flux.data(), wave.data(), nwave, nspec, time, expose, nsub, fwhm);

    // Restore the GIL
    py::gil_scoped_acquire acquire;

    // Write the modified flux data back into the Data object
    if (!update_data(flux.data(), data)) {
        throw std::runtime_error("doppler.comdat: Failed to update data");
    }

    return py::none();
}


py::none doppler_comdef(py::object& map) {
    // Check if map is None
    if (map.is_none()) {
        throw std::invalid_argument("doppler.comdef: map is NULL");
    }

    // Read the number of pixels in the map
    size_t nipix;
    if (!npix_map(map, nipix)) {
        throw std::runtime_error("doppler.comdef: Failed to calculate number of pixels in map");
    }

    // Allocate memory for input and output arrays
    std::vector<float> input(nipix);
    std::vector<float> output(nipix);
    std::vector<Nxyz> nxyz;
    std::vector<DefOpt> def;
    std::vector<double> vxy, vz, bias, fwhmxy, fwhmz, squeeze, sqfwhm;
    std::vector<std::vector<double>> wavel;
    std::vector<std::vector<float>> gamma, scale;
    std::vector<Itype> itype;
    double tzero, period, quad, vfine, sfac;

    // Read the map
    if (!read_map(map, input.data(), nxyz, vxy, vz, wavel, gamma, scale, itype,
                  def, bias, fwhmxy, fwhmz, squeeze, sqfwhm,
                  tzero, period, quad, vfine, sfac)) {
        throw std::runtime_error("doppler.comdef: Failed to read map");
    }

    // Process each image in the map
    float* iptr = input.data();
    float* optr = output.data();

    for (size_t nim = 0; nim < nxyz.size(); nim++) {
        size_t npix = nxyz[nim].ntot();

        if (def[nim] == 1) {
            // Uniform default
            double ave = 0.0;
            for (size_t ip = 0; ip < npix; ip++) {
                ave += iptr[ip];
            }
            ave /= npix;
            ave *= bias[nim];
            for (size_t ip = 0; ip < npix; ip++) {
                optr[ip] = static_cast<float>(ave);
            }
        } else if (def[nim] == 2 || def[nim] == 3) {
            // Gaussian defaults
            double fx = fwhmxy[nim] / vxy[nim];
            double fy = fwhmxy[nim] / vxy[nim];
            double fz = (nxyz[nim].nz > 1) ? fwhmz[nim] / vz[nim] : 0.0;

            if (def[nim] == 3 && squeeze[nim] > 0.0) {
                // GAUSS3D with squeeze
                std::vector<float> tptr(npix);

                size_t nx = nxyz[nim].nx;
                size_t ny = nxyz[nim].ny;
                size_t nz = nxyz[nim].nz;
                size_t nxy = nx * ny;
                float sigma = sqfwhm[nim] / EFAC;

                // Loop over x-y
                for (size_t iy = 0; iy < ny; iy++) {
                    for (size_t ix = 0; ix < nx; ix++) {
                        size_t ixyp = iy * nx + ix;

                        // Compute weighted mean vz
                        double vzmean = 0.0, znorm = 0.0;
                        float vslice = -vz[nim] * (nz - 1) / 2;

                        for (size_t iz = 0, ip = ixyp; iz < nz; iz++, ip += nxy, vslice += vz[nim]) {
                            znorm += iptr[ip];
                            vzmean += iptr[ip] * vslice;
                        }
                        vzmean /= znorm;

                        // Construct Gaussian in vz
                        double znorm2 = 0.0;
                        vslice = -vz[nim] * (nz - 1) / 2;

                        for (size_t iz = 0, ip = ixyp; iz < nz; iz++, ip += nxy, vslice += vz[nim]) {
                            tptr[ip] = std::exp(-std::pow(vslice / sigma, 2) / 2.0);
                            znorm2 += tptr[ip];
                        }

                        // Scale
                        double sfac = znorm / znorm2 * squeeze[nim];
                        for (size_t iz = 0, ip = ixyp; iz < nz; iz++, ip += nxy) {
                            tptr[ip] *= sfac;
                        }
                    }
                }

                // Blur in vx-vy
                std::vector<float> ttptr(npix);
                gaussdef(tptr.data(), nxyz[nim], fx, fy, 0.0, ttptr.data());

                // Blur the original over all dimensions
                gaussdef(iptr, nxyz[nim], fx, fy, fz, optr);

                // Combine blurred original and squeezed Gaussian
                for (size_t ip = 0; ip < npix; ip++) {
                    optr[ip] = (1 - squeeze[nim]) * optr[ip] + ttptr[ip];
                }
            } else {
                // GAUSS2D or GAUSS3D without squeeze
                gaussdef(iptr, nxyz[nim], fx, fy, fz, optr);
            }

            // Apply bias
            if (bias[nim] != 1.0) {
                for (size_t ip = 0; ip < npix; ip++) {
                    optr[ip] *= bias[nim];
                }
            }
        }

        // Move to the next image
        iptr += npix;
        optr += npix;
    }

    // Write the modified image back into the map
    update_map(output.data(), map);

    return py::none();
}


py::none doppler_datcom(py::object& data, py::object& map) {
    // Check if data or map is None
    if (data.is_none()) {
        throw std::invalid_argument("doppler.datcom: data is NULL");
    }
    if (map.is_none()) {
        throw std::invalid_argument("doppler.datcom: map is NULL");
    }

    // Read sizes
    size_t nipix, ndpix;
    if (!npix_map(map, nipix)) {
        throw std::runtime_error("doppler.datcom: Failed to calculate number of pixels in map");
    }
    if (!npix_data(data, ndpix)) {
        throw std::runtime_error("doppler.datcom: Failed to calculate number of pixels in data");
    }

    // Allocate memory for Map data
    std::vector<float> image(nipix);
    std::vector<Nxyz> nxyz;
    std::vector<DefOpt> def;
    std::vector<double> vxy, vz, bias, fwhmxy, fwhmz, squeeze, sqfwhm;
    std::vector<std::vector<double>> wavel;
    std::vector<std::vector<float>> gamma, scale;
    std::vector<Itype> itype;
    double tzero, period, quad, vfine, sfac;

    // Read the Map
    if (!read_map(map, image.data(), nxyz, vxy, vz, wavel, gamma, scale, itype,
                  def, bias, fwhmxy, fwhmz, squeeze, sqfwhm,
                  tzero, period, quad, vfine, sfac)) {
        throw std::runtime_error("doppler.datcom: Failed to read map");
    }

    // Allocate memory for Data data
    std::vector<size_t> nwave, nspec;
    std::vector<std::vector<double>> time;
    std::vector<std::vector<float>> expose;
    std::vector<std::vector<int>> nsub;
    std::vector<double> fwhm;
    std::vector<float> flux(ndpix);
    std::vector<float> ferr(ndpix);
    std::vector<double> wave(ndpix);

    // Read the Data
    if (!read_data(data, flux.data(), ferr.data(), wave.data(), nwave,
                   nspec, time, expose, nsub, fwhm)) {
        throw std::runtime_error("doppler.datcom: Failed to read data");
    }

    // Release the GIL for computation
    py::gil_scoped_release release;

    // Perform the data-to-image transformation
    tr(image.data(), nxyz, vxy, vz, wavel, gamma, scale, itype, tzero, period, quad,
       vfine, sfac, flux.data(), wave.data(), nwave, nspec, time, expose, nsub, fwhm);

    // Restore the GIL
    py::gil_scoped_acquire acquire;

    // Write the modified image back into the map
    if (!update_map(image.data(), map)) {
        throw std::runtime_error("doppler.datcom: Failed to update map");
    }

    return py::none();
}


py::none doppler_memit(py::object& map, py::object& data, int niter, float caim, float tlim = 1.e-4, float rmax = 0.2) {
    // Validate arguments
    if (map.is_none()) {
        throw std::invalid_argument("doppler.memit: map is NULL");
    }
    if (data.is_none()) {
        throw std::invalid_argument("doppler.memit: data is NULL");
    }
    if (niter < 1) {
        throw std::invalid_argument("doppler.memit: niter < 1");
    }
    if (caim <= 0.0) {
        throw std::invalid_argument("doppler.memit: caim <= 0");
    }
    if (tlim >= 0.5) {
        throw std::invalid_argument("doppler.memit: tlim >= 0.5");
    }
    if (rmax >= 0.5) {
        throw std::invalid_argument("doppler.memit: rmax >= 0.5");
    }

    // Read the number of image and data pixels
    size_t nipix, ndpix;
    if (!npix_map(map, nipix)) {
        throw std::runtime_error("doppler.memit: Failed to calculate number of pixels in map");
    }
    if (!npix_data(data, ndpix)) {
        throw std::runtime_error("doppler.memit: Failed to calculate number of pixels in data");
    }

    // Compute size needed for MEM buffer
    const size_t MXBUFF = Mem::memsize(nipix, ndpix);

    // Allocate memory for MEM buffer
    std::vector<float> mem_buffer(MXBUFF);
    Mem::Gbl::st = mem_buffer.data();
    std::cerr << "Allocated " << MXBUFF * sizeof(float) << " bytes to MEM buffer" << std::endl;

    // Set pointer offsets
    Mem::memcore(MXBUFF, nipix, ndpix);

    // Read the map into the MEM buffer
    if (!read_map(map, Mem::Gbl::st + Mem::Gbl::kb[0],
                  Dopp::nxyz, Dopp::vxy, Dopp::vz,
                  Dopp::wavel, Dopp::gamma, Dopp::scale,
                  Dopp::itype, Dopp::def, Dopp::bias,
                  Dopp::fwhmxy, Dopp::fwhmz, Dopp::squeeze, Dopp::sqfwhm,
                  Dopp::tzero, Dopp::period, Dopp::quad, Dopp::vfine,
                  Dopp::sfac)) {
        throw std::runtime_error("doppler.memit: Failed to read map");
    }

    // Allocate memory for wavelengths
    std::vector<double> wave(ndpix);
    Dopp::wave = wave.data();

    // Read the data and errors into the MEM buffer
    if (!read_data(data, Mem::Gbl::st + Mem::Gbl::kb[20],
                   Mem::Gbl::st + Mem::Gbl::kb[21], Dopp::wave,
                   Dopp::nwave, Dopp::nspec, Dopp::time,
                   Dopp::expose, Dopp::nsub, Dopp::fwhm)) {
        throw std::runtime_error("doppler.memit: Failed to read data");
    }

    // Process errors
    float* eptr = Mem::Gbl::st + Mem::Gbl::kb[21];
    for (size_t np = 0; np < ndpix; np++) {
        if (eptr[np] > 0.0) {
            eptr[np] = 2.0 / std::pow(eptr[np], 2) / ndpix;
        }
    }

    float c, test, acc = 1.0, cnew, s, rnew, snew, sumf;
    int mode = 10;
    for (size_t nd = 0; nd < Dopp::def.size(); nd++) {
        if (Dopp::def[nd] > 1 || Dopp::bias[nd] != 1.0) {
            mode = 30;
        }
    }

    std::cerr << "mode = " << mode << std::endl;

    // Release the GIL for computation
    py::gil_scoped_release release;

    for (int it = 0; it < niter; it++) {
        std::cerr << "\nIteration " << it + 1 << std::endl;
        if (mode == 30) {
            std::cerr << "Computing default ..." << std::endl;
            float* iptr = Mem::Gbl::st + Mem::Gbl::kb[0];
            float* optr = Mem::Gbl::st + Mem::Gbl::kb[19];
            for (size_t nim = 0; nim < Dopp::nxyz.size(); nim++) {
                size_t npix = Dopp::nxyz[nim].ntot();
                if (Dopp::def[nim] == UNIFORM) {
                    double ave = 0.0;
                    for (size_t np = 0; np < npix; np++) {
                        ave += iptr[np];
                    }
                    ave /= npix;
                    ave *= Dopp::bias[nim];
                    for (size_t np = 0; np < npix; np++) {
                        optr[np] = static_cast<float>(ave);
                    }
                } else if (Dopp::def[nim] == GAUSS2D || Dopp::def[nim] == GAUSS3D) {
                    double fwhmx = Dopp::fwhmxy[nim] / Dopp::vxy[nim];
                    double fwhmy = Dopp::fwhmxy[nim] / Dopp::vxy[nim];
                    double fwhmz = (Dopp::nxyz[nim].nz > 1) ? Dopp::fwhmz[nim] / Dopp::vz[nim] : 0.0;

                    if (Dopp::def[nim] == GAUSS3D && Dopp::squeeze[nim] > 0.0) {
                        std::vector<float> tptr(npix);
                        gaussdef(tptr.data(), Dopp::nxyz[nim], fwhmx, fwhmy, 0.0, optr);
                    } else {
                        gaussdef(iptr, Dopp::nxyz[nim], fwhmx, fwhmy, fwhmz, optr);
                    }
                }
                iptr += npix;
                optr += npix;
            }
        }

        Mem::memprm(mode, 20, caim, rmax, 1.0, acc, c, test, cnew, s, rnew, snew, sumf);
        if (test < tlim && c <= caim) {
            break;
        }
    }

    // Restore the GIL
    py::gil_scoped_acquire acquire;

    // Write the modified image back into the map
    if (!update_map(Mem::Gbl::st + Mem::Gbl::kb[0], map)) {
        throw std::runtime_error("doppler.memit: Failed to update map");
    }

    return py::none();
}
