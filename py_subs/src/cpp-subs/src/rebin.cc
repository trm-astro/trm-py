#include "trm/subs.h"
#include "trm/poly.h"

/** General purpose rebin function, equivalent to F77 function of same name with
 * a few modifications.
 * \param input the array to be rebinned.
 * \param infirst the first pixel to bother with in the input array (minimum 0)
 * \param inlast one more than the last pixel to deal with in the input array (maximum = number of pixels)
 * \param inpoly polynomial describing scale of input data
 * \param output the array to hold the output. Will be set to zero where there is no overlap.
 * \param outfirst the first pixel to bother with in the output array (minimum 0)
 * \param outlast one more than the last pixel to deal with in the output array (maximum = number of pixels)
 * \param outpoly polynomial describing the scale of the output array.
 * \param mode AVERAGE to average values, INTEGRATE to sum.
 * \param type LINEAR for linear rebin
 */

void Subs::rebin(const Buffer1D<double>& input, int infirst, int inlast, const Poly& inpoly, 
		 Buffer1D<double>& output, int outfirst, int outlast, const Poly& outpoly, 
		 REBIN_MODE mode, REBIN_TYPE type){
  

} 


/** General purpose rebin function, equivalent to F77 function of same name with
 * a few modifications. This one rebins a data array and its variances at the same time
 * \param indat the data array to be rebinned.
 * \param invar the variances array to be rebinned.
 * \param infirst the first pixel to bother with in the input array (minimum 0)
 * \param inlast one more than the last pixel to deal with in the input array (maximum = number of pixels)
 * \param inpoly polynomial describing scale of input data
 * \param outdat the array to hold the rebinned data. Will be set to zero where there is no overlap.
 * \param outvar the array to hold the rebinned variances. Will be set to zero where there is no overlap.
 * \param outfirst the first pixel to bother with in the output array (minimum 0)
 * \param outlast one more than the last pixel to deal with in the output array (maximum = number of pixels)
 * \param outpoly polynomial describing the scale of the output array.
 * \param mode AVERAGE to average values, INTEGRATE to sum.
 */

void Subs::rebin(const Buffer1D<double>& indat, const Buffer1D<float>& invar, int infirst, int inlast, const Poly& inpoly, 
		 Buffer1D<double>& outdat, Buffer1D<float>& outvar, int outfirst, int outlast, const Poly& outpoly, 
		 REBIN_MODE mode, REBIN_TYPE type){
  
    if(indat.size() != invar.size() || outdat.size() != outvar.size())
	throw Subs_Error("rebin: conflicting numbers of data & variance pixels: " + Subs::str(indat.size()) + " " + 
			 Subs::str(invar.size()) + " " + Subs::str(outdat.size()) + " " + Subs::str(outvar.size()));
  
    // Check validity of input ranges
    if(infirst < 0 || inlast > indat.size() || outfirst < 0 || outlast > outdat.size() || infirst >= inlast || outfirst >= outlast)
	throw Subs_Error("rebin: invalid pixels ranges " + Subs::str(infirst) + " " + Subs::str(inlast) + " " + 
			 Subs::str(indat.size()) + " " + Subs::str(outfirst) + " " + Subs::str(outlast) + " " + Subs::str(outdat.size()));
  
    if(inpoly.size() == 0 || outpoly.size() == 0)
	throw Subs_Error("rebin: invalid input and/or output polys, npoly = " + Subs::str(inpoly.size()) + " " + Subs::str(outpoly.size())); 
  
    // First compute the extreme range of the indat array. The scale is assumed to be monotonic.
    // Variables starting 'in' refer to the indat array, 'pix' means they are pixel values.
    double in_pix_leftmost   = infirst-0.5;
    double in_pix_rightmost  = inlast-0.5;
    double in_leftmost       = inpoly.get_value(in_pix_leftmost);
    double in_rightmost      = inpoly.get_value(in_pix_rightmost);
    double in_min            = std::min(in_leftmost, in_rightmost);
    double in_max            = std::max(in_leftmost, in_rightmost);
    double in_pix_right_old  = 0.;
    bool increase = outpoly.get_value((outfirst+outlast-1)/2.) > 0.;
    for(int iout=0; iout<outdat.size(); iout++){
    
	if(iout >= outfirst && iout < outlast){
      
	    // OK we are on a valid output pixel. Compute its boundary values.
	    double out_left   = outpoly.get_value(iout-0.5);
	    double out_right  = outpoly.get_value(iout+0.5);
	    if(out_left > out_right) std::swap(out_left, out_right);
      
	    // Now check that there is any overlap at all with the indat array
	    if(out_left > in_max || out_right < in_min){
		outdat[iout] = 0.;
		outvar[iout] = 0.;
	    }else{
	
		// Yes there is overlap, we now need to compute the pixel boundaries in
		// the indat array equivalent to those from the output array.
		double in_pix_left;
		if(increase && iout > outfirst){
		    in_pix_left = in_pix_right_old;
		}else{
		    if(in_leftmost < in_rightmost && out_left < in_leftmost){
			in_pix_left = in_pix_leftmost;
		    }else if(in_leftmost > in_rightmost && out_left < in_rightmost){
			in_pix_left = in_pix_rightmost;
		    }else{
			// Linearly interpolate to get a start then refine
			in_pix_left = linterp(in_leftmost, in_pix_leftmost, in_rightmost, in_pix_rightmost, out_left);
			in_pix_left = inpoly.get_x(out_left, in_pix_left, 1.e-4);
		    }
		}
		double in_pix_right;
		if(!increase && iout > outfirst){
		    in_pix_right = in_pix_right_old;
		}else{
		    if(in_leftmost < in_rightmost && out_right > in_rightmost){
			in_pix_right = in_pix_rightmost;
		    }else if(in_leftmost > in_rightmost && out_right > in_leftmost){
			in_pix_right = in_pix_leftmost;
		    }else{
			// Linearly interpolate to get a start then refine
			in_pix_right = linterp(in_leftmost, in_pix_leftmost, in_rightmost, in_pix_rightmost, out_right);
			in_pix_right = inpoly.get_x(out_right, in_pix_right, 1.e-4);
		    }
		}
	
		// Save value for (true) right-hand edge of output pixel to save time in steps after the first
		if(increase)
		    in_pix_right_old = in_pix_right;
		else
		    in_pix_right_old = in_pix_left;
	
		// Finally we have the boundaries of the output pixel in terms of the indat pixels
		if(in_pix_left > in_pix_right) std::swap(in_pix_left, in_pix_right);
	
		// Therefore, we now sum over these pixels starting with partial ones.
		double sumdat = 0., sumvar = 0.;
		int first = int(nint(in_pix_left));
		int last  = int(nint(in_pix_right));
		double part_first = first+0.5-in_pix_left;
		double part_last  = in_pix_right-last+0.5;
	
		sumdat += part_first*indat[first];
		sumdat += part_last*indat[last];
	
		sumvar += sqr(part_first)*invar[first];
		sumvar += sqr(part_last)*invar[last];
	
		for(int j=first+1; j<last; j++){
		    sumdat += indat[j];
		    sumvar += invar[j];
		}
	
		// Normalise by number of input pixels or not as the case may be
		if(mode == AVERAGE){
		    double lpix = in_pix_right-in_pix_left; 
		    sumdat /= lpix;
		    sumvar /= sqr(lpix);
		}
		outdat[iout] = sumdat;
		outvar[iout] = sumvar;
	    }
	}else{
	    outdat[iout] = 0.;
	    outvar[iout] = 0.;
	}
    }
} 



