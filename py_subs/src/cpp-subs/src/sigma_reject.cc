#include "trm/subs.h"

/** Carries out sigma clipping. There are two possible methods: one in which pixels are rejected one
 * at a time, worst first, with the mean and rms recomputed after each rejection, or all pixels out of 
 * range can be rejected at one go, the mean and rms computed, and then another go, etc. The first method is safe
 * but potentially slow as it has a section which scales with N**2. The more cavalier method is faster but
 * can lead to the rejection of good data.
 * \param data    array of values
 * \param n       number of values
 * \param thresh  threshold number of sigma to reject at
 * \param careful controls type of clipping. If true, then pixels are rejected one by one, if false they are
 * rejected in whole swathes per cycle. 
 * \param rawmean mean of all values
 * \param rawrms  RMS of all values (sample variance, N-1)
 * \param mean    mean after clipping
 * \param rms     rms after clipping
 * \param nrej    number of rejects
 */

void Subs::sigma_reject(const float* data, int n, float thresh, bool careful,
			double& rawmean, double& rawrms, double& mean, double& rms, int& nrej){

  static Buffer1D<bool> ok;

  if(n < 1){
    rawrms = rms = 0.;
    return;
  }
  nrej = 0;
  int ntot = n;
  ok.resize(n); // mask used during rejection to keep track of bad pixels
  ok = true;
  double s1=0.;
  for(int i=0; i<n; i++)
    s1   += data[i];

  mean = rawmean = s1/n;
  if(n < 2){
    rawrms = rms = 0.;
    return;
  }

  double s2 = 0.;
  for(int i=0; i<n; i++)
    s2  += sqr(data[i]-rawmean);

  rms = rawrms = sqrt(s2/(n-1));
  if(n < 3) return;

  // Now rejection of values
  if(careful){

    // Only one per cycle (the worst)
    float resid;
    for(;;){
      float worst = -1.;
      int iworst = 0;
      for(int i=0; i<n; i++){
	if(ok[i] && (resid = fabs(data[i]-mean)) > worst){
	  worst  = resid;
	  iworst = i;
	}
      }

      if(worst > thresh*rms){

	// Now must update sums. Save some time by not recomputing them but
	// just altering them while being careful to avoid round off.
	s1 -= data[iworst];
	double newmean = s1/(ntot-1);
	s2 -= (data[iworst]-mean)*(data[iworst]-newmean);

	ok[iworst] = false;
	mean = newmean;
	nrej++;
	ntot--;
	if(ntot > 1){
	  rms  = sqrt(s2/(ntot-1));
	}else{
	  rms  = 0.;
	  break;
	}
      }else{
	break;
      }
    }

  }else{

    // Reject all points above threshold

    for(;;){
      int nr = 0;
      s1=0.;
      float cut = rms*thresh;
      for(int i=0; i<n; i++){
	if(ok[i]){
	  if(fabs(data[i]-mean) > cut){
	    ok[i] = false;
	    nr++;
	    ntot--;
	  }else{
	    s1 += data[i];
	  }
	}
      }

      if(ntot == 0){
	nrej += nr;
	return;
      }

      mean =  s1/ntot;
      if(ntot < 2){
	rms = 0.;
	return;
      }

      s2 = 0.;
      for(int i=0; i<n; i++)
	if(ok[i]) s2  += sqr(data[i]-mean);

      rms = sqrt(s2/(ntot-1));
      if(ntot < 3) return;

      if(nr == 0) break;
      nrej += nr;
    }
  }
}
