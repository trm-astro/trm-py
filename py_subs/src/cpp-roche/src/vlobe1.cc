/*

!!begin 
!!title  Computes primary's Roche lobe in velocity coordinates 
!!author T.R.Marsh
!!descr  computes primary's Roche lobe in velocity coordinates 
!!root   vlobe1
!!index  vlobe1.cc
!!class  Functions
!!css    style.css
!!head1  vlobe1 - computes primary star's Roche lobe in velocity

!!emph{vlobe1} returns arrays vx and vy for plotting an equatorial section
of the Roche lobe of the primary star in a binary of mass ratio q = M2/M1
in Doppler coordinates. The arrays start and end at the inner Lagrangian 
point and march around uniformly in azimuth looking from the centre of 
mass of the primary star. n is the number of points and must be at least 3. 

See also: !!ref{vlobe2.html}{vlobe2}, !!ref{lobe1.html}{lobe1},
!!ref{lobe2.html}{lobe2}

!!head2 Function call 

void Roche::vlobe1(double q, float *vx, float *vy, int n)

!!head2 Arguments


!!table
!!arg{q}{mass ratio = M2/M1}
!!arg{vx[n]}{array of vx values}
!!arg{vy[n]}{array of vy values}
!!arg{n}{number of x and y values}
!!table

!!end

*/

#include <stdlib.h>
#include <math.h>
#include "trm/roche.h"

void Roche::vlobe1(double q, float *vx, float *vy, int n){

  int i;
  float tvx, tvy, mu;

  // Call lobe1 then transform appropriately

  lobe1(q, vx, vy, n);

  mu = q/(1.0+q);
  for(i=0;i<n;i++){
    tvx = - vy[i];
    tvy =   vx[i] - mu;
    vx[i] = tvx;
    vy[i] = tvy;
  }
  return;
}

