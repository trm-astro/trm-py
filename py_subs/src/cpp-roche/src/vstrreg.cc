/*

!!begin 
!!title  Computes path of gas stream in velocity coordinates
!!author T.R.Marsh
!!descr  computes path of gas stream in velocity coordinates
!!root   vstrreg
!!index  vstrreg.cc
!!class  Functions
!!css    style.css
!!head1  Computes path of gas stream in velocity coordinates, spaced regularly in radius

!!emph{vstrreg} returns x and y arrays of velocities along stream stepped
regularly in terms of L1 distance.

!!head2 Function call

void Roche::vstrreg(double q, double step, float *vx, float *vy, int n, int type);

!!head2 Arguments

!!table

!!arg{q}{mass ratio = M2/M1. Stream flows from star 2 to 1.}
!!arg{step}{step (fraction of RL1).}
!!arg{vx[]}{array of x velocities returned.}
!!arg{vy[]}{array of y velocities returned.}
!!arg{n}{number of points to compute.}
!!arg{type}{type of velocity, see !!ref{vtrans.html}{vtrans} for supported types.}
!!table

!!end

*/

#include <cstdlib>
#include <cmath>
#include <string>
#include <iostream>
#include "trm/roche.h"
#include "trm/subs.h"

void Roche::vstrreg(double q, double step, float vx[], float vy[], int n, int type){

    const double TLOC =  1.e-8;
    const double RLOC =  1.e-8;
    int np;
    double rl1, tvx, tvy, rend, rnext;
    Subs::Vec3 r, v, rm, vm;
    int decr;

    if(n < 2)
	throw Roche::Roche_Error("Roche::vstrreg: need at least 2 points in vstrreg\n");
  
    rl1 = xl1(q);

    /* Store L1 as first point */

    vtrans(q, type, rl1, 0., 0., 0., tvx, tvy);
    vx[0] = tvx;
    vy[0] = tvy;
    np    = 1;
    rnext = rl1*(1.-step);
    decr  = 1;

    /* Initialise stream */
    strinit(q, r, v);
  
    while(np < n){

	/* Advance one step */
	stradv(q, r, v, rnext, RLOC, 1.e-3);
	vtrans(q, type, r.x(), r.y(), v.x(), v.y(), tvx, tvy);
	vx[np] = tvx;
	vy[np] = tvy;
	np++;
	rnext = decr ? rnext - rl1*step : rnext + rl1*step;

	/* Locate and store next turning point */
	rm = r;
	vm = v;
	strmnx(q, rm, vm, TLOC);
	rend = rm.length();

	/* Loop over all radii wanted before next turning point */
	while(np < n && ((decr && rnext > rend) || (!decr && rnext < rend))){
	    stradv(q, r, v, rnext, RLOC, 1.e-3);
	    vtrans(q, type, r.x(), r.y(), v.x(), v.y(), tvx, tvy);
	    vx[np] = tvx;
	    vy[np] = tvy;
	    np++;
	    rnext = decr ? rnext - rl1*step : rnext + rl1*step;
	}

	/* Change direction of search, and move it to start at turning point */
	rnext = decr ? rnext + rl1*step : rnext - rl1*step;
	r  = rm;
	v  = vm;
	decr = !decr;
    }    
    return;
}







