/*
!!begin
!!title  Velocity transformations
!!author T.R.Marsh
!!date   13 October 2000
!!root   vtrans
!!index  vtrans.cc
!!descr  Computes velocity transforms for rotating frame velocities
!!css   style.css
!!class  Functions
!!head1 vtrans - computes velocity transforms

vtrans computes two velocity transforms, (1) a straight transform 
from rotating to inertial frame and (2) an inertial frame velocity 
in the disc.

!!head2 Function call 

void Roche::vtrans(double q, int type, double x, double y, double vx, double vy,
            double &tvx, double &tvy){

!!head2 Arguments

!!table
!!arg{q}{ mass ratio = M2/M1}
!!arg{type}{1 for rotating->inertial, 2 for position to disc, 3 for rotating.}
!!arg{x}{x position (units of separation)}
!!arg{y}{y position (units of separation)}
!!arg{vx}{x velocity (omega*a = 1 units)}
!!arg{vy}{y velocity (omega*a = 1 units)}
!!arg{tvx}{transformed x velocity}
!!arg{tvy}{transformed y velocity}
!!table

When translating to inertial, the accretor velocity is added.
If you want the velocity relative to this you must add
mu = q/(1+q) to tvy before using it.

!!end

*/

#include <cstdlib>
#include <cmath>
#include <string>
#include <iostream>
#include "trm/roche.h"

void Roche::vtrans(double q, int type, double x, double y, double vx, double vy, double &tvx, double &tvy){

    double mu, rad, vkep;
    
    mu = q/(1.0+q);
    switch(type){
	case 1:
	    tvx = vx - y;
	    tvy = vy + x - mu;
	    break;
	case 2:
	    rad  = sqrt(x*x+y*y);
	    vkep = 1.0/sqrt((1.0+q)*rad);
	    tvx = -vkep*y/rad;
	    tvy =  vkep*x/rad-mu;
	    break;
	case 3:
	    tvx = vx;
	    tvy = vy;
	    break;
	default:
	    throw Roche_Error("Error in vtrans: did not recognize type = " + Subs::str(type) + 
			      ". Only 1 or 2 supported.");
    }
    return;
}
