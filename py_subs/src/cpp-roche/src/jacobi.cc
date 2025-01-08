
/*

!!begin
!!title  Jacobi constant
!!author T.R. Marsh
!!descr  returns Jacobi constant
!!index  jacobi.cc
!!root   jacobi
!!class  Functions
!!css    style.css
!!head1  returns Jacobi constant

!!emph{jacobi} returns the Jacobi constant corresponding to a particular
position and velocity.

!!head2 Function call

double Roche::jacobi(double q, coord r, coord v)

!!head2 Arguments

!!table

!!arg{q}{mass ratio}

!!arg{r, v}{position and velocity}

!!table


!!end

*/

#include <math.h>
#include "trm/roche.h"
#include "trm/subs.h"

double Roche::jacobi(double q, const Subs::Vec3& r, const Subs::Vec3& v){
    double f1, f2, yzsq;
    
    f1 = 1.0/(1.0+q);
    f2 = f1*q;
    
    yzsq = Subs::sqr(r.y()) + Subs::sqr(r.z());
    double jac = (Subs::sqr(v.x())+Subs::sqr(v.y())+Subs::sqr(v.z())-Subs::sqr(r.y())-Subs::sqr(r.x()-f2))/2.0 - 
	f1/sqrt(Subs::sqr(r.x())+yzsq) - f2/sqrt(Subs::sqr(r.x()-1.0)+yzsq);
    return jac;
}
