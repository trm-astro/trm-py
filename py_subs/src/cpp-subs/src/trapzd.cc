/*

!!begin 
!!title  Trapezoidal integration 
!!author T.R.Marsh
!!root   trapzd
!!index  trapzd.cc 
!!descr  trapezoidal integration workhorse 
!!css   style.css
!!head1 trapzd trapezoidal integration workhorse.

!!emph{trapzd} carries out steps of trapezoidal integration at a level
defined by a number n. It should be called in a series of steps
starting at n=0. For n=0, it evaluates the function at
either end. Next it adds a contribution at the half way point.
Then it splits each half by half again etc. For n>0, the
number of function evaluations is 2**(n-1).

See !!ref{qtrap.html}{qtrap} or !!ref{qsimp.html}{qsimp} for the
use of trapzd in integration.
 
!!head2 Function call

template <class X> X trapzd(X (*func)(X x), X a, X b, int n);

!!head2 Argument list

!!table
!!arg{ func}{a 1D function}
!!arg{ a   }{lower limit to integrate from.}
!!arg{ b   }{upper limit to integrate to.}
!!arg{ n   }{step number.}
!!table

!!end

*/

