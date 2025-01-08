// !!begin qtrap qtrap simple trapezoidal integrator !!title qtrap
// !!author T.R.Marsh
// !!head1 qtrap simple trapezoidal integrator
//
// qtrap carries out trapezoidal integration of a function
// until a specified level of accuracy has been reached. It does
// this by repeatedly calling !!ref{trapzd.html}{trapzd} until
// the integral changes by less than a specified fraction of its
// previous value or until a maximum counter has been reached.
//
// !!head2 Function call
//
// template <class X> X qtrap(X (*func)(X x), X a, X b, X eps, int nmax);
//
// !!head2 Argument list
//
// !!table
// !!arg func !! a 1D function 
// !!arg a    !! lower limit to integrate from
// !!arg b    !! upper limit to integrate to
// !!arg eps  !! accuracy
// !!arg nmax !! maximum sub-division factor. 
// !!table
//
// !!end















