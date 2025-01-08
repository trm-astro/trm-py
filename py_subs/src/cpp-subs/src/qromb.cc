/*
 
!!begin qromb qromb romberg trapezoidal integrator !!title qromb
!!author T.R.Marsh
!!head1 qromb romberg trapezoidal integrator

qromb carries out trapezoidal integration of a function
until a specified level of accuracy has been reached. It does
this by repeatedly calling !!ref{trapzd.html}{trapzd} until
the integral changes by less than a specified fraction of its
previous value, until a maximum counter has been reached or
if the integral stays at zero twice in a row. It does this
by extrapolating the value and seeing how the extrapolation
changes with successively smaller steps. It is good for
relatively smooth function.

qromb also takes an argument to force it to take a minimum
number of steps to avoid it stopping too early.

See also: !!ref{trapzd.html}{trapzd}, !!ref{qtrap.html}{qtrap}, and
!!ref{qsimp.html}{qsimp}.

!!head2 Function calls

template <class X> X qromb(X (*f)(X x), X a, X b, X eps, int nmin,
int nmax, bool print);

template <class Func, class X> X qromb(const Func& f, X a, X b, X eps, 
int nmin, int nmax, bool print);

!!head2 Argument list

!!table
!!arg{func}{a 1D function or a function object. f(x) is
the value of the function}
!!arg{a}{lower limit to integrate from}
!!arg{b}{upper limit to integrate to}
!!arg{eps}{accuracy}
!!arg{nmin}{number of points to extrapolate (e.g. 5)}
!!arg{nmax}{maximum sub-division factor (should be > nmin)}
!!arg{print}{print diagnostic information}
!!table

!!end

*/


