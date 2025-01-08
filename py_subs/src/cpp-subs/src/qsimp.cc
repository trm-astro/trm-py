/*
 
!!begin 
!!title  Simpson integrator
!!author T.R.Marsh
!!descr  simpson trapezoidal integration
!!css   style.css
!!root   qsimp
!!index  qsimp.cc
!!head1  simpson trapezoidal integration

!!emph{qsimp} carries out trapezoidal integration of a function
until a specified level of accuracy has been reached. It does
this by repeatedly calling !!ref{trapzd.html}{trapzd} until
the integral changes by less than a specified fraction of its
previous value, until a maximum counter has been reached or
if the integral stays at zero twice in a row.

!!emph{qsimp} also takes an argument to force it to take a minimum
number of steps to avoid it stopping too early.

See also: !!ref{trapzd.html}{trapzd}, !!ref{qtrap.html}{qtrap}

!!head2 Function call

template <class X> X qsimp(X (*func)(X x), X a, X b, X eps, 
int nmin, int nmax, bool print);

!!head2 Argument list

!!table
!!arg{func  }{ a 1D function }
!!arg{a     }{ lower limit to integrate from}
!!arg{b     }{ upper limit to integrate to}
!!arg{eps   }{ accuracy}
!!arg{nmin  }{ mininum number to reach when calling trapzd}
!!arg{nmax  }{ maximum sub-division factor. }
!!arg{print }{ print diagnostic info }
!!table

!!end



