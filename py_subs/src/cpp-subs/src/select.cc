/*

!!begin 
!!title  select
!!author T.R. Marsh 
!!date   18 May 2001
!!descr  finds the k-th smallest element of an array
!!css   style.css
!!root   select
!!class  Functions
!!head1  select the k-th smallest element of an array

!!ref{select} finds the k-th smallest element of an array. Although
this can be done by sorting, it turns out that there is a speed
advantage in just going for the one element of interest.

!!head2 Function call

template <class T> 
T select(T arr[], unsigned long int n, unsigned long int k)

!!head2 Arguments

!!table
!!arg{arr}{Array to select from. !!emph{It is returned in a scrambled order!}}
!!arg{n}{Number of elements}
!!arg{n}{The element to choose, starting with k=0 as the smallest,
and running up to n-1 at the largest.}
!!table

!!head2 Related commands

!!ref{quicksort.html}{quicksort}, !!ref{heapsort.html}{heapsort}.

!!head2 Timing

Takes 0.4 secs to select two percentiles of a 2 million array
of floats (300MHz Pentium II, g++ -O). Linear scaling with n 
plus small overhead.

!!end

*/






