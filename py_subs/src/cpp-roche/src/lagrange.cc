/*

!!begin
!!title  Lagrange
!!author T.R. Marsh
!!descr  prints out Lagrangian points
!!index  lagrange
!!root   lagrange
!!class  Programs
!!css    style.css
!!head1  Program to print out Lagrangian point distances

!!emph{lagrange} prints out the positions of the first three Lagrangian
points (units of separation).

!!head2 Arguments

!!table

!!arg{q}{mass ratio = M2/M1}

!!table

Distance referred to M1. C of M of M2 = 1.

!!end

*/

#include <cstdlib>
#include <iostream>
#include "trm/roche.h"

int main (int argc, char *argv[]){
  if(argc < 2)
    std::cerr << "usage: lagrange q\n";
  double q = atof(argv[1]);
  std::cout << "xL1 = " << Roche::xl1(q) 
	    << ", xL2 = " << Roche::xl2(q) 
	    << ", xL3 = " << Roche::xl3(q) << std::endl;
  exit(EXIT_SUCCESS);
}


