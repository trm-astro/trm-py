#include <iostream>
#include "trm/binary_star.h"

std::ostream& Subs::operator<<(std::ostream& ost, const Binary& binary){
  ost << (Star&)binary << "\n" << (Ephem&)binary;
  return ost;
}

void Subs::Binary::set(const Star& star, const Ephem& eph){
  *this = Binary(star,eph);
}
