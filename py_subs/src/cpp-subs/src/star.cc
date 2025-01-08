#include <iostream>
#include "trm/star.h"

/**
 * ASCII input of star data
 */
std::istream& Subs::operator>>(std::istream& ist, Star& star){

  if(!ist) return ist;

  do{
    if(!getline(ist,star.nam)) return ist;
  }
  while(star.nam.find_first_not_of(" ") == std::string::npos || star.nam[0] == '#');

  std::string::size_type first = 0;
  if(star.nam[first] == ' ') first = star.nam.find_first_not_of(" ");
  std::string::size_type last = star.nam.length()-1;
  if(star.nam[last] == ' ') last = star.nam.find_last_not_of(" ");
  star.nam = star.nam.substr(first,last-first+1);

  ist >> (Position&)star;
  return ist;
}

std::ostream& Subs::operator<<(std::ostream& ost, const Star& star){
  ost << "   Star: " << star.name() << "\n" << (Position&)star;
  return ost;
}
