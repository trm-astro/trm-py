/*

!!sphinx

weekday -- returns day of the week
==================================

*weekday* returns the day of the week for a given date, e.g.::

  weekday "25 Dec 2000"

!!sphinx

*/

#include <string>
#include "trm/subs.h"
#include "trm/date.h"

int main(int argc, char* argv[]){
  try{
    if(argc != 2) throw std::string("usage: weekday date");
    Subs::Date date(argv[1]);
    std::cout << date << " ---> " << date.day_of_week() << std::endl;
  }
  catch(const std::string& message){
    std::cerr << message << std::endl;
  }
}
