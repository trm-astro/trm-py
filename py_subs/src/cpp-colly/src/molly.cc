#include "molly.h"

std::ostream& operator<<(std::ostream& ost, const Molly_Error& obj){
  ost << obj.message;
  return ost;
}
