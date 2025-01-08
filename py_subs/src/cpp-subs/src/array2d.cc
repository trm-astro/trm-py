#include "trm/array2d.h"

namespace Array2D_Error {
  std::ostream& operator<<(std::ostream& ost, const Array2D_Error::Error& bd){
    ost << bd.s;
    return ost;
  }
}
