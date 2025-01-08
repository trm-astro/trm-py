#include "trm/array1d.h"

// Specialisations for float vectors

template <>
void pgline(const Array1D<float>& vec){
  if(vec.size()){
    float x[vec.size()];
    for(size_t i=0; i<vec.size(); i++) x[i] = float(i+1);
    cpgline(vec.size(),x,vec);
  }
}

template <> 
void pgline(const Array1D<float>& vec1, const Array1D<float>& vec2){
  if(vec1.size()){
    if(vec1.size() == vec2.size()){
      cpgline(vec1.size(),vec1,vec2);
    }else{
      throw Array1D_Error::Size_Mismatch();
    }
  }
}


