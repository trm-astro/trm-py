#include <iostream>
int main(){

#ifdef __BIG_ENDIAN__ 
  std::cout << "big-endian defined" << std::endl;
#else
  std::cout << "big-endian undefined" << std::endl;
#endif

#ifdef __LITTLE_ENDIAN__ 
  std::cout << "little-endian defined" << std::endl;
#else
  std::cout << "little-endian undefined" << std::endl;
#endif

}
