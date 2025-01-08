#include <iostream>

int main(){

  unsigned char c1 = 1, c2 = 2;
  unsigned short int usi;

  usi = ((0xffff & c1) << 8) | c2;
  std::cerr << "one way   = " << usi << std::endl; 

  usi = ((0xffff & c2) << 8) | c1;
  std::cerr << "the other = " << usi << std::endl; 

  usi = (c2 << 8) | c1;
  std::cerr << "an  other = " << usi << std::endl; 

  std::cerr << (c2 << 8) << std::endl; 

}
