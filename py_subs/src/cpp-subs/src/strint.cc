//
// routine to convert an integer into a string padded initially
// with zeroes.
//
// num    -- the number (should be >=0)
// nd     -- number of digits to write it into
// intstr -- character pointer (must have enough memory allocated)

#include <cstdlib>
#include <iostream>

void strint(const unsigned int num, const unsigned int nd, char* intstr){

  unsigned int tnum = 10, md=1;
  while(tnum <= num){
    tnum *= 10;
    md++;
  }
  if(md > nd){
    std::cerr << "num too large in strint!!\n";
    exit(EXIT_FAILURE);
  }
  tnum /= 10;
  unsigned int ld, nc=0, snum=num;
  for(unsigned int i = 0; i < nd-md; i++){
    intstr[nc++] = '0';
  }
  while(tnum > 0){
    ld = snum / tnum;
    switch(ld){
    case 0:
      intstr[nc] = '0';
      break;
    case 1:
      intstr[nc] = '1';
      break;
    case 2:
      intstr[nc] = '2';
      break;
    case 3:
      intstr[nc] = '3';
      break;
    case 4:
      intstr[nc] = '4';
      break;
    case 5:
      intstr[nc] = '5';
      break;
    case 6:
      intstr[nc] = '6';
      break;
    case 7:
      intstr[nc] = '7';
      break;
    case 8:
      intstr[nc] = '8';
      break;
    case 9:
      intstr[nc] = '9';
      break;
    default:
      std::cerr << "Disaster in strint!!" << std::endl;
      std::cerr << "ld = " << ld << ", snum = " << snum << 
	", tnum = " << tnum << std::endl;
      exit(EXIT_FAILURE);
    }
    snum -= ld*tnum;
    tnum /= 10;
    nc++;
  }
  intstr[nc] = 0;
}
