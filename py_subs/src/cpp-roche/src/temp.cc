#include <stdio.h>
#include "trm/roche.h"

int main(){
  double q = 0.7;
  int i;
  float vx[100], vy[100];

  vstrreg(q,vx,vy,100,1,0.01);

  for(i=0; i<100; i++){
    printf("%d %f %f\n",i,vx[i],vy[i]);
  }
}
