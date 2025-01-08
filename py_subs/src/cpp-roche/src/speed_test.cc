#include <cstdlib>
#include <iostream>
#include "trm/subs.h"
#include "trm/vec3.h"
#include "trm/roche.h"

int main(){

    Subs::Vec3 r;
    double q = 0.2;
    double iangle = 87.;
    int N = 200000;
    Subs::INT4 seed = -12345;
    double ffac = 0.9, ingress, egress;
    for(int i=0; i<N; i++){
	r.x() = -0.5 + Subs::ran1(seed);
	r.y() = -0.5 + Subs::ran1(seed);
	r.z() = -0.5 + Subs::ran1(seed);
	bool ecl = Roche::ingress_egress(q, Roche::SECONDARY, 1.0, ffac, iangle, 3.e-8, r,  ingress, egress);
	//	std::cout << r << " " << ecl << " " << ingress << " " << egress << std::endl;
    }
    
}
