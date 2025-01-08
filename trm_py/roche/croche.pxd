from libcpp.string cimport string
from libcpp cimport bool

from csubs cimport Vec3

cdef extern from "trm/roche.h" namespace "Roche":
    enum STAR:
        PRIMARY, SECONDARY
        
    #Computes the position of a point on the Roche distorted surface
    void roche_face "Roche::face" \
        (double q, STAR star, double spin, \
        Vec3& dirn, double rref, double pref, double acc,\
        Vec3& pvec, Vec3& dvec, double& r, double& g) except+
    
    #Advances a free particle orbit in a Roch potential   
    double roche_stradv "Roche::stradv" \
        (double q, Vec3&, Vec3&, double, double, double) except+
    
    #Initialises gas stream    
    void roche_strinit "Roche::strinit" (double q, Vec3& r, Vec3& v)

    #Computes radius of reference sphere and Roche potential for
    #Roche-distorted star
    void roche_ref_sphere "Roche::ref_sphere" \
        (double q, STAR star, double spin, double ffac,\
        double& rref, double& pref)
    
    #Returns the earth vector for a given inclination angle and orbital
    #phase
    Vec3 roche_set_earth "Roche::set_earth" (double iangle, double phase)
    Vec3 roche_set_earth "Roche::set_earth" \
        (double cosi, double sini, double phase)
        
    #Is a point in eclipse or not?
    bint roche_fblink "Roche::fblink" \
        (double q, STAR star, double spin, double ffac, double acc,\
        const Vec3& earth, const Vec3& p)
    
    #Computes ingress & egress phases in Roche geometry    
    bint roche_ingress_egress "Roche::ingress_egress" \
        (double q, STAR star, double spin, double ffac, \
        double iangle, double delta, const Vec3& r,double& \
        ingress, double& egress)

    #Calculates primary star's Roche lobe in orbital plane
    void roche_lobe1 "Roche::lobe1" \
        (double q, float *x, float *y, int n)    
    #Calculates secondary star's Roche lobe in orbital plane, 
    void roche_lobe2 "Roche::lobe2" \
        (double q, float *x, float *y, int n) 
    
    #Calculates arbitrary Roche equipotential around the Primary in the
    #orbital plane
    void roche_flobe1 "Roche::flobe1" \
        (double q, float *x, float *y, int n, double pot)
        
    #Computes L1 point distance from primary
    double roche_xl1 "Roche::xl1" (double q)
    #Computes L2 point distance from primary
    double roche_xl2 "Roche::xl2" (double q)
    #Computes L3 point distance from primary
    double roche_xl3 "Roche::xl3" (double q)
    #Computes L1 point distance from primary when primary is asynchronous
    double roche_xl11 "Roche::xl11" (double q, double spin)
    #Computes L1 point distance from primary when secondary is asynchronous
    double roche_xl12 "Roche::xl12" (double q, double spin)
    
    #Computes Roche potential
    double roche_rpot  "Roche::rpot"  (double q, Vec3& p)
    #Computes asynchronous Roche potential, star 1
    double roche_rpot1 "Roche::rpot1" (double q, double spin, Vec3& p)
    #Computes asynchronous Roche potential, star 2
    double roche_rpot2 "Roche::rpot2" (double q, double spin, Vec3& p)
    
    #Computes derivative of Roche potential
    Vec3 roche_drpot  "Roche::drpot"  (double q, Vec3& p)
    #Computes derivative of asynchronous Roche potential, star 1
    Vec3 roche_drpot1 "Roche::drpot1" (double q, double spin, Vec3& p)
    #Computes derivative of asynchronous Roche potential, star 2
    Vec3 roche_drpot2 "Roche::drpot2" (double q, double spin, Vec3& p)
    
    #Computes shadows of Roche lobe in equatorial plane
    void roche_shadow "Roche::roche_shadow" \
        (double q, double iangle, double phi, double dist, \
        double acc, float x[], float y[], bool shade[], int n);

    #Calculates a gas stream
    void roche_streamr "Roche::streamr" \
        (double q, double rad, float *x, float *y, int n) except +
    #Calculates a gas stream
    void roche_stream  "Roche::stream" \
        (double q, double step, float *x, float *y, int n) except +

    #Integrates gas stream    
    void roche_gsint "Roche::gsint" \
        (double q, Vec3 &r, Vec3 &v, double ttry, \
        double &tdid, double &tnext, double &time, double eps) except +
     
    #Calculates a gas stream in velocity coords   
    void roche_vtrans "Roche::vtrans" \
        (double q, int type, double x, double y, \
        double vx, double vy, double &tvx, double &tvy) except +
    
    #Locate next point at minimum or maximum distance from the accretor
    void roche_strmnx "Roche::strmnx" \
        (double q, Vec3 &r, Vec3 &v, double acc)
    
    #Calculates primary star's Roche lobe in orbital plane, velocity
    void roche_vlobe1 "Roche::vlobe1" \
        (double q, float *vx, float *vy, int n)

    #Calculates secondary star's Roche lobe in orbital plane, velocity
    void roche_vlobe2 "Roche::vlobe2" \
        (double q, float *vx, float *vy, int n)
       
    #Gas stream at regular radius steps in velocity coords
    void roche_vstrreg "Roche::vstrreg" \
        (double q, double step, float *vx, float *vy, int n, int type)
        
    #Works out Jacobi constant
    double roche_jacobi "Roche::jacobi" \
        (double q, Vec3& r, Vec3& v)
        
    cdef cppclass Roche_Error:
        Roche_Error()
        Roche_Error(string)
        
