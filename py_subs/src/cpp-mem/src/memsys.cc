//
// C++ maximum entropy code.
//
// T.R.Marsh July 2000
//

#include <cstdlib>
#include <cmath>
#include <iostream>
#include "trm/memsys.h"

namespace Mem {

    // Global variables; defined here, declared in the user header
    // file memsys.h to prevent multiple declarations.

    namespace Gbl {
        int nj = 0, mj = 0, nk = 0, mk = 0;
        int ka[40] = {0}, kb[40] = {0}, kc[40] = {0}, kd[40] = {0};
        int l0 = 0, l1 = 0, m0 = 0, m10 = 0, m11 = 0, m20 = 0, m21 = 0, m3 = 0;
        float *st = nullptr;
        char pr[41] = {};
    }
    // Function declarations.

    void meinit();

    void memex();

    void mema(const float acc , float &c);

    void mema0(const int mk, const float d[], const float f[], const int m2,
               const float e[], const float acc, float &x);

    void mema1(const int mk, const float d[], const float f[], const int m2,
               const float e[], const float acc, const float x, float gc[],
               float &c);

    void memb(const float def, float sdd[6][6], float &s, float &sdef,
              float &test);

    void memb0(const int mj, const float f[], const int m3,
               const float wtfile[], float &sumf, float &sumflf);

    void memb1(const int mj, const float f[], const float gc[], const int m3,
               const float wtfile[], const int m1, const float amodel[],
               const float def, float fgs[], float fgc[],
               float &sdef, float &s11, float &s12,
               float &s13, float &s22, float &s23,
               float &s33);

    void memc (const float acc, float cdd[6][6]);

    void memc1(const int mk, const float a1[], const float a2[],
               const float a3[], const int m2, const float e[],
               const float acc, float &c11, float &c12, float &c13,
               float &c22, float &c23, float &c33);

    void memcc(const float acc, const float w[3]);

    void memcc1(const int mk, const float a1[], const float a2[],
                const float a3[], const int m2, const float e[],
                const float acc, const float w1, const float w2,
                const float w3, float b[]);

    void memd(float sdd[6][6]);

    void memd1(const int mj, const float f[], const float a2[],
               const float a3[], const float b[], const int m3,
               const float wtfile[], float a4[], float &s14,
               float &s24, float &s34, float &s44);

    void meme (const float acc, float cdd[6][6]);

    void meme1(const int mk, const float a1[], const float a2[],
               const float a3[], const float a4[], const int m2,
               const float e[], const float acc, float &c14,
               float &c24, float &c34, float &c44);

    void memee(const float acc, const float w[4]);

    void memee1(const int mk, const float a1[], const float a2[],
                const float a3[], const float a4[], const int m2,
                const float e[], const float acc, float b[],
                const float w1, const float w2, const float w3,
                const float w4);

    void memf(const float w[4], float &sumf);

    void memf1(const int mj, const float a1[], const float a2[],
               const float a3[], const float a4[], const int m3,
               const float wtfile[], const float w1,
               const float w2, const float w3, const float w4,
               float f[], float &sum);

    void memf2(const int mj, const float scale, const float old[], float f[]);

    void memff(const float w[5], float &sumf);

    void memff1(const int mj, const float a1[], const float a2[],
                const float a3[], const float a4[], const float a5[],
                const int m3, const float wtfile[], const float w1,
                const float w2, const float w3, const float w4,
                const float w5, float f[], float &sum);

    void memff2(const int mj, const float scale, const float old[], float f[]);

    void memj(const int m, const int n);

    void memj1(const int m, const float a[], float b[]);

    void memk(const int m, const int n);

    void memk1(const int m, const float a[], float b[]);

    void meml(const int nsrch, const float sdda[6][6], const float cdda[6][6],
              const float sa, const float ca, const float aima,
              const float rmaxa, const float sdefa, float &snewa,
              float &cnewa, float &rnewa, float wa[6]);

    void meml1(const int nsrch, const float sdda[6][6], const float cdda[6][6],
               const float sa, const float ca, const float aima,
               const float rmaxa, const float sdefa,
               double sdd[6][6], double cdd[6][6], double gsd0[],
               double gcd0[], double &s, double &c, double &aim,
               double &sdef, double &rates, double &al0);

    void meml2(const int nsrch, const double sdd[6][6],
               const double cdd[6][6], const double gsd0[],
               const double gcd0[]);

    void meml22(const int nsrch, double sdd[6][6], double cdd[6][6],
                double gsd0[], double gcd0[]);

    void meml3(const int nsrch, const int k, const double sdd[6][6],
               const double cdd[6][6], int &ndim, double eval[6],
               double vecu[6][6], double &pen);

    void meml33(const int nd, const double amat[6][6], double val[6],
                double vec[6][6]);

    void meml4(const int nsrch, const int ndim, const double gsd0[6],
               const double gcd0[6], const double eval[6],
               const double vecu[6][6], const double c,
               const double aim, double gsd[6], double gcd[6],
               double &caim);

    void meml5(const int ndim, const double gsd[6],
               const double gcd[6], const double eval[6],
               const double c, const double caim,
               const double pen0, const double rates,
               const double al0, double &aaa, double &pp);

    void meml55(const int ndim, const double gsd[6],
                const double gcd[6], const double eval[6],
                const double c, const double caim, const double p,
                const double rates, double &a1, double &al,
                double &a2, double &aa);

    void meml6(const int nsrch, const int ndim, const double gsd[6],
               const double gcd[6], const double eval[6],
               const double vecu[6][6], const double s,
               const double c, const double pdist, const double alpha,
               double gqd[6], double xu[6], double w[6],
               double& snew, double& cnew, double& cpen,
               double& d2s);

    void meml66(const int nsrch, const double sdd[6][6], double w[6]);

    void meml7(const int nsrch, const double w[6], const double snew,
               const double cnew, const double d2s, const float sdefa,
               float wa[6], float &snewa, float &cnewa, float &rnewa);

    void meml8(const int nsrch, const int ndim, const double sdd[6][6],
               const double eval[6], const double gsd[6],
               const double gcd[6], const double gqd[6],
               const double xu[6], const float w[6],
               const float rnewa, const float s, const double c,
               const double caim, const float snew,
               const double cnew, const double cpen,
               const double alpha, const double pdist,
               const double pen0);

    void memop(const int k, const int l);

    void memp(float sdd[6][6]);

    void memp1(const int mj, const float f[], const float a2[],
               const float a3[], const float a4[],
               const float b[], const int m3, const float wtfile[],
               float a5[], float &s15, float &s25, float &s35,
               float &s45, float &s55);

    void memprm(const int method, const int level, const float aim,
                const float rmax, const float def, const float acc,
                float &c, float &test, float &cnew, float &s,
                float &rnew, float &snew, float &sumf);

    void memq(const float acc, float cdd[6][6]);

    void memq1(const int mk, const float a1[], const float a2[],
               const float a3[], const float a4[], const float a5[],
               const int m2, const float e[], const float acc,
               float &c15, float &c25, float &c35, float &c45, float &c55);

    void memqq(const float acc, const float w[5]);

    void memqq1(const int mk, const float a1[], const float a2[],
                const float a3[], const float a4[], const float a5[],
                const float w1, const float w2, const float w3,
                const float w4, const float w5, const int m2,
                const float e[], const float acc, float b[]);

    void memqq2(const int mk, const float a4[], const int m2,
                const float e[], const float acc, float b[]);

    void memr(const float w[5], float sdd[6][6]);

    void memr1(const int mj, const float a1[], const float a2[],
               const float a3[], const float a4[], const float a5[],
               const int m3, const float wtfile[], const float w1,
               const float w2, const float w3, const float w4, const float w5,
               float a6[], float &s14, float &s24, float &s34, float &s44);

    void memtr(const int k, const int l);

    //extern void opus(const int k, const int l);

    //extern void tropus(const int k, const int l);

    void udiag();

    void uinit();

    void uread(const int i);

    void ureset(const int i);

    void uwrite(const int i);

    void memcore(size_t mxbuf, size_t nmod, size_t ndat);

    size_t memsize(size_t nmod, size_t ndat);

    // Function definitions

    void memprm(const int method, const int level, const float aim,
                const float rmax, const float def, const float acc,
                float &c, float &test, float &cnew, float &s,
                float &rnew, float &snew, float &sumf){
        //                                                                  *
        //                  --------------------------------                *
        //                  "memsys"  maximum entropy system                *
        //                  --------------------------------                *
        //                                                                  *
        //               John Skilling                                      *
        //               dept. applied maths. and theoretical physics       *
        //               silver street, cambridge, england                  *
        //                                                                  *
        //         (c)   copyright 1985                                     *
        //               maximum entropy data consultants ltd               *
        //               33 north end, meldreth, royston, england           *
        // ******************************************************************
        //
        // purpose:
        //   perform one iterate of maximum entropy and update map
        //
        //   enter with <0> = map f, <20> = data d
        //
        //   and optionally:
        //
        //        <19> = default model ,  <18> = weights,
        //        <21> = accuracies = (2/variances)
        //
        //   exit with <0> = new map.
        //
        // ******************************************************************
        //
        //   entropy s := -1 - total( <0> * log <0>/e<19> ) / total<19>
        //
        //
        //   constraint c :=  sum (opus<0> - <20>)**2  * <21> / 2
        //
        // ******************************************************************
        // parameters:
        //   argument type  i/o  dimension   description
        //
        //    method  i     i      -         flags for procedural method
        //    level   i     i      -         flags for level of diagnostics
        //    aim     r     i      -         ultimate target of constraint
        //    rmax    r     i      -         dimensionless distance limit
        //    def     r     i      -         default level (if constant)
        //    acc     r     i      -         accuracy      (if constant)
        //
        //    c       r       o    -         constraint value on entry
        //    test    r       o    -         gradient direction misfit on entry
        //    cnew    r       o    -         predicted constraint value on exit
        //    s       r       o    -         entropy value on entry
        //    rnew    r       o    -         dimensionless distance moved
        //    snew    r       o    -         predicted entropy value on exit
        //    sumf    r       o    -         map total on exit
        //
        //
        // globals:
        //   variable  common block  i/o  description
        //
        //    Gbl::nj     mecomp     i    number of buffers in map-space area
        //    Gbl::mj     mecomp     i    map-space buffer length in words
        //    Gbl::nk     mecomp     i    number of buffers in data-space area
        //    Gbl::mk     mecomp     i    data-space buffer length in words
        //
        //    Gbl::ka     mecomp     i    allocation pointers for areas
        //    Gbl::kb     mecomp     i    base addresses in core for buffers
        //    Gbl::kc     mecomp     (o) core addresses for current buffers
        //    Gbl::kd     mecomp     (o) disc or upper core addresses of buffers
        //
        //    Gbl::iin    mecomp     i    fortran stream for user input (if any)
        //
        //    Gbl::l0     mecomp       o  flag for progress diagnostics
        //    Gbl::l1     mecomp       o  flag for numerical diagnostics
        //    Gbl::m0     mecomp       o  flag for number of search directions
        //    Gbl::m10    mecomp       o  flag for type of entropy
        //    Gbl::m11    mecomp       o  flag for constant total(map)
        //    Gbl::m20    mecomp       o  flag for constant or variable errors
        //    Gbl::m21    mecomp       o  flag for treatment of outliers
        //    Gbl::m3     mecomp       o  flag for weighting of map cells
        //    Gbl::m4     mecomp       o  flag for constant total(map)
        //
        //    Gbl::st     mecoms     i o  storage vector
        //
        //    Gbl::pr     mecomc      (o) read/write diagnostic flags
        //
        //
        // external calls:
        //            map-data
        //    memop     f d     call user transform routine opus
        //    memtr     f d     call user transpose routine tropus
        //
        //    mema        d     constraint value and gradient
        //    memb      f       entropy and contravariant gradients
        //    memc        d     scalars for first three search directions
        //    memcc       d     start calculating fourth search direction
        //    memd      f       fourth search direction
        //    meme        d     scalars for fourth search direction
        //    memee       d     start calculating fifth search direction
        //    memp      f       fifth search direction
        //    memq        d     scalars for fifth search direction
        //    memqq       d     revised fourth search direction transform
        //    memr      f       revised fourth search direction
        //
        //    meml       -      "control" - map increment coefficients
        //
        //    memf      f       new map
        //    memff     f       new map
        //
        //
        // key:
        //
        //    method = 10000*m4 + 1000*m3 + 100*m2 + 10*m1 + m0
        //
        //        m0 = 0 is full recursive generation of search directions
        //        m0 = 3 is three search directions (plus the current map)
        //        m0 = 4 is four  search directions (plus the current map)
        //
        //        m1 = 0 is constant default value def
        //        m1 = 1 is default given by exp(total(p log f))
        //        m1 = 3 is variable default levels in <20>
        //        m1 = 5 is "m1=0" with total(f) held fixed
        //        m1 = 6 is "m1=1" with total(f) held fixed
        //        m1 = 8 is "m1=3" with total(f) held fixed
        //
        //        m2 = 0 is chisquared with variable accuracies in <22>
        //        m2 = 1 is chisquared with constant accuracy acc
        //        m2 = 5 is "m2=0" but with outliers suppressed
        //        m2 = 6 is "m2=1" but with outliers suppressed
        //
        //        m3 = 0 is all cells weighted equally
        //        m3 = 1 is variable weights in <19>
        //
        //   level = 10*Gbl::l1 + Gbl::l0
        //
        //        Gbl::l0 = 0 gives no progress diagnostics
        //        Gbl::l0 = 1 gives names of main routines
        //        Gbl::l0 = 2 gives additionally read/write flowchart
        //        Gbl::l0 = 3 additionally calls udiag after each main routine
        //
        //        Gbl::l1 = 0 gives no numerical diagnostics
        //        Gbl::l1 = 1 gives s,test,c and snew,rnew,cnew
        //        Gbl::l1 = 2 gives additionally snew,rnew,cnew in each minor
        //                    loop
        //        Gbl::l1 = 3 gives additionally technical diagnostics from
        //                    "control"
        //        Gbl::l1 = 4 gives additionally the subspace scalars
        //
        //
        // flowchart:            map              data         scalars
        //   memex       r...................  ..w.......
        //   mema        ....................  rrr....w..   c
        //   tropus      .....w..............  .......r..
        //   memb        rww..r............rr  ..........   sdd11,12,13,22,23,33
        //                                                  s, sdef, sumf, test
        //   opus        .r..................  ...w......
        //   opus        ..r.................  ....w.....
        //   memc        ....................  .rrrr.....   cdd11,12,13,22,23,33
        //   control                                        w1,2,3
        //   memcc       ....................  .rrrr..w..
        //   tropus      .....w..............  .......r..
        //   memd        rrrw.r............r.  ..........   sdd14,24,34,44
        //   opus        ...r................  .....w....
        //   meme        ....................  .rrrrr....   cdd14,24,34,44
        //   control                                        w1,2,3,4, cnew
        // exit to memf ?
        //   memee       ....................  .rrrrr.w..
        //   tropus      .....w..............  .......r..
        //   memp        rrrrwr............r.  ..........   sdd15,25,35,45,55
        //   opus        ....r...............  ......w...
        //   memq        ....................  .rrrrrr...   cdd15,25,35,45,55
        // repeat
        //     control                                      w1,2,3,4,5, cnew
        //   exit to memff?
        //     memqq     ....................  ..rrrrrw..
        //               ....................  .....w.r..
        //               ....................  .r...r.w..
        //     memr      rrrrrw............r.  ..........   sdd14,24,34,44
        //               ...w.r..............  ..........
        //     tropus    .....w..............  .......r..
        //     memp      rrrrwr............r.  ..........   sdd15,25,35,45,55
        //     opus      ....r...............  ......w...
        //     memq      ....................  .rrrrrr...   cdd15,25,35,45,55
        // until false                                      cdd14,24,34,44
        //
        //   memf        rrrr.w............r.  ..........
        //               w....r..............  ..........   sumf
        // return
        //
        //   memff       rrrrrw............r.  ..........
        //               w....r..............  ..........   sumf
        // return
        //
        // history:
        //   john skilling    8 nov 1985     initial release
        //                    5 dec 1985     correction to meml22
        //                    9 dec 1985     overflow risk reduced in meml3
        //                   10 dec 1985     "entry" removed for ibm3081 linker
        //                   10 dec 1985     udiag input different for ibm3081
        //                   10 dec 1985     meml2 format altered for ibm3081
        //                   10 dec 1985     meml8 format altered for ibm3081
        //
        // notes:
        //    (1)  mem is the simpler call, but memprm is provided for users who
        //         wish to access the main diagnostic variables.
        //    (2)  in the definition of s, "total" means "sum weighted by <18>"
        //    (3)  areas <5> and <27> are workspace
        //
        //-----------------------------------------------------------------------

        float sdd[6][6], cdd[6][6], w[6];

        meinit();

        // method switches

        Gbl::m0  = method;
        int m1   = Gbl::m0/10;
        int m2   = m1/10;
        Gbl::m3  = m2/10;
        int m4   = Gbl::m3/10;
        Gbl::m0 -= 10*m1;
        m1      -= 10*m2;
        m2      -= 10*Gbl::m3;
        Gbl::m3 -= 10*m4;

        if((Gbl::m0 != 0 && Gbl::m0!=3 && Gbl::m0!=4) ||
           (m1!=0 && m1!=1 && m1!=3 &&
            m1!=5 && m1!=6 && m1!=8) ||
           (m2!=0 && m2!=1 && m2!=5 && m2!=6) ||
           (Gbl::m3!=0 && Gbl::m3!=1) ||
           (m4!=0 && m4!=1) ){
            std::cerr << "method = " << method << "not implemented" << std::endl;
            exit(EXIT_FAILURE);
        }

        Gbl::m11=m1/5;
        Gbl::m10=m1-Gbl::m11*5;
        Gbl::m21=m2/5;
        Gbl::m20=m2-Gbl::m21*5;

        // diagnostic switches

        Gbl::l0=level;
        Gbl::l1=Gbl::l0/10;
        Gbl::l0=Gbl::l0-10*Gbl::l1;
        if(   (Gbl::l0!=0 && Gbl::l0!=1 && Gbl::l0!=2 && Gbl::l0!=3) ||
              (Gbl::l1!=0 && Gbl::l1!=1 && Gbl::l1!=2 && Gbl::l1!=3
               && Gbl::l1!=4) ){
            std::cerr << "level = " << level << "   not implemented." << std::endl;
            exit(EXIT_FAILURE);
        }

        if(Gbl::l0>=1)
            std::cerr << " mem december 1985  method = " << m4 <<
                Gbl::m3 << m2 << m1 << Gbl::m0 << std::endl;

        // maximum entropy

        // Next line carries out first model --> data opus comp.
        // Stores result in <22>
        memex();

        // Computes c = reduced chi**2, and residual vector <27> =
        // <21>*(<22>-<20>) (i.e "accuracies" times (predicted data - data))
        mema(acc, c);
        if(Gbl::l1>=1) std::cerr << "       " <<
                           "                                            C    === "
                                 << c << std::endl;

        // Applies tropus to <27> --> <5>. <5> thus contains the chi**2
        // gradient one-form
        memtr(27,5);

        float sdef;
        // Computes the entropy and chi**2 gradient vectors which are stored
        // in <1> and <2>. These are the basic search directions for the
        // algorithm called f(grad(S)) and f(grad(C)) in Skilling & Bryan
        // they are called fgs and fgc in memb1
        memb(def, sdd, s, sdef, test);

        sumf=sdd[0][0];
        if(Gbl::l1>=1)
            std::cerr << "      S    === " << s << "    TEST === " << test << std::endl;
        memop(1,23);
        memop(2,24);

        // memc computes quadratic part of chi**2 sub-space
        memc(acc, cdd);
        meml(3,sdd,cdd,s,c,aim,rmax,sdef,snew,cnew,rnew,w);
        float speed0=rnew/2.;
        if(Gbl::l1>=2)
            std::cerr << "      SNEW === " << snew << "    DIST === " << rnew
                      << "    CNEW === " << cnew << std::endl;
        memcc(acc,w);
        memtr(27,5);
        memd(sdd);
        memop(3,25);
        meme(acc, cdd);
        meml(4,sdd,cdd,s,c,aim,rmax,sdef,snew,cnew,rnew,w);
        float speed=rnew/3.;
        if(Gbl::l1>=2)
            std::cerr << "      SNEW === " << snew << "    DIST === " << rnew
                      << "    CNEW === " << cnew << std::endl;
        if(speed<speed0 || Gbl::m0==3){
            if(Gbl::l1==1)
                std::cerr << "      SNEW === " << snew << "    DIST === " << rnew
                          << "    CNEW === " << cnew << std::endl;
            memf(w,sumf);
            return;
        }
        speed0=speed;
        memee(acc,w);
        memtr(27,5);
        memp(sdd);
        memop(4,26);
        memq(acc, cdd);
        int ndir=3;

        for(;;){
            meml(5,sdd,cdd,s,c,aim,rmax,sdef,snew,cnew,rnew,w);
            ndir++;
            speed=rnew/float(ndir);
            if(Gbl::l1>=2)
                std::cerr << "      SNEW === " << snew << "    DIST === " << rnew
                          << "    CNEW === " << cnew << std::endl;
            if(speed<speed0 || Gbl::m0==4){
                if(Gbl::l1==1)
                    std::cerr << "      SNEW === " << snew << "    DIST === " << rnew
                              << "    CNEW === " << cnew << std::endl;
                memff(w,sumf);
                return;
            }
            speed0=speed;
            memqq(acc,w);
            memr(w,sdd);
            memtr(27,5);
            memp(sdd);
            memop(4,26);
            memq(acc,cdd);
            cdd[0][3]=sdd[0][4];
            cdd[1][3]=sdd[1][4];
            cdd[2][3]=sdd[2][4];
            cdd[3][3]=sdd[3][4];
        }
    }

    void mema(const float acc , float &c){

        // purpose:
        //   constraint value and covariant gradient (data space)
        //
        //                  c  := ( <22> - <20> )**2 * <21> /2
        //
        //                <27> := ( <22> - <20> ) * <21>
        //
        // parameters:
        //   argument type  i/o  dimension   description
        //    acc     r    (i)     -         accuracy (alternative to <21>)
        //    c       r       o    -         constraint value
        //
        // globals:
        //   variable  common block  i/o  description
        //    st          mecoms     i o  storage vector
        //    nk,mk       mecomp     i    sizes (data space)
        //    ka,kb,kc,kd mecomp     i    initialised pointers
        //    l0          mecomp     i    diagnostics
        //    m20         mecomp     i    accuracy switch (0 =<21>, 1 =acc)
        //    m21         mecomp     i    outlier switch (0=normal, 1=suppressed)
        //
        // external calls:
        //    mema0       optional vector arithmetic for rms residual
        //    mema1       vector arithmetic
        //
        // history:
        //   john skilling    8 nov 1985     initial release
        //
        // notes:
        //    (1) constraint is chisquared, possibly corrected for outliers
        //

        if(Gbl::l0>=1) std::cerr << " mema" << std::endl;

        // set x = twice (6 * rmsresidual )**2  (or infinity if not needed)

        float x;
        if(Gbl::m21==0){
            x=1.e35;
        }else{
            float data = float(Gbl::nk*Gbl::mk);
            x=0.;
            for(int i=0;i<Gbl::nk;i++){
                uread(20);
                uread(22);
                if(Gbl::m20==0) uread(21);
                mema0(Gbl::mk,Gbl::st+Gbl::kc[20],Gbl::st+Gbl::kc[22],
                      Gbl::m20,Gbl::st+Gbl::kc[21],acc,x);
            }
            uinit();
            x *= 36./(36.+data);
        }

        c=0.;
        for(int i=0;i<Gbl::nk;i++){
            uread(20);
            uread(22);
            if(Gbl::m20==0) uread(21);
            mema1(Gbl::mk, Gbl::st+Gbl::kc[20], Gbl::st+Gbl::kc[22],
                  Gbl::m20, Gbl::st+Gbl::kc[21], acc, x,
                  Gbl::st+Gbl::kc[27], c);
            uwrite(27);
        }
        c *= 0.5;
        uinit();
        return;
    }

    void mema0(const int mk, const float d[], const float f[], const int m2,
               const float e[], const float acc, float &x){
        float a, r;
        float z1=0.;
        for(int i=0; i<mk; i++){
            if(m2==0){
                a=e[i];
            }else{
                a=acc;
            }
            if(a > 0.){
                r = f[i]-d[i];
                z1 += a*r*r;
            }
        }
        x += z1;
    }

    void mema1(int mk, const float d[], const float f[], const int m2, const float e[],
               const float acc, const float x, float gc[], float &c){
        // comes back with c = 2*(reduced chi-squared), gc=2*(f-d)/sigma**2/ndata
        float z1=0., a, r, y;

        for(int i=0; i<mk; i++){
            if(m2==0){
                a=e[i];
            }else{
                a=acc;
            }
            if(a > 0.){
                r = f[i]-d[i];
                y = r*(gc[i] = a*r);
                if(y > x){
                    gc[i] *= sqrt(x/y);
                    y=x;
                }
                z1 += y;
            }
        }
        c += z1;
    }


    void memb(const float def, float sdd[6][6], float &s, float &sdef,
              float &test){

        // purpose:
        //   contravariant second and third search directions
        //   ( f grad s  and   f grad c )  and initial entropy scalars
        //
        //   entropy  s  :=  - 1 - total( <0> * log <0>/e<19> ) / total<19>
        //
        //         sdef  :=  total<19>
        //
        //           <1> := ( log <19> - log <0> ) * metric
        //
        //           <2> := <5> * metric
        //
        //         sdd   := <i> * <j> / metric for ij = 00,01,02,11,12,22
        //
        //         total<i> = sum( <i> * <18> ), metric = <0> / <18>
        //
        //
        // parameters:
        //   argument type  i/o  dimension   description
        //    def     r    (i)     -         default level (if needed)
        //    sdd     r       o   6,6        entropy curvature scalars
        //    s       r       o    -         entropy value (dimensionless)
        //    sdef    r       o    -         total model
        //    test    r       o    -         1 - cos(angle(grad s, grad c))
        //
        // globals:
        //   variable  common block  i/o  description
        //    st          mecoms     i o  storage vector
        //    nj,mj       mecomp     i    sizes (map space)
        //    ka,kb,kc,kd mecomp     i    initialised pointers
        //    l0          mecomp     i    diagnostics
        //    m10         mecomp     i    model switch (0=def, 1=free, 3=<20>)
        //    m11         mecomp     i    total(map) switch (0 =off, 1 =on)
        //    m3          mecomp     i    weight switch (0 for none, 1 =<19>)
        //
        // external calls:
        //    memb0       optional vector arithmetic to determine default level
        //    memb1       vector arithmetic
        //
        // history:
        //   john skilling    8 nov 1985     initial release
        //
        // notes:
        //    (1)  s is calculated dimensionless, but gradients and curvatures
        //         are dimensional (multiplied by sdef). "control" corrects this.
        //    (2)  the factor 0.99999 guards against accidentally making test
        //         imaginary when projecting onto fixed total(map)
        //

        if(Gbl::l0 >= 1) std::cerr << " memb " << std::endl;

        // set defa = default value

        float defa;
        if(Gbl::m10!=1){
            defa=def;
        }else{
            float sumf=0., sumflf=0.;
            for(int i=0; i<Gbl::nj; i++){
                uread(0);
                if(Gbl::m3==1) uread(18);
                memb0(Gbl::mj,Gbl::st+Gbl::kc[0],Gbl::m3,Gbl::st+Gbl::kc[18],
                      sumf,sumflf);
            }
            uinit();
            defa=exp(sumflf/sumf);
        }

        sdef=0.0;
        sdd[0][0]=0.;
        sdd[0][1]=0.;
        sdd[0][2]=0.;
        sdd[1][1]=0.;
        sdd[1][2]=0.;
        sdd[2][2]=0.;

        for(int i=0; i<Gbl::nj; i++){
            uread(0);
            uread(5);
            if(Gbl::m10==3) uread(19);
            if(Gbl::m3 ==1) uread(18);
            memb1(Gbl::mj, Gbl::st+Gbl::kc[0], Gbl::st+Gbl::kc[5],
                  Gbl::m3, Gbl::st+Gbl::kc[18], Gbl::m10,
                  Gbl::st+Gbl::kc[19], defa, Gbl::st+Gbl::kc[1],
                  Gbl::st+Gbl::kc[2], sdef,
                  sdd[0][0],sdd[0][1],sdd[0][2],sdd[1][1],sdd[1][2],sdd[2][2]);
            uwrite(1);
            uwrite(2);
        }
        uinit();

        s=(sdd[0][0]+sdd[0][1])/sdef-1.0;
        float ss, sc, cc;
        if(Gbl::m11==0){
            ss=sdd[1][1];
            sc=sdd[1][2];
            cc=sdd[2][2];
        }else{
            ss=sdd[1][1]-0.99999*sdd[0][1]*sdd[0][1]/sdd[0][0];
            sc=sdd[1][2]-0.99999*sdd[0][1]*sdd[0][2]/sdd[0][0];
            cc=sdd[2][2]-0.99999*sdd[0][2]*sdd[0][2]/sdd[0][0];
        }

        test=1.-sc/(sqrt(ss+1.e-35)*sqrt(cc+1.e-35));
        return;
    }
 
    void memb0(const int mj, const float f[], const int m3, 
               const float wtfile[], float &sumf, float &sumflf){

        float z1=0., z2=0., weight;

        for(int i=0; i<mj; i++){
            if(m3==1){
                weight=wtfile[i];
            }else{
                weight=1.;
            }
            z1 += weight*f[i];
            if(f[i]>0.) z2 += weight*f[i]*log(f[i]);
        }
        sumf   += z1;
        sumflf += z2;
        return;
    }
 
    void memb1(const int mj, const float f[], const float gc[], const int m3,
               const float wtfile[], const int m1, float const amodel[], 
               const float def, float fgs[], float fgc[],
               float &sdef, float &s11, float &s12,
               float &s13, float &s22, float &s23,
               float &s33){

        float z0=0., z1=0., z2=0., z3=0., z4=0., z5=0., z6=0.;
        float weight, deflt, gs, amet, x, y;

        for(int i=0; i<mj; i++){
            if(m3==1){
                weight=wtfile[i];
            }else{
                weight=1.;
            }
            if(m1>=3){
                deflt=amodel[i];
            }else{
                deflt=def;
            }
            if(f[i]>0. && deflt>0.){
                gs=weight*log(deflt/f[i]);
            }else{
                gs=0.;
            }
            amet = f[i]/weight;
            // x,y are the contra-variant components of the entropy and chi**2
            // gradients
            x    = amet*gs;
            y    = amet*gc[i];
            z0  += weight*deflt;
            z1  += weight*f[i];
            z2  += weight*x;
            z3  += weight*y;
            z4  += x*gs;
            z5  += y*gs;
            z6  += y*gc[i];
            fgs[i] = x;
            fgc[i] = y;
            if(std::abs(gs) > 1.e25 || std::abs(gc[i]) > 1.e25 || std::abs(y) > 1.e25){
                std::cerr << "weight,x,y,gs,gc,i,deflt,f="
                          << weight << ", "
                          << x << ", "
                          << y << ", "
                          << gs << ", "
                          << gc[i] << ", "
                          << i << ", "
                          << deflt << ", "
                          << f[i] << std::endl;
            }
        }

        sdef += z0;
        s11  += z1;
        s12  += z2;
        s13  += z3;
        s22  += z4;
        s23  += z5;
        s33  += z6;
        return;
    }

    void memc (const float acc, float cdd[6][6]){

        // purpose:
        //   constraint scalars of first three search directions
        //
        //   cdd   := <22+i> * <21> * (22+j>      for     ij = 00,01,02,11,12,22
        //
        // parameters:
        //   argument type  i/o  dimension   description
        //    acc     r    (i)     -         accuracy (alternative to <21>)
        //    cdd     r       o   6,6        constraint curvature scalars
        //
        // globals:
        //   variable  common block  i/o  description
        //    st          mecoms     i    storage vector
        //    nk,mk       mecomp     i    sizes (data space)
        //    ka,kb,kc,kd mecomp     i    initialised pointers
        //    l0          mecomp     i    diagnostics
        //    m20         mecomp     i    accuracy switch (0 =<21>, 1 =acc)
        //    m21         mecomp     i    outlier switch (0=normal, 1=suppressed)
        //
        // external calls:
        //    memc1       vector arithmetic
        //
        // history:
        //   john skilling    8 nov 1985     initial release
        //
        // notes:
        //

        if(Gbl::l0>=1) std::cerr << "  mem" << std::endl;
        cdd[0][0]=0.;
        cdd[0][1]=0.;
        cdd[0][2]=0.;
        cdd[1][1]=0.;
        cdd[1][2]=0.;
        cdd[2][2]=0.;
        for(int i=0; i<Gbl::nk; i++){
            uread(22);
            uread(23);
            uread(24);
            if(Gbl::m20==0) uread(21);
            memc1( Gbl::mk, Gbl::st+Gbl::kc[22], Gbl::st+Gbl::kc[23], 
                   Gbl::st+Gbl::kc[24], Gbl::m20, Gbl::st+Gbl::kc[21], acc,
                   cdd[0][0], cdd[0][1], cdd[0][2], cdd[1][1], 
                   cdd[1][2], cdd[2][2]);
        }
        uinit();
        return;
    }

 
    void memc1(const int mk, const float a1[], const float a2[], 
               const float a3[], const int m2, const float e[], 
               const float acc, float &c11, float &c12, float &c13,
               float &c22, float &c23, float &c33){

        float z1=0., z2=0., z3=0., z4=0., z5=0., z6=0.;
        float a;

        for(int i=0; i<mk; i++){
            if(m2==0){
                a=e[i];
            }else{
                a=acc;
            }
            if(a > 0.){
                z1 += a1[i]*a1[i]*a;
                z2 += a1[i]*a2[i]*a;
                z3 += a1[i]*a3[i]*a;
                z4 += a2[i]*a2[i]*a;
                z5 += a2[i]*a3[i]*a;
                z6 += a3[i]*a3[i]*a;
            }
        }
        c11 += z1;
        c12 += z2;
        c13 += z3;
        c22 += z4;
        c23 += z5;
        c33 += z6;
        return;
    }
 
 
    void memcc(const float acc, const float w[3]){

        // purpose:
        //   start generating fourth search direction
        //
        //           <27> := ( w1*<22> + w2*<23> + w3*<24> ) * <21>
        //
        //
        // parameters:
        //   argument type  i/o  dimension   description
        //    acc     r    (i)     -         accuracy (alternative to <21>)
        //    w       r     i      3         map increment coefficients
        //
        // globals:
        //   variable  common block  i/o  description
        //    st          mecoms     i o  storage vector
        //    nk,mk       mecomp     i    sizes (data space)
        //    ka,kb,kc,kd mecomp     i    initialised pointers
        //    l0          mecomp     i    diagnostics
        //    m20         mecomp     i    accuracy switch (0 =<21>, 1 =acc)
        //    m21         mecomp     i    outlier switch (0=normal, 1=suppressed)
        //
        // external calls:
        //    memcc1       vector arithmetic
        //
        // history:
        //   john skilling    8 nov 1985     initial release
        //
        // notes:
        //
  
        if(Gbl::l0>=1) std::cerr << " memcc " << std::endl;
        for(int i=0; i<Gbl::nk; i++){
            uread(22);
            uread(23);
            uread(24);
            if(Gbl::m20==0) uread(21);
            memcc1(Gbl::mk, Gbl::st+Gbl::kc[22], Gbl::st+Gbl::kc[23], 
                   Gbl::st+Gbl::kc[24], Gbl::m20, Gbl::st+Gbl::kc[21],
                   acc, w[0], w[1], w[2], Gbl::st+Gbl::kc[27]);
            uwrite(27);
        }
        uinit();
        return;
    }

    void memcc1(const int mk, const float a1[], const float a2[], 
                const float a3[], const int m2, const float e[], 
                const float acc, const float w1, const float w2, 
                const float w3, float b[]){

        float a;
        for(int i=0; i<mk; i++){
            if(m2==0){
                a=e[i];
            }else{
                a=acc;
            }
            if(a > 0.) b[i] = a*(w1*a1[i]+w2*a2[i]+w3*a3[i]);
        }
        return;
    }
 
    void memd(float sdd[6][6]){

        // purpose:
        //   contravariant fourth search direction and its entropy scalars
        //
        //          <3> := <5> * metric    ,    metric = <0> / <18>
        //
        //        sdd   := <5> * <i>      for     i = 0,1,2,3
        //           i4
        //
        // parameters:
        //   argument type  i/o  dimension   description
        //    sdd     r       o   6,6        entropy curvature scalars
        //
        // globals:
        //   variable  common block  i/o  description
        //    st          mecoms     i o  storage vector
        //    nj,mj       mecomp     i    sizes (map space)
        //    ka,kb,kc,kd mecomp     i    initialised pointers
        //    l0          mecomp     i    diagnostics
        //    m3          mecomp     i    weight switch (0 for none, 1 =<19>)
        //
        // external calls:
        //    memd1       vector arithmetic
        //
        // history:
        //   john skilling    8 nov 1985     initial release
        //
        // notes:
        //

        if(Gbl::l0>=1) std::cerr << " memd" << std::endl;
        sdd[0][3]=0.;
        sdd[1][3]=0.;
        sdd[2][3]=0.;
        sdd[3][3]=0.;
        for(int i=0; i<Gbl::nj; i++){
            uread(0);
            uread(1);
            uread(2);
            uread(5);
            if(Gbl::m3==1) uread(18);
            memd1(Gbl::mj, Gbl::st+Gbl::kc[0], Gbl::st+Gbl::kc[1], 
                  Gbl::st+Gbl::kc[2], Gbl::st+Gbl::kc[5], Gbl::m3,
                  Gbl::st+Gbl::kc[18], Gbl::st+Gbl::kc[3],
                  sdd[0][3], sdd[1][3], sdd[2][3], sdd[3][3]);
            uwrite(3);
        }
        uinit();
        return;
    }
 
    void memd1(const int mj, const float f[], const float a2[], 
               const float a3[], const float b[],
               const int m3, const float wtfile[], float a4[], 
               float &s14,float &s24,float &s34, float &s44){

        float z1=0.,z2=0.,z3=0.,z4=0.,weight,amet,x;

        for(int i=0; i<mj; i++){
            if(m3==1){
                weight=wtfile[i];
            }else{
                weight=1.;
            }
            amet=f[i]/weight;
            x=amet*b[i];
            z1 += f[i]*b[i];
            z2 += a2[i]*b[i];
            z3 += a3[i]*b[i];
            z4 += x*b[i];
            a4[i] = x;
        }
        s14 += z1;
        s24 += z2;
        s34 += z3;
        s44 += z4;
        return;
    }
 
    void meme (const float acc, float cdd[6][6]){

        // purpose:
        //   constraint scalars of fourth search direction
        //
        //   cdd   := <22+i> * <21> * (22+j>      for     ij = 03,13,23,33
        //      ij
        //
        // parameters:
        //   argument type  i/o  dimension   description
        //    acc     r    (i)     -         accuracy (alternative to <21>)
        //    cdd     r       o   6,6        constraint curvature scalars
        //
        // globals:
        //   variable  common block  i/o  description
        //    st          mecoms     i    storage vector
        //    nk,mk       mecomp     i    sizes (data space)
        //    ka,kb,kc,kd mecomp     i    initialised pointers
        //    l0          mecomp     i    diagnostics
        //    m20         mecomp     i    accuracy switch (0 =<21>, 1 =acc)
        //    m21         mecomp     i    outlier switch (0=normal, 1=suppressed)
        //
        // external calls:
        //    meme1       vector arithmetic
        //
        // history:
        //   john skilling    8 nov 1985     initial release
        //
        // notes:
        //

        if(Gbl::l0>=1) std::cerr << " meme" << std::endl;
        cdd[0][3]=0.;
        cdd[1][3]=0.;
        cdd[2][3]=0.;
        cdd[3][3]=0.;
        for(int i=0; i<Gbl::nk; i++){
            uread(22);
            uread(23);
            uread(24);
            uread(25);
            if(Gbl::m20==0) uread(21);
            meme1(Gbl::mk, Gbl::st+Gbl::kc[22], Gbl::st+Gbl::kc[23], 
                  Gbl::st+Gbl::kc[24], Gbl::st+Gbl::kc[25],
                  Gbl::m20, Gbl::st+Gbl::kc[21], acc, cdd[0][3],
                  cdd[1][3], cdd[2][3], cdd[3][3]);
        }
        uinit();
        return;
    }
 
    void meme1(const int mk, const float a1[], const float a2[], 
               const float a3[], const float a4[], const int m2, 
               const float e[], const float acc, float &c14,
               float &c24, float &c34, float &c44){

        float z1=0.,z2=0.,z3=0.,z4=0., x, a;

        for(int i=0; i<mk; i++){
            if(m2==0){
                a=e[i];
            }else{
                a=acc;
            }
            if(a > 0.){
                x   = a4[i]*a;
                z1 += a1[i]*x;
                z2 += a2[i]*x;
                z3 += a3[i]*x;
                z4 += a4[i]*x;
            }
        }
        c14 += z1;
        c24 += z2;
        c34 += z3;
        c44 += z4;
        return;
    }
 
    void memee(const float acc, const float w[4]){

        // purpose:
        //   start generating fifth search direction
        //
        //       <27> := ( w1*<22> + w2*<23> + w3*<24> + w4*<25> ) * <21>
        //
        //
        // parameters:
        //   argument type  i/o  dimension   description
        //    acc     r    (i)     -         accuracy (alternative to <21>)
        //    w       r     i      4         map increment coefficients
        //
        // globals:
        //   variable  common block  i/o  description
        //    st          mecoms     i o  storage vector
        //    nk,mk       mecomp     i    sizes (data space)
        //    ka,kb,kc,kd mecomp     i    initialised pointers
        //    l0          mecomp     i    diagnostics
        //    m20         mecomp     i    accuracy switch (0 =<21>, 1 =acc)
        //    m21         mecomp     i    outlier switch (0=normal, 1=suppressed)
        //
        // external calls:
        //    memee1       vector arithmetic
        //
        // history:
        //   john skilling    8 nov 1985     initial release
        //
        // notes:
        //

        if(Gbl::l0>=1) std::cerr << " memee" << std::endl;
        for(int i=0; i<Gbl::nk; i++){
            uread(22);
            uread(23);
            uread(24);
            uread(25);
            if(Gbl::m20==0) uread(21);
            memee1(Gbl::mk, Gbl::st+Gbl::kc[22], Gbl::st+Gbl::kc[23],
                   Gbl::st+Gbl::kc[24], Gbl::st+Gbl::kc[25],
                   Gbl::m20, Gbl::st+Gbl::kc[21], acc,
                   Gbl::st+Gbl::kc[27], w[0], w[1], w[2], w[3]);
            uwrite(27);
        }
        uinit();
        return;
    }

    void memee1(const int mk, const float a1[], const float a2[], 
                const float a3[], const float a4[], const int m2, 
                const float e[], const float acc, float b[],
                const float w1, const float w2, const float w3, 
                const float w4){

        float a;
        for(int i=0; i<mk; i++){
            if(m2==0){
                a=e[i];
            }else{
                a=acc;
            }
            if(a > 0.) b[i] = a*(w1*a1[i]+w2*a2[i]+w3*a3[i]+w4*a4[i]);
        }
        return;
    }
 
    void memp(float sdd[6][6]){

        // purpose:
        //   contravariant fifth search direction and its entropy scalars
        //
        //          <4> := <5> * metric    ,    metric = <0> / <18>
        //
        //        sdd   := <5> * <i>      for     i = 0,1,2,3,4
        //           i5
        //
        // parameters:
        //   argument type  i/o  dimension   description
        //    sdd     r       o   6,6        entropy curvature scalars
        //
        // globals:
        //   variable  common block  i/o  description
        //    st          mecoms     i o  storage vector
        //    nj,mj       mecomp     i    sizes (map space)
        //    ka,kb,kc,kd mecomp     i    initialised pointers
        //    l0          mecomp     i    diagnostics
        //    m3          mecomp     i    weight switch (0 for none, 1 =<18>)
        //
        // external calls:
        //    memp1       vector arithmetic
        //
        // history:
        //   john skilling    8 nov 1985     initial release
        //
        // notes:
        //

        if(Gbl::l0>=1) std::cerr << " memp" << std::endl;
        sdd[0][4]=0.;
        sdd[1][4]=0.;
        sdd[2][4]=0.;
        sdd[3][4]=0.;
        sdd[4][4]=0.;
        for(int i=0; i<Gbl::nj; i++){
            uread(0);
            uread(1);
            uread(2);
            uread(3);
            uread(5);
            if(Gbl::m3==1) uread(18);
            memp1(Gbl::mj, Gbl::st+Gbl::kc[0], Gbl::st+Gbl::kc[1], 
                  Gbl::st+Gbl::kc[2], Gbl::st+Gbl::kc[3], Gbl::st+Gbl::kc[5],
                  Gbl::m3, Gbl::st+Gbl::kc[18], Gbl::st+Gbl::kc[4], sdd[0][4],
                  sdd[1][4], sdd[2][4], sdd[3][4], sdd[4][4]);
            uwrite(4);
        }
        uinit();
        return;
    }
 
    void memp1(const int mj, const float f[], const float a2[], 
               const float a3[], const float a4[],
               const float b[], const int m3, const float wtfile[], 
               float a5[], float &s15, float &s25, float &s35, 
               float &s45, float &s55){

        float z1=0., z2=0., z3=0., z4=0., z5=0., weight, amet, x;

        for(int i=0; i<mj; i++){
            if(m3==1){
                weight=wtfile[i];
            }else{
                weight=1.;
            }
            amet = f[i]/weight;
            x    = amet*b[i];
            z1  += f[i]*b[i];
            z2  += a2[i]*b[i];
            z3  += a3[i]*b[i];
            z4  += a4[i]*b[i];
            z5  += x*b[i];
            a5[i]=x;
        }
        s15 += z1;
        s25 += z2;
        s35 += z3;
        s45 += z4;
        s55 += z5;
        return;
    }
 
    void memq(const float acc , float cdd[6][6]){

        // purpose:
        //   constraint scalars of fifth search direction
        //
        //   cdd   := <22+i> * <21> * (22+j>      for     ij = 04,14,24,34,44
        //      ij
        //
        // parameters:
        //   argument type  i/o  dimension   description
        //    acc     r    (i)     -         accuracy (alternative to <21>)
        //    cdd     r       o   6,6        constraint curvature scalars
        //
        // globals:
        //   variable  common block  i/o  description
        //    st          mecoms     i    storage vector
        //    nk,mk       mecomp     i    sizes (data space)
        //    ka,kb,kc,kd mecomp     i    initialised pointers
        //    l0          mecomp     i    diagnostics
        //    m20         mecomp     i    accuracy switch (0 =<21>, 1 =acc)
        //    m21         mecomp     i    outlier switch (0=normal, 1=suppressed)
        //
        // external calls:
        //    memq1       vector arithmetic
        //
        // history:
        //   john skilling    8 nov 1985     initial release
        //
        // notes:
        //
        //-----------------------------------------------------------------------

        if(Gbl::l0>=1) std::cerr << " memq" << std::endl;
        cdd[0][4]=0.;
        cdd[1][4]=0.;
        cdd[2][4]=0.;
        cdd[3][4]=0.;
        cdd[4][4]=0.;
        for(int i=0; i<Gbl::nk; i++){
            uread(22);
            uread(23);
            uread(24);
            uread(25);
            uread(26);
            if(Gbl::m20==0) uread(21);
            memq1(Gbl::mk, Gbl::st+Gbl::kc[22], Gbl::st+Gbl::kc[23],
                  Gbl::st+Gbl::kc[24], Gbl::st+Gbl::kc[25], Gbl::st+Gbl::kc[26],
                  Gbl::m20, Gbl::st+Gbl::kc[21], acc, cdd[0][4], cdd[1][4],
                  cdd[2][4], cdd[3][4], cdd[4][4]);
        }
        uinit();
        return;
    }
 
    void memq1(const int mk, const float a1[], const float a2[], 
               const float a3[], const float a4[], const float a5[], 
               const int m2, const float e[], const float acc, 
               float &c15, float &c25, float &c35, float &c45, float &c55){

        float a, z1=0., z2=0., z3=0., z4=0., z5=0., x;

        for(int i=0; i<mk; i++){
            if(m2==0){
                a=e[i];
            }else{
                a=acc;
            }
            if(a > 0.){
                x   = a5[i]*a;
                z1 += a1[i]*x;
                z2 += a2[i]*x;
                z3 += a3[i]*x;
                z4 += a4[i]*x;
                z5 += a5[i]*x;
            }
        }
        c15 += z1;
        c25 += z2;
        c35 += z3;
        c45 += z4;
        c55 += z5;
    }
 
    void memqq(const float acc, const float w[5]){

        // purpose:
        //   set data-space version of revised fourth search direction
        //   and start generating revised fifth search direction
        //
        //       <25> :=  w1*<22> + w2*<23> + w3*<24> + w4*<25> + w5*<26>
        //
        //       <27> :=  <25> * <21>
        //
        // parameters:
        //   argument type  i/o  dimension   description
        //    acc     r    (i)     -         accuracy (alternative to <21>)
        //    w       r     i      5         map increment coefficients
        //
        // globals:
        //   variable  common block  i/o  description
        //    st          mecoms     i o  storage vector
        //    nk,mk       mecomp     i    sizes (data space)
        //    ka,kb,kc,kd mecomp     i    initialised pointers
        //    l0          mecomp     i    diagnostics
        //    m20         mecomp     i    accuracy switch (0 =<21>, 1 =acc)
        //    m21         mecomp     i    outlier switch (0=normal, 1=suppressed)
        //
        // external calls:
        //    memqq1       vector arithmetic for new <25>
        //    memk         copy vector
        //    memqq2       vector arithmetic for new <27>
        //
        // history:
        //   john skilling    8 nov 1985     initial release
        //
        // notes:
        //

        if(Gbl::l0>=1) std::cerr << " memqq" << std::endl;
        for(int i=0; i<Gbl::nk; i++){
            uread(22);
            uread(23);
            uread(24);
            uread(25);
            uread(26);
            memqq1(Gbl::mk, Gbl::st+Gbl::kc[22], Gbl::st+Gbl::kc[23], 
                   Gbl::st+Gbl::kc[24], Gbl::st+Gbl::kc[25], Gbl::st+Gbl::kc[26], 
                   w[0], w[1], w[2], w[3], w[4], Gbl::m20, Gbl::st+Gbl::kc[21], 
                   acc, Gbl::st+Gbl::kc[27]);
            uwrite(27);
        }
        uinit();
        memk(27,25);
        for(int i=0; i<Gbl::nk; i++){
            uread(25);
            if(Gbl::m20==0) uread(21);
            memqq2(Gbl::mk, Gbl::st+Gbl::kc[25], Gbl::m20, Gbl::st+Gbl::kc[21],
                   acc, Gbl::st+Gbl::kc[27]);
            uwrite(27);
        }
        uinit();
    }
 
    // modification here to allow negative errors to be ignored (TRM)

    void memqq1(const int mk, const float a1[], const float a2[], 
                const float a3[], const float a4[], const float a5[], 
                const float w1, const float w2, const float w3,
                const float w4, const float w5, const int m2, 
                const float e[], const float acc, float b[]){
        float a;
        for(int i=0; i<mk; i++){
            if(m2==0){
                a=e[i];
            }else{
                a=acc;
            }
            if(a > 0.)
                b[i] =w1*a1[i]+w2*a2[i]+w3*a3[i]+w4*a4[i]+w5*a5[i];
        }
    }
 
    void memqq2(const int mk, const float a4[], const int m2, 
                const float e[], const float acc, float b[]){

        float a;

        for(int i=0; i<mk; i++){
            if(m2==0){ 
                a=e[i];
            }else{
                a=acc;
            }
            if(a > 0.) b[i] =a*a4[i];
        }
    }
 
    void memr(const float w[5], float sdd[6][6]){

        // purpose:
        //   revised fourth search direction and its entropy scalars
        //
        //          <3> := w1*<0> + w2*<1> + w3*<2> + w4*<3> +w5*<4>
        //
        //        sdd   := <3> * <i> / metric     for     i = 0,1,2,3
        //           i4
        //                             metric = <0> / <18>
        //
        // parameters:
        //   argument type  i/o  dimension   description
        //    w       r     i      5         map increment coefficients
        //    sdd     r       o   6,6        entropy curvature scalars
        //
        // globals:
        //   variable  common block  i/o  description
        //    st          mecoms     i o  storage vector
        //    nj,mj       mecomp     i    sizes (map space)
        //    ka,kb,kc,kd mecomp     i    initialised pointers
        //    l0          mecomp     i    diagnostics
        //    m3          mecomp     i    weight switch (0 for none, 1 =<18>)
        //
        // external calls:
        //    memr1       vector arithmetic
        //    memj        copy vector
        //
        // history:
        //   john skilling    8 nov 1985     initial release
        //
        // notes:
        //

        if(Gbl::l0>=1) std::cerr << " memr" << std::endl;
        sdd[0][3]=0.;
        sdd[1][3]=0.;
        sdd[2][3]=0.;
        sdd[3][3]=0.;
        for(int i=0; i<Gbl::nj; i++){
            uread(0);
            uread(1);
            uread(2);
            uread(3);
            uread(4);
            if(Gbl::m3==1) uread(18);
            memr1(Gbl::mj, Gbl::st+Gbl::kc[0], Gbl::st+Gbl::kc[1], 
                  Gbl::st+Gbl::kc[2], Gbl::st+Gbl::kc[3], Gbl::st+Gbl::kc[4],
                  Gbl::m3, Gbl::st+Gbl::kc[18], w[0], w[1], w[2], w[3], w[4],
                  Gbl::st+Gbl::kc[5], sdd[0][3], sdd[1][3], sdd[2][3], sdd[3][3]);
            uwrite(5);
        }
        uinit();
        memj(5,3);
    }
 
    void memr1(const int mj, const float a1[], const float a2[], 
               const float a3[], const float a4[], const float a5[], 
               const int m3, const float wtfile[], const float w1,
               const float w2, const float w3, const float w4, const float w5, 
               float a6[], float &s14, float &s24, float &s34, float &s44){

        float z1=0., z2=0., z3=0., z4=0., weight, x, zmet;

        for(int i=0; i<mj; i++){
            x=w1*a1[i]+w2*a2[i]+w3*a3[i]+w4*a4[i]+w5*a5[i];
            if(m3==1){
                weight=wtfile[i];
            }else{
                weight=1.;
            }
            if(a1[i]>0.){
                zmet=weight/a1[i];
                z1 += x*weight;
                z2 += x*a2[i]*zmet;
                z3 += x*a3[i]*zmet;
                z4 += x*x*zmet;
            }
            a6[i]=x;
        }
        s14 += z1;
        s24 += z2;
        s34 += z3;
        s44 += z4;
    }
 
    void memf(const float w[4], float &sumf){

        // purpose:
        //   new map from four search directions (including old map)
        //
        //             <0>  :=  w1*<0> + w2*<1> + w3*<2> + w4*<3> ,
        //                      possibly rescaled to constant total(map)
        //
        //             sumf :=  total(map)
        //
        // parameters:
        //   argument type  i/o  dimension   description
        //    w       r     i      4         map increment coefficients
        //    sumf    r       o    -         map total
        //
        // globals:
        //   variable  common block  i/o  description
        //    st          mecoms     i    storage vector
        //    nj,mj       mecomp     i    sizes (map space)
        //    ka,kb,kc,kd mecomp     i    initialised pointers
        //    l0          mecomp     i    diagnostics
        //    m11         mecomp     i    total(map) switch (0 =off, 1 =on)
        //    m3          mecomp     i    weight switch (0 for none, 1 =<18>)
        //
        // external calls:
        //    memf1       vector arithmetic
        //    memf2       copy (and rescale)
        //
        // history:
        //   john skilling    8 nov 1985     initial release
        //
        // notes:
        //    (1) map components f may not decrease by more than 80 percent,
        //        or below 1.e-12
        //

        float sum=0.;
  
        if(Gbl::l0>=1) std::cerr << " memf" << std::endl;
        for(int i=0; i<Gbl::nj; i++){
            uread(0);
            uread(1);
            uread(2);
            uread(3);
            if(Gbl::m3==1) uread(18);
            memf1(Gbl::mj, Gbl::st+Gbl::kc[0], Gbl::st+Gbl::kc[1], 
                  Gbl::st+Gbl::kc[2], Gbl::st+Gbl::kc[3], Gbl::m3,
                  Gbl::st+Gbl::kc[18], w[0], w[1], w[2], w[3],
                  Gbl::st+Gbl::kc[5],sum);
            uwrite(5);
        }
        uinit();
        float scale;
        if(Gbl::m11==0){
            sumf=sum;
            scale=1.;
        }else{
            scale=sumf/sum;
        }
        for(int i=0; i<Gbl::nj; i++){
            uread(5);
            memf2(Gbl::mj, scale, Gbl::st+Gbl::kc[5], Gbl::st+Gbl::kc[0]);
            uwrite(0);
        }
        uinit();
    }
 
    void memf1(const int mj, const float a1[], const float a2[], 
               const float a3[], const float a4[], const int m3, 
               const float wtfile[], const float w1,
               const float w2, const float w3, const float w4, 
               float f[], float &sum){

        float z1=0., weight, x, d;
        for(int i=0; i<mj; i++){
            if(m3==1){
                weight=wtfile[i];
            }else{
                weight=1.;
            }
            x=a1[i];
            d=w1*a1[i]+w2*a2[i]+w3*a3[i]+w4*a4[i];
            if(x>0.) x=std::max(std::max(0.2f*x,x+d),1.e-12f);
            z1 += x*weight;
            f[i]=x;
        }
        sum += z1;
    }
 
    void memf2(const int mj, const float scale, const float old[], float f[]){
        for(int i=0; i<mj; i++){
            f[i]=old[i]*scale;
        }
    }

    void memff(const float w[5], float &sumf){

        // purpose:
        //   new map from five search directions (including old map)
        //
        //         <0>  :=  w1*<0> + w2*<1> + w3*<2> + w4*<3> + w5*<4> ,
        //                  possibly rescaled to constant total(map)
        //
        //         sumf :=  total(map)
        //
        // parameters:
        //   argument type  i/o  dimension   description
        //    w       r     i      5         map increment coefficients
        //    sumf    r       o    -         map total
        //
        // globals:
        //   variable  common block  i/o  description
        //    st          mecoms     i o  storage vector
        //    nj,mj       mecomp     i    sizes (map space)
        //    ka,kb,kc,kd mecomp     i    initialised pointers
        //    l0          mecomp     i    diagnostics
        //    m11         mecomp     i    total(map) switch (0 =off, 1 =on)
        //    m3          mecomp     i    weight switch (0 for none, 1 =<18>)
        //
        // external calls:
        //    memff1      vector arithmetic
        //    memff2      copy (and rescale)
        //
        // history:
        //   john skilling    8 nov 1985     initial release
        //
        // notes:
        //    (1) map components f may not decrease by more than 80 percent,
        //        or below 1.e-12
        //

        float sum=0.;
        if(Gbl::l0>=1) std::cerr << " memff" << std::endl;
        for(int i=0; i<Gbl::nj; i++){
            uread(0);
            uread(1);
            uread(2);
            uread(3);
            uread(4);
            if(Gbl::m3==1) uread(18);
            memff1(Gbl::mj, Gbl::st+Gbl::kc[0], Gbl::st+Gbl::kc[1],
                   Gbl::st+Gbl::kc[2], Gbl::st+Gbl::kc[3], Gbl::st+Gbl::kc[4],
                   Gbl::m3, Gbl::st+Gbl::kc[18], w[0], w[1], w[2], w[3], w[4],
                   Gbl::st+Gbl::kc[5], sum);
            uwrite(5);
        }
        uinit();
        float scale;
        if(Gbl::m11==0){
            sumf=sum;
            scale=1.;
        }else{
            scale=sumf/sum;
        }
        for(int i=0; i<Gbl::nj; i++){
            uread(5);
            memff2(Gbl::mj, scale, Gbl::st+Gbl::kc[5], Gbl::st+Gbl::kc[0]);
            uwrite(0);
        }
        uinit();
    }
 
    void memff1(const int mj, const float a1[], const float a2[], 
                const float a3[], const float a4[], const float a5[], 
                const int m3, const float wtfile[], const float w1,
                const float w2, const float w3, const float w4, 
                const float w5, float f[], float &sum){

        float z1=0., weight, x, d;
        for(int i=0; i<mj; i++){
            if(m3==1){
                weight=wtfile[i];
            }else{
                weight=1.;
            }
            x=a1[i];
            d=w1*a1[i]+w2*a2[i]+w3*a3[i]+w4*a4[i]+w5*a5[i];
            if(x>0.) x=std::max(std::max(0.2f*x,x+d),1.e-12f);
            z1 += x*weight;
            f[i]=x;
        }
        sum += z1;
    }

    void memff2(const int mj, const float scale, const float old[], float f[]){
        for(int i=0;i<mj;i++){
            f[i]=old[i]*scale;
        }
    }
 
    void memj(const int m, const int n){

        // purpose:
        //   copy map-space file.          <n> := <m>
        //
        // parameters:
        //   argument type  i/o  dimension   description
        //    m       i     i      -         input  area number
        //    n       i     i      -         output area number
        //
        // globals:
        //   variable  common block  i/o  description
        //    st          mecoms     i o  storage vector
        //    nj,mj       mecomp     i    sizes (map space)
        //    ka,kb,kc,kd mecomp     i    initialised pointers
        //
        // external calls:
        //    memj1       vector arithmetic
        //
        // history:
        //   john skilling    8 nov 1985     initial release
        //
        // notes:
        //    (1)  no diagnostics produced
        //
  
        for(int i=0; i<Gbl::nj; i++){
            uread(m);
            memj1(Gbl::mj, Gbl::st+Gbl::kc[m], Gbl::st+Gbl::kc[n]);
            uwrite(n);
        }
        uinit();
    }
 
    void memj1(const int m, const float a[], float b[]){
        for(int i=0; i<m; i++){
            b[i]=a[i];
        }
    }
 
    void memk(const int m, const int n){

        //
        // purpose:
        //   copy data-space file.          <n> := <m>
        //
        // parameters:
        //   argument type  i/o  dimension   description
        //    m       i     i      -         input  area number
        //    n       i     i      -         output area number
        //
        // globals:
        //   variable  common block  i/o  description
        //    st          mecoms     i o  storage vector
        //    nk,mk       mecomp     i    sizes (data space)
        //    ka,kb,kc,kd mecomp     i    initialised pointers
        //
        // external calls:
        //    memk1       vector arithmetic
        //
        // history:
        //   john skilling    8 nov 1985     initial release
        //
        // notes:
        //    (1)  no diagnostics produced
        //

        for(int i=0; i<Gbl::nk; i++){
            uread(m);
            memk1(Gbl::mk, Gbl::st+Gbl::kc[m], Gbl::st+Gbl::kc[n]);
            uwrite(n);
        }
        uinit();
    }

    void memk1(const int m, const float a[], float b[]){
        for(int i=0; i<m; i++){
            b[i]=a[i];
        }
    }
 
    void memex(){

        //
        // purpose:
        //   model the actual behaviour of the experiment.
        //
        //                 <22> :=  actual opus transform <0>
        //
        // linear experiments are treated by a call to the
        //   differential response routine opus.
        //
        // nonlinear experiments should transform an image on area 0
        // to mock data on file 22, and may possibly need to adjust
        // area 20, area 21, and the required value of c.
        //
        //
        // parameters:
        //   argument type  i/o  dimension   description
        //     -
        //
        // globals:
        //   variable  common block  i/o  description
        //     -
        //
        // external calls:
        //    memop        driver routine for opus
        //
        // history:
        //   john skilling    8 nov 1985     initial release
        //
        // notes:
        //    (1) this implementation is for linear experiments.
        //

        memop(0,22);

    }

    void memop(const int k, const int l){

        // purpose:
        //   call user's transform routine opus and re-initialise pointers
        //
        //                 <l> :=  opus <k>
        //
        // parameters:
        //   argument type  i/o  dimension   description
        //    k       i     i      -         input  area number  (map space)
        //    l       i     i      -         output area number (data space)
        //
        // globals:
        //   variable  common block  i/o  description
        //    l0          mecomp     i    diagnostic stream and flag
        //    pr          mecomc       o  read/write flags
        //
        // external calls:
        //    opus        user routine
        //
        // history:
        //   john skilling    8 nov 1985     initial release
        //
        // notes:
        //

        Gbl::pr[k] = 'r';
        Gbl::pr[l] = 'w';
        
        // Extern function call to user defined routine
        // Usually defined by doppler.cc during the python build
        opus(k,l);

        // re-initialise to protect the user
        uinit();
    }

    void memtr(const int k, const int l){

        // purpose:
        //   call user's transform routine tropus and re-initialise pointers
        //
        //                 <l> :=  tropus <k>
        //
        // parameters:
        //   argument type  i/o  dimension   description
        //    k       i     i      -         input  area number (data space)
        //    l       i     i      -         output area number  (map space)
        //
        // globals:
        //   variable  common block  i/o  description
        //    l0          mecomp     i    diagnostic stream and flag
        //    pr          mecomc       o  read/write flags
        //
        // external calls:
        //    tropus        user routine
        //
        // history:
        //   john skilling    8 nov 1985     initial release
        //
        // notes:

        Gbl::pr[k] = 'r';
        Gbl::pr[l] = 'w';

        // Extern function call to user defined routine
        // Usually defined by doppler.cc during the python build
        tropus(k,l);

        // re-initialise to protect the user
        uinit();
    }

    void meml(const int nsrch, const float sdda[6][6],
              const float cdda[6][6], const float sa,
              const float ca, const float aima,
              const float rmaxa, const float sdefa,
              float &snewa, float &cnewa, float &rnewa,
              float wa[6]){

        // purpose:
        //                          " control "
        //
        //      maximise quadratic model of entropy s (=objective function)
        //      over the target value (or an approach to it) of the quadratic
        //      model of the constraint c within a given distance limit.
        //      enter with models
        //       s(w) = s + w.gsd - w.sdd.w/2  ,  c(w) = c + w.gcd + w.cdd.w/2
        //      exit with map increments w and diagnostics.
        //
        //
        // parameters:
        //   argument type  i/o  dimension   description
        //    nsrch   i     i      -         number of search directions
        //    sdda    r     i   nsrch,nsrch  entropy metric
        //    cdda    r     i   nsrch,nsrch  constraint matrix
        //    sa      r     i      -         current value of entropy
        //    ca      r     i      -         current value of constraint
        //    aima    r     i      -         ultimate target of constraint
        //    rmaxa   r     i      -         dimensionless distance limit
        //    sdefa   r     i      -         scale factor for entropy
        //    snewa   r       o    -         predicted value of entropy
        //    cnewa   r       o    -         predicted value of constraint
        //    rnewa   r       o    -         dimensionless distance moved
        //    wa      r       o    6         map increment components
        //
        // globals:
        //   variable  common block  i/o  description
        //    m11       mecomp       i    total(map) switch (0 =off, 1 =on)
        //    l1        mecomp       i    level of numerical diagnostics
        //
        // external calls:
        //    meml1     copy to internal (double precision) variables
        //    meml2       optional diagnostics of input model "scalars"
        //    meml22      optional orthogonalisation if total(map) held fixed
        //    meml3     simultaneous diagonalisation of entropy and constraint
        //    meml4     transform to eigenvector coords
        //    meml5     * central control * find alpha and distance penalty *
        //    meml6     transform back to original coordinates
        //    meml66      optional coefficient of "map" if total(map) held fixed
        //    meml7     copy to external (single precision) parameters
        //    meml8       optional numerical diagnostics
        //
        // history:
        //   john skilling    8 nov 1985     initial release
        //
        // notes:
        //   (1) there is some protection against a non-convex
        //       constraint function, but this version of control
        //       is not designed to handle such problems robustly
        //   (2) input parameters are preserved
        //   (3) only the upper triangles j>i of sdda and cdda are read
        //   (4) the map total sumf and the gradients of s and c are extracted
        //       from sdda, which involves the following assumptions:
        //      a) scalar products in map-space have been calculated using
        //         the entropy curvature -gradgrads
        //      b) the first  search direction is f=map
        //      c) the second search direction is contravariant grads
        //      d) the third  search direction is contravariant gradc
        //   (5) if total(map) is held fixed, the central control routines use
        //       only directions 2,3,... , orthogonalised to direction 1
        //   (6) the distance limit is   w.sdd.w /sumf  <  rmax**2
        //   (7) all matrices and vectors are stored in arrays dimensioned 6
        //
        // notation:
        //       s = entropy
        //       c = constraint (e.g. chisquared)
        //       prefix g = gradient, as in gs and gc
        //       suffix d = downstairs covariant index, as in gsd, gcd, sdd, cdd
        //       suffix u = upstairs contravariant index, as in xu
        //       suffix a = single-precision parameter for external use only
        //

        double sdd[6][6], cdd[6][6], gsd0[6], gcd0[6], s, c, aim;
        double sdef, rates, pen0, al0, eval[6], vecu[6][6], gsd[6];
        double gcd[6], caim, aaa, pp, gqd[6], xu[6], w[6], snew, cnew;
        double cpen, d2s;

        if(Gbl::l0>=1) std::cerr << " control" << std::endl;
        if( nsrch<1 || nsrch>6 ){
            std::cerr << "  invalid nsrch =" << nsrch << std::endl;
            exit(EXIT_FAILURE);
        }

        // k = first relevant search direction

        int k = Gbl::m11;

        // copy to internal (double precision) variables

        meml1(nsrch,sdda,cdda,sa,ca,aima,rmaxa,sdefa,
              sdd,cdd,gsd0,gcd0,s,c,aim,sdef,rates,al0);

        // optional diagnostics of input model "scalars"

        if(Gbl::l1>=4) meml2(nsrch,sdd,cdd,gsd0,gcd0);

        // orthogonalise?

        if(k==1) meml22(nsrch,sdd,cdd,gsd0,gcd0);

        // diagonalise

        int ndim;
        meml3(nsrch-k,k,sdd,cdd,ndim,eval,vecu,pen0);

        // transform to eigenvector coords

        meml4(nsrch-k,ndim,gsd0+k,gcd0+k,eval,vecu,c,aim,gsd,gcd,caim);

        // central control * find alpha and distance penalty 

        meml5(ndim,gsd,gcd,eval,c,caim,pen0,rates,al0, aaa,pp);

        // transform back to original coordinates

        meml6(nsrch-k,ndim,gsd,gcd,eval,vecu,s,c,pp,aaa,
              gqd,xu,w+k,snew,cnew,cpen,d2s);

        // orthogonalised?

        if(k==1) meml66(nsrch,sdd,w);

        // copy to external (single precision) parameters

        meml7(nsrch,w,snew,cnew,d2s,sdefa,wa,snewa,cnewa,rnewa);

        // optional numerical diagnostics

        if(Gbl::l1>=3) meml8(nsrch,ndim,sdd,eval,gsd,gcd,gqd,xu,wa,
                             rnewa,sa,c,caim,snewa,cnew,cpen,aaa,pp,pen0);
    }

    void meml1(const int nsrch, const float sdda[6][6], 
               const float cdda[6][6], const float sa, 
               const float ca, const float aima,
               const float rmaxa, const float sdefa,
               double sdd[6][6], double cdd[6][6], double gsd0[],
               double gcd0[], double &s, double &c, double &aim,
               double &sdef, double &rates, double &al0){

        // purpose:
        //   copy to internal (double precision) variables
        //
        // parameters:
        //   argument type  i/o  dimension   description
        //    nsrch   i     i      -         number of original search dirns
        //    sdda    r     i   nsrch,nsrch  entropy metric
        //    cdda    r     i   nsrch,nsrch  constraint matrix
        //    sa      r     i      -         current value of entropy
        //    ca      r     i      -         current value of constraint
        //    aima    r     i      -         ultimate target of constraint
        //    rmaxa   r     i      -         dimensionless distance limit
        //    sdefa   r     i      -         scale factor for entropy
        //    sdd     r*8     o nsrch,nsrch  entropy metric
        //    cdd     r*8     o nsrch,nsrch  constraint matrix
        //    gsd0    r*8     o   nsrch      components of gsd
        //    gcd0    r*8     o   nsrch      components of gcd
        //    s       r*8     o    -         current value of entropy
        //    c       r*8     o    -         current value of constraint
        //    aim     r*8     o    -         ultimate target of constraint
        //    sdef    r*8     o    -         scale factor for entropy
        //    rates   r*8     o    -         maximum distance-squared (*sdef)
        //    al0     r*8     o    -         suggested first choice of alpha
        //
        // globals:
        //   variable  common block  i/o  description
        //    -
        //
        // external calls:
        //    -
        //
        // history:
        //   john skilling    8 nov 1985     initial release
        //
        // notes:
        //    (1) internal calculations use s scaled by sdefa
        //

        s     = double(sa*sdefa);
        c     = double(ca);
        aim   = double(aima);
        sdef  = double(sdefa);
        rates = sdef*double(rmaxa*rmaxa);

        for(int i=0; i<nsrch; i++){
            for(int j=i; j<nsrch; j++){
                sdd[i][j] = double(sdda[i][j]);
                cdd[i][j] = double(cdda[i][j]);
                sdd[j][i] = double(sdda[i][j]);
                cdd[j][i] = double(cdda[i][j]);
            }
            gsd0[i] = sdd[i][1];
            gcd0[i] = sdd[i][2];
        }
        al0 =sqrt(sdd[2][2])/sqrt(sdd[1][1]+rates);
    }

    void meml2(const int nsrch, const double sdd[6][6],
               const double cdd[6][6], const double gsd0[],
               const double gcd0[]){

        // purpose:
        //   optional diagnostics of input "scalars"
        //
        // parameters:
        //   argument type  i/o  dimension   description
        //    nsrch   i     i      -         number of original search dirns
        //    sdd     r*8   i   nsrch,nsrch  entropy metric
        //    cdd     r*8   i   nsrch,nsrch  constraint matrix
        //    gsd0    r*8   i     nsrch      components of gsd
        //    gcd0    r*8   i     nsrch      components of gcd
        //
        // globals:
        //   variable  common block  i/o  description
        //    -
        //
        // external calls:
        //    -
        //
        // history:
        //   john skilling    8 nov 1985     initial release
        //                   10 dec 1985     format altered for ibm3081
        //
        // notes:
        //

        std::cerr << " gsd0";
        for(int j=0;j<nsrch;j++)
            std::cerr << " " << gsd0[j];
        std::cerr << std::endl;

        std::cerr << " gcd0";
        for(int j=0;j<nsrch;j++)
            std::cerr << " " << gcd0[j];
        std::cerr << std::endl;

        for(int i=0;i<nsrch;i++){
            std::cerr << " sdd";
            for(int j=0;j<nsrch;j++)
                std::cerr << " " << sdd[i][j];
            std::cerr << std::endl;
        }
        for(int i=0;i<nsrch;i++){
            std::cerr << " cdd";
            for(int j=0;j<nsrch;j++)
                std::cerr << " " << cdd[i][j];
            std::cerr << std::endl;
        }
    }

    void meml22(const int nsrch, double sdd[6][6],
                double cdd[6][6], double gsd0[],
                double gcd0[]){

        // purpose:
        //    orthogonalise search directions 2,3,...,nsrch to direction 1 (=map)
        //    if total(map) held fixed.
        //
        // parameters:
        //   argument type  i/o  dimension   description
        //    nsrch   i     i      -         number of original search dirns
        //    sdd     r*8   i o nsrch,nsrch  entropy metric
        //    cdd     r*8   i o nsrch,nsrch  constraint matrix
        //    gsd0    r*8   i o   nsrch      components of gsd
        //    gcd0    r*8   i o   nsrch      components of gcd
        //
        // globals:
        //   variable  common block  i/o  description
        //    -
        //
        // external calls:
        //    -
        //
        // history:
        //   john skilling    8 nov 1985     initial release
        //                    5 dec 1985     correction
        //
        // notes:
        //    (1) the factors of 0.99999 guard against accidentally making
        //        matrices non-positive
        //

        double x, y;

        for(int i=1; i<nsrch; i++){
            x=0.99999*sdd[0][i]/sdd[0][0];
            for(int j=1; j<nsrch; j++){
                y=0.99999*sdd[0][j]/sdd[0][0];
                sdd[i][j] -= y*sdd[0][i];
                cdd[i][j] -= x*cdd[0][j]+y*cdd[0][i]-x*y*cdd[0][0];
            }
            gsd0[i] -= x*gsd0[0];
            gcd0[i] -= x*gcd0[0];
        }
    }

    void meml3(const int nsrch, const int k, const double sdd[6][6],
               const double cdd[6][6], int &ndim, double eval[6],
               double vecu[6][6], double &pen){

        //
        // purpose:
        //   simultaneous diagonalisation of entropy and constraint matrices
        //                             l       l                     l
        //   solves  cdd(i,j) * vecu(j)  = eval  * sdd(i,j) * vecu(j)
        //   for eigenvalues eval and eigenvectors vecu, for l=1,2,...,ndim
        //
        // parameters:
        //   argument type  i/o  dimension   description
        //    nsrch   i     i      -         number of original search dirns
        //    k       i     i      -         offset inside sdd, cdd
        //    sdd     r*8   i   nsrch,nsrch  entropy metric    (symmetric)
        //    cdd     r*8   i   nsrch,nsrch  constraint matrix (symmetric)
        //    ndim    i     o    -         number usefully independent dirns
        //    eval    r*8   o   ndim       eigenvalues of cdd over sdd metric
        //    vecu    r*8   o nsrch,ndim   eigenvectors of cdd over sdd metric
        //    pen     r*8   o    -         uncertainty in lowest eval
        //
        // globals:
        //   variable  common block  i/o  description
        //    -
        //
        // external calls:
        //    meml33    diagonalise symmetric matrix
        //
        // history:
        //   john skilling    8 nov 1985     initial release
        //                    9 dec 1985     overflow risk reduced
        //
        // notes:
        //   (1) eigenvalues eval are returned in increasing order
        //   (2) eigenvectors vecu are normalised to  v.sdd.v=1
        //   (3) input matrices sdd and cdd are preserved
        //   (4) ndim usually returns as the full dimension nsrch of the metric
        //       sdd, but if sdd is ill-conditioned the routine will use only
        //       the subspace spanned by its larger eigenvalues, and ndim will
        //       return with an appropriately smaller value.
        //

        double amet[6][6],amat[6][6],w1[6][6],w2[6][6],svec[6][6];
        double sval[6], cval[6], x, y;

        // eps is related to machine accuracy

        const double err=3.0e-5;

        // copy normalised input matrices to amat=cdd, amet=sdd

        for(int i=0; i<nsrch; i++){
            for(int j=0; j<nsrch; j++){
                x=sqrt(sdd[i+k][i+k])*sqrt(sdd[j+k][j+k]);
                amat[i][j] = cdd[i+k][j+k]/(x+1.0e-18);
                amet[i][j] = sdd[i+k][j+k]/(x+1.0e-18);
            }
        }

        // entropy eigenvalues sval and eigenvectors svec
        // ndim = number of usefully independent directions

        meml33(nsrch,amet,sval,svec);

        x = err*sval[nsrch-1];
        int l;
        for(l=0; sval[l]<x; l++);
        ndim = nsrch - l;

        // rotate cdd to sdd eigenvector space

        int m;
        for(int i=0; i<nsrch; i++){
            for(int j=0; j<ndim; j++){
                m=j+l;
                x=0.0;
                for(int p=0; p<nsrch; p++)
                    x += amat[i][p]*svec[p][m];
                w2[i][j] = x;
            }
        }

        for(int i=0; i<ndim; i++){
            m = i+l;
            for(int j=0; j<ndim; j++){
                x=0.0;
                for(int p=0; p<nsrch; p++)
                    x += svec[p][m]*w2[p][j];
                w1[i][j] = x;
            }
        }

        // cdd is now w1, sdd is now diag(sval)
        // squeeze cdd=w1 to isotropise sdd

        for(int i=0; i<ndim; i++)
            for(int j=0; j<ndim; j++)
                w1[i][j] /= sqrt(sval[i+l])*sqrt(sval[j+l]);

        // cdd=w1 eigenvalues cval and eigenvectors w2
        meml33(ndim,w1,cval,w2);

        // complete squeeze of w2 back to sdd eigenvector space
        for(int i=0; i<ndim; i++)
            for(int j=0; j<ndim; j++)
                w2[i][j] /= sqrt(sval[i+l]);

        // rotate w2 to original space
        for(int j=0; j<ndim; j++){
            for(int i=0; i<nsrch; i++){
                x=0.0;
                for(int p=0; p<ndim; p++){
                    x += svec[i][p+l]*w2[p][j];
                }
                w1[i][j] = x/(sqrt(sdd[i+k][i+k])+1.0e-35);
            }
        }

        // polish eigenvalues eval and recover eigenvectors vecu

        for(int j=0; j<ndim; j++){
            x=0.0;
            y=0.0;
            for(int p=0; p<nsrch; p++){
                for(int q=0; q<nsrch; q++){
                    x += w1[p][j]*cdd[p+k][q+k]*w1[q][j];
                    y += w1[p][j]*sdd[p+k][q+k]*w1[q][j];
                }
            }
            eval[j] = x/y;

            for(int p=0; p<nsrch; p++)
                vecu[p][j] = w1[p][j]/sqrt(y);
        }

        // uncertainty in eval[0]  (care with overflow risk)

        x=0.0;
        for(int p=0; p<nsrch; p++){
            for(int q=0; q<nsrch; q++){
                y  = vecu[p][0]*cdd[p+k][q+k];
                x += fabs(y*vecu[q][0]);
                y  = vecu[p][0]*sdd[p+k][q+k];
                x += fabs(eval[0]*y*vecu[q][0]);
            }
        }
        pen=err*x/nsrch;
    }


    void meml33(const int nd, const double amat[6][6], double val[6],
                double vec[6][6]){

        //
        // purpose:
        //   diagonalisation of symmetric matrix
        //                                 k      k         k
        //   solves      amat(i,j) * vec(j)  = val  * vec(i)
        //   for eigenvalues val and eigenvectors vec, for k=1,2,...,nd
        //
        // parameters:
        //   argument type  i/o  dimension   description
        //    nd      i     i      -         number of dimensions
        //    amat    r*8   i     nd,nd      matrix (from 6*6 array)
        //    val     r*8     o    nd        eigenvalues
        //    vec     r*8     o   nd,nd      eigenvectors
        //
        // globals:
        //   variable  common block  i/o  description
        //    -
        //
        // external calls:
        //    -
        //
        // history:
        //   john skilling    8 nov 1985     initial release
        //
        //
        // notes:
        //   (1) eigenvalues val are returned in increasing order
        //   (2) eigenvectors vec are normalised to v.v=1
        //   (3) input matrix amat is preserved
        //   (4) only the upper triangle j>i of amat(i,j) is read
        //
        //   (5) algorithm is to repeatedly square amat until the largest
        //       eigenvalue dominates, then to subtract off that eigenvector,
        //       and repeat until all the eigenvectors have been found.
        //       this nd**4 algorithm would be too slow for large matrices,
        //       but is robust for small.
        //
        // NB TRM: I don't believe the docs are quite right here: this
        // routine starts by adding a constant times the identity
        // matrix to the input matrix to ensure that the result is
        // positive definite. It does so even if the input is already
        // positive definite. The answers returned are the eigenvalues
        // and vectors of this modified matrix, but with the
        // eigenvalues corrected back by the applied constant. I don't
        // understand the rationale behind the constant chosen
        // (04/01/2019), but it is probably one of those slightly
        // heuristic procedures associated with names like Powell,
        // Fletcher etc. It's a bit odd though.

        double c, z[6][6], p[6][6], q[6][6], y[6], x[6];

        // eps is related to machine accuracy
        const double eps=3.0e-16;

        // copy amat to a positive definite matrix z (by adding c to eigenvals)
        double a=0.0, b=0.0, t=0.0;
        for(int i=0; i<nd; i++){
            a += amat[i][i];
            for(int j=0; j<nd; j++){
                c = amat[i][j];
                if(fabs(c)<1.0e-18) c=0.0;
                b += c*c;
            }
            t += 1.0;
        }
        c = std::max(b-a*a/t,0.0);
        c = a/t-sqrt(c)-eps*sqrt(b);

        for(int i=0; i<nd; i++){
            for(int j=i; j<nd; j++){
                z[i][j] = amat[i][j];
                z[j][i] = amat[i][j];
            }
            z[i][i] -= c;
        }

        // lmax = maximum number of inner loop iterates

        t = -log(eps)/(eps*log(2.0));
        t = log(t)/log(2.0);
        int lmax = int(t+2.0);

        int n = nd-1;

        // outer loop over n for successively smaller eigenvalues

        for(;;){
            t=0.0;
            for(int i=0; i<nd; i++){
                for(int j=0; j<nd; j++)
                    p[i][j] = z[i][j];
                t += p[i][i];
            }

            int l = 0;

            // inner loop over l for squaring p and setting t=trace

            do{
                l++;

                t = 1.0/t;
                for(int i=0; i<nd; i++){
                    for(int j=0; j<nd; j++){
                        q[i][j] = t*p[i][j];
                        if(fabs(q[i][j]) < 1.0e-18) q[i][j] =0.0;
                    }
                }

                t=0.0;
                for(int i=0; i<nd; i++){
                    for(int j=i; j<nd; j++){
                        a=0.0;
                        for(int k=0; k<nd; k++)
                            a += q[i][k]*q[k][j];
                        if(fabs(a)<1.0e-18) a=0.0;
                        p[i][j] = a;
                        p[j][i] = a;
                    }
                    t += p[i][i];
                }
            }while( t < 1.0-eps*10.0 && l <= lmax );

            // end inner loop when p is dyadic
            //
            // k = largest column = estimate of current largest eigenvector

            int k=0;
            a=0.0;
            for(int i=0; i<nd; i++){
                if(p[i][i] > a){
                    a = p[i][i];
                    k = i;
                }
            }

            // p(p(largest column)) = better estimate x of eigenvector

            for(int i=0; i<nd; i++){
                a=0.0;
                for(int j=0; j<nd; j++)
                    a += p[i][j]*p[j][k];

                if(fabs(a) < 1.0e-18) a=0.0;
                y[i]=a;
            }

            t=0.0;
            for(int i=0; i<nd; i++){
                a=0.0;
                for(int j=0; j<nd; j++)
                    a += p[i][j]*y[j];
                if(fabs(a)<1.0e-18) a=0.0;
                t += a*a;
                x[i]=a;
            }

            // repeat..  orthogonalise x to previous eigenvectors

            k=nd-1;
            if(k != n){
                do{
                    a=0.0;
                    for(int i=0; i<nd; i++)
                        a += x[i]*vec[i][k];
                    for(int i=0; i<nd; i++){
                        b=x[i]-a;
                        if(fabs(b)<1.0e-18) b=0.0;
                        x[i]=b;
                    }
                    k--;
                }while(k>n);
            }

            // ..while
            // normalise eigenvector x

            a=0.0;
            for(int i=0; i<nd; i++)
                a += x[i]*x[i];

            a =1.0/sqrt(a);

            // copy eigenvector x into output array vec

            for(int i=0; i<nd; i++){
                x[i] *= a;
                vec[i][n] = x[i];
            }

            // set eigenvalue val directly from the eigenvector

            a=0.0;
            for(int i=0; i<nd; i++){
                for(int j=0; j<nd; j++){
                    b = x[i]*x[j];
                    if(fabs(b)<1.0e-18) b=0.0;
                    a += z[i][j]*b;
                }
            }
            // note the addition of the constant 'c' at this
            // point -- means the eigenvectors do not match
            // the eigen values
            val[n] = a+c;

            // finish ?   (if full set of eigenvectors has been found)
            if(--n < 0)  return;

            // otherwise, remove dyadic  x.xtranspose  from matrix z
            for(int i=0; i<nd; i++){
                for(int j=i; j<nd; j++){
                    b =x[i]*x[j];
                    if(fabs(b)<1.0e-18) b=0.0;
                    z[i][j] -= a*b;
                    if(fabs(z[i][j])<1.0e-18) z[i][j] =0.0;
                    z[j][i] = z[i][j];
                }
            }
        }
    }


    void meml4(const int nsrch, const int ndim, const double gsd0[6],
               const double gcd0[6], const double eval[6],
               const double vecu[6][6], const double c,
               const double aim, double gsd[6], double gcd[6],
               double &caim){

        //
        // purpose:
        //   transfer problem to eigenvector coordinates
        //
        // parameters:
        //   argument type  i/o  dimension   description
        //    nsrch   i     i      -         number of original search dirns
        //    ndim    i     i      -         number of independent eigenvectors
        //    gsd0    r*8   i     nsrch      components of gsd in original coords
        //    gcd0    r*8   i     nsrch      components of gcd in original coords
        //    eval    r*8   i     ndim       eigenvalues of c matrix
        //    vecu    r*8   i   nsrch,ndim   eigenvector matrix from (6,6) array
        //    c       r*8   i      -         current value of constraint
        //    aim     r*8   i      -         ultimate target of constraint
        //    gsd     r*8     o   ndim       components of gsd in diagonal coords
        //    gcd     r*8     o   ndim       components of gcd in diagonal coords
        //    caim    r*8     o    -         realistic target value of constraint
        //
        // globals:
        //   variable  common block  i/o  description
        //    -
        //
        // external calls:
        //    -
        //
        // history:
        //   john skilling    8 nov 1985     initial release
        //
        // notes:
        //   (1) caim is set not more than 4/5 of the way down to the lowest
        //       value of constraint accessible in the entire subspace.
        //

        // rewrite vectors gsd,gcd in eigenvector coordinates

        for(int j=0; j<ndim; j++){
            gcd[j]=0.0;
            gsd[j]=0.0;
            for(int k=0; k<nsrch; k++){
                gcd[j] += gcd0[k]*vecu[k][j];
                gsd[j] += gsd0[k]*vecu[k][j];
            }
        }

        // metric sdd is now unit, matrix cdd is now diagonal(eval)
        // set caim

        double cmin = c;
        for(int k=0; k<ndim; k++)
            cmin -= 0.5*gcd[k]*gcd[k]/std::max(eval[k],1.0e-12);

        caim = std::max((4.0*cmin+c)/5.0, aim);
    }

    void meml5(const int ndim, const double gsd[6],
               const double gcd[6], const double eval[6],
               const double c, const double caim,
               const double pen0, const double rates,
               const double al0, double &aaa, double &pp){

        //
        // purpose:
        //   find smallest acceptable distance penalty pp and its alpha aaa
        //   of the final subspace solution which maximises s over the target c
        //   subject to the distance constraint
        //
        // parameters:
        //   argument type  i/o  dimension   description
        //    ndim    i     i      -         number of dimensions in subspace
        //    gsd     r*8   i     ndim       components of gsd in diagonal coords
        //    gcd     r*8   i     ndim       components of gcd in diagonal coords
        //    eval    r*8   i     ndim       eigenvalues of c matrix (now diag.)
        //    c       r*8   i      -         current value of constraint
        //    caim    r*8   i      -         realistic target value of constraint
        //    pen0    r*8   i      -         penalty offset (usually small)
        //    rates   r*8   i      -         maximum limit of distance**2 moved
        //    al0     r*8   i      -         suggested starting value of alpha
        //    aaa     r*8     o    -         alpha of final subspace solution
        //    pp      r*8     o    -         penalty of final subspace solution
        //
        // globals:
        //   variable  common block  i/o  description
        //    -
        //
        // external calls:
        //    meml55    alpha-loop
        //
        // history:
        //   john skilling    8 nov 1985     initial release
        //
        // notes:
        //   (1) there is some protection against negative eigenvales eval of the
        //       constraint function, but this version of control is not
        //       designed to handle such non-convex problems robustly
        //   (2) distance penalty p increases the effective constraint value cpen
        //       and not the entropy
        //

        // eps related to machine accuracy

        const double eps=3.0e-16;

        double eps4 = sqrt(sqrt(eps));

        // ee is eigenvalue norm

        double ee=1.0e-35;

        for(int k=0; k<ndim; k++)
            ee += fabs(eval[k])/float(k+1);

        // initialise penalty loop

        double pen =0.0;
        double p1  =0.0;
        double p2  =1.0e30;
        pp  = 0.0;
        aaa = 0.0;
        int loop=0;

        // begin penalty loop

        double a1, a2, al, aa, p;
        do{
            if(loop==1)  pen=(2.0*p1+1.0+p1/p2)/(1.0+(2.0+p1)/p2);
            loop=1;

            // initialise alpha loop

            al=al0;
            a1=al*eps;
            a2=al/eps;
            aa=0.0;
            p=pen0+0.1*pen*ee;

            // call alpha-loop routine to find correct value of alpha

            meml55(ndim,gsd,gcd,eval,c,caim,p,rates,a1,al,a2,aa);

            // check acceptability of alpha returned from its loop

            if( aa==a2 ){

                // penalty large enough. store it and prepare to decrease it

                aaa=aa;
                pp=pen;
                p2=pen;
            }else{

                // penalty too small. prepare to increase it

                p1=pen;
            }
        }while( p2-p1 > eps4*pen );

        // exit with smallest acceptable distance penalty and its alpha

        pp=pen0+0.1*pp*ee;
    }


    void meml55(const int ndim, const double gsd[6],
                const double gcd[6], const double eval[6],
                const double c, const double caim, const double p,
                const double rates, double &a1, double &al,
                double &a2, double &aa){

        //
        // purpose:
        //      find the value of alpha which maximises s over the target value
        //      of c within the distance limit in the diagonalised subspace.
        //
        // parameters:
        //   argument type  i/o  dimension   description
        //    ndim    i     i      -         number of dimensions in subspace
        //    gsd     r*8   i     ndim       components of gsd in diagonal coords
        //    gcd     r*8   i     ndim       components of gcd in diagonal coords
        //    eval    r*8   i     ndim       eigenvalues of c matrix (now diag.)
        //    c       r*8   i      -         current value of constraint
        //    caim    r*8   i      -         target value (within distance limit)
        //    p       r*8   i      -         distance penalty
        //    rates   r*8   i      -         maximum limit of distance**2 moved
        //    a1      r*8   i o    -         lower limit on alpha in binary chop
        //    al      r*8   i o    -         alpha
        //    a2      r*8   i o    -         upper limit on alpha in binary chop
        //    aa      r*8   i o    -         smallest a2 within distance limit
        //
        // globals:
        //   variable  common block  i/o  description
        //    -
        //
        // external calls:
        //    -
        //
        // history:
        //   john skilling    8 nov 1985     initial release
        //
        // notes:
        //   (1) alpha is found to full machine accuracy
        //

        // set extreme values

        double cmin=std::min(c,caim);
        double cmax=std::max(c,caim);

        // revise eigenvalues by distance penalty

        double e[6];
        for(int k=0; k<ndim; k++)
            e[k] = std::max( eval[k]+p ,0.0e0);

        // begin central alpha-chop loop

        double da, d2s, cpen, cnew, gqd, xu;
        do{

            // set alpha between a1 and a2 (allowing for large dynamic range)
            // and da for use in terminating the loop

            al=(2.0*a1*(a2/al)+a1+a2)/(2.0+(a1+a2)/al);
            da=a2-a1;

            // solve b(..)x(.)=gq(.) for x(.) and get corresponding new
            // statistics cnew, cpen and squared-distance d2s

            d2s=0.0;
            cpen=c;
            cnew=c;
            for(int k=0; k<ndim; k++){
                gqd   = al*gsd[k]-gcd[k];
                xu    = gqd/(e[k]+al);
                d2s  += xu*xu;
                cpen += xu*(gcd[k]+0.5*e[k]*xu);
                cnew += xu*(gcd[k]+0.5*eval[k]*xu);
            }

            // alpha control

            if( cpen > cmax || (cnew > cmin && (d2s-rates)*(c-caim)<=0.0)){

                // reduce alpha and store upper limit in a2 if distance is
                // acceptable

                a2=al;
                if(d2s<=rates) aa=al;

            }else{

                // increase alpha

                a1=al;
            }
        }while(a2-a1<da);

        // end alpha loop when full machine accuracy is reached
        //  and da can no longer decrease
        //
        // exit with aa=a2 if the result is within the distance limit
    }

    void meml6(const int nsrch, const int ndim, const double gsd[6],
               const double gcd[6], const double eval[6],
               const double vecu[6][6], const double s,
               const double c, const double pdist, const double alpha,
               double gqd[6], double xu[6], double w[6],
               double& snew, double& cnew, double& cpen,
               double& d2s){

        //
        // purpose:
        //   transfer solution back to original coordinates
        //
        // parameters:
        //   argument type  i/o  dimension   description
        //    nsrch   i     i      -         number of original search dirns
        //    ndim    i     i      -         number of independent eigenvectors
        //    gsd     r*8   i     ndim       components of gsd in diagonal coords
        //    gcd     r*8   i     ndim       components of gcd in diagonal coords
        //    eval    r*8   i     ndim       eigenvalues of c matrix
        //    vecu    r*8   i   nsrch,ndim   eigenvector matrix from (6,6) array
        //    s       r*8   i      -         current value of entropy
        //    c       r*8   i      -         current value of constraint
        //    pdist   r*8   i      -         distance penalty
        //    alpha   r*8   i      -         alpha
        //    gqd     r*8     o   ndim       alpha*gsd - gcd
        //    xu      r*8     o   ndim       map increment in diagonal coords
        //    w       r*8     o   nsrch      map increment in original coords
        //    snew    r*8     o    -         predicted value of entropy
        //    cnew    r*8     o    -         predicted value of constraint
        //    cpen    r*8     o    -         cnew modified by distance penalty
        //    d2s     r*8     o    -         distance**2 actually moved  (*sumf)
        //
        // globals:
        //   variable  common block  i/o  description
        //    -
        //
        // external calls:
        //    -
        //
        // history:
        //   john skilling    8 nov 1985     initial release
        //
        // notes:
        //   (1) the output parameters gqd, xu, cpen, cnew, d2s
        //       are only returned for diagnostics
        //
        // set map increment in diagonal coords and set scalar diagnostics

        double e;
        d2s=0.0;
        cnew = cpen = c;
        snew=s;
        for(int k=0; k<ndim; k++){
            e      = std::max( eval[k]+pdist ,0.0e0);
            gqd[k] = alpha*gsd[k]-gcd[k];
            xu[k]  = gqd[k]/(e+alpha);
            d2s   += xu[k]*xu[k];
            cpen  += xu[k]*(gcd[k]+0.5*e*xu[k]);
            cnew  += xu[k]*(gcd[k]+0.5*eval[k]*xu[k]);
            snew  += xu[k]*(gsd[k]-0.5*xu[k]);
        }

        // rewrite map increment in original coords

        for(int j=0; j<nsrch; j++){
            w[j]=0.0;
            for(int k=0; k<ndim; k++)
                w[j] += vecu[j][k]*xu[k];
        }
    }


    void meml66(const int nsrch, const double sdd[6][6], double w[6]){

        //
        // purpose:
        //          w[0] := coefficient of "map" if total(map) held fixed
        //
        // parameters:
        //   argument type  i/o  dimension   description
        //    nsrch   i     i      -         number of original search dirns
        //    sdd     r*8   i   nsrch,nsrch  entropy metric
        //    w       r*8   i o   nsrch      map increment in original coords
        //
        // globals:
        //   variable  common block  i/o  description
        //    -
        //
        // external calls:
        //    -
        //
        // history:
        //   john skilling    8 nov 1985     initial release
        //

        double x=0.;
        for(int i=1; i<nsrch; i++)
            x += sdd[0][i]*w[i];

        w[0] = -x/sdd[0][0];
    }

    void meml7(const int nsrch, const double w[6], const double snew,
               const double cnew, const double d2s, const float sdefa,
               float wa[6], float &snewa, float &cnewa,
               float &rnewa){

        //
        // purpose:
        //   copy to external (single precision) parameters
        //
        // parameters:
        //   argument type  i/o  dimension   description
        //    nsrch   i     i      -         number of original search dirns
        //    w       r*8   i     nsrch      map increment in original coords
        //    snew    r*8   i      -         predicted value of entropy
        //    cnew    r*8   i      -         predicted value of constraint
        //    d2s     r*8   i      -         distance moved (*sumf)
        //    sdefa   r     i      -         scale factor for entropy
        //    wa      r       o    -         map increment components
        //    snewa   r       o    -         predicted value of entropy
        //    cnewa   r       o    -         predicted value of constraint
        //    rnewa   r       o    -         dimensionless distance moved
        //
        // globals:
        //   variable  common block  i/o  description
        //    -
        //
        // external calls:
        //    -
        //
        // history:
        //   john skilling    8 nov 1985     initial release
        //
        // notes:
        //    (1) internal calculations use s scaled by sumdef
        //

        snewa = float(snew)/sdefa;
        cnewa = float(cnew);
        rnewa = sqrt(float(d2s)/sdefa);
        for(int i=0; i<nsrch; i++)
            wa[i] =float(w[i]);
    }


    void meml8(const int nsrch, const int ndim, const double sdd[6][6],
               const double eval[6], const double gsd[6],
               const double gcd[6], const double gqd[6],
               const double xu[6], const float w[6],
               const float rnewa, const float s, const double c,
               const double caim, const float snew,
               const double cnew, const double cpen,
               const double alpha, const double pdist,
               const double pen0){

        //
        // purpose:
        //   optional diagnostics
        //
        // parameters:
        //   argument type  i/o  dimension   description
        //    nsrch   i     i      -         number of original search dirns
        //    ndim    i     i      -         number of independent eigenvectors
        //    sdd     r*8   i   nsrch,nsrch  entropy metric in original coords
        //    eval    r*8   i     ndim       eigenvalues of c matrix
        //    gsd     r*8   i     ndim       components of gsd in diagonal coords
        //    gcd     r*8   i     ndim       components of gcd in diagonal coords
        //    gqd     r*8   i     ndim       alpha*gsd - gcd
        //    xu      r*8   i     ndim       map increment in diagonal coords
        //    w       r     i     nsrch      map increment in original coords
        //    rnewa   r     i      -         distance moved (dimensionless)
        //    s       r     i      -         current value of entropy
        //    c       r*8   i      -         current value of constraint
        //    caim    r*8   i      -         realistic target value of constraint
        //    snew    r     i      -         predicted value of entropy
        //    cnew    r*8   i      -         predicted value of constraint
        //    cpen    r*8   i      -         cnew modified by distance penalty
        //    alpha   r*8   i      -         alpha
        //    pdist   r*8   i      -         distance penalty
        //    pen0    r*8   i      -         penalty offset
        //
        // globals:
        //   variable  common block  i/o  description
        //    -
        //
        // external calls:
        //    -
        //
        // history:
        //   john skilling    8 nov 1985     initial release
        //                   10 dec 1985     format altered for ibm3081
        //

        double cmin=std::min(c,caim);
        double cmax=std::max(c,caim);
        std::cerr << " rcap  " << rnewa << " " << cnew
                  << " " << alpha << " " << pdist-pen0 << std::endl;
        std::cerr << " <<<<c " << cmin << " " << cnew
                  << " " << cpen << " " << cmax << " " << c << std::endl;
        std::cerr << " new s " << snew << " " << s << std::endl;

        std::cerr << " a.a ";
        for(int i=0;i<nsrch;i++) std::cerr << sdd[i][i] << " ";
        std::cerr << std::endl;

        std::cerr << " eval ";
        for(int i=0;i<ndim;i++) std::cerr << eval[i] << " ";
        std::cerr << std::endl;

        std::cerr << " gsd ";
        for(int i=0;i<ndim;i++) std::cerr << gsd[i] << " ";
        std::cerr << std::endl;

        std::cerr << " gcd ";
        for(int i=0;i<ndim;i++) std::cerr << gcd[i] << " ";
        std::cerr << std::endl;

        std::cerr << " gqd ";
        for(int i=0;i<ndim;i++) std::cerr << gqd[i] << " ";
        std::cerr << std::endl;

        std::cerr << " xu ";
        for(int i=0;i<ndim;i++) std::cerr << xu[i] << " ";
        std::cerr << std::endl;

        std::cerr << " w ";
        for(int i=0;i<nsrch;i++) std::cerr << w[i] << " ";
        std::cerr << std::endl;

        std::cerr << std::endl;
        return;
    }

    void uread(const int i){

        //
        // purpose:
        //   read a buffer from given area and increment pointers appropriately
        //
        // parameters:
        //   argument type  i/o  dimension   description
        //    i       i     i      -         input area number
        //
        // globals:
        //   variable  common block  i/o  description
        //    st          mecoms       o  storage vector
        //    mj,mk       mecomp     i    buffer sizes (map and data space)
        //    ka          mecomp     i    allocation pointers
        //    kb          mecomp      -   base addresses
        //    kc          mecomp     i o  core addresses
        //    kd          mecomp     i o  disc or upper core addresses
        //    pr          mecomc       o  read/write flags
        //
        // external calls:
        //    uget        get buffer from disc
        //
        // history:
        //   john skilling    8 nov 1985     initial release
        //

        Gbl::pr[i] = 'r';
        int m;
        if( i < 20){
            m = Gbl::mj;
        }else{
            m = Gbl::mk;
        }
        if(Gbl::ka[i]>0){

            // disc

            //      uget(i,Gbl::st+Gbl::kc[i],m);
            std::cerr << "Disc I/O not enabled" << std::endl;
            exit(EXIT_FAILURE);
        }else{

            // core

            Gbl::kc[i] = Gbl::kd[i];
            Gbl::kd[i] = Gbl::kd[i]+m;
        }
    }

    void uwrite(const int i){

        //
        // purpose:
        //   write a buffer to given area and increment pointers appropriately
        //
        // parameters:
        //   argument type  i/o  dimension   description
        //    i       i     i      -         output area number
        //
        // globals:
        //   variable  common block  i/o  description
        //    st          mecoms     i    storage vector
        //    mj,mk       mecomp     i    buffer sizes (map and data space)
        //    ka          mecomp     i    allocation pointers
        //    kb          mecomp      -   base addresses
        //    kc          mecomp     i o  core addresses
        //    kd          mecomp     i o  disc or upper core addresses
        //    pr          mecomc       o  read/write flags
        //
        // external calls:
        //    uput        put buffer to disc
        //
        // history:
        //   john skilling    8 nov 1985     initial release
        //
        // notes:
        //

        Gbl::pr[i] ='w';
        int m;
        if(i < 20){
            m = Gbl::mj;
        }else{
            m = Gbl::mk;
        }

        if(Gbl::ka[i]>0){

            // disc

            //      uput(i,Gbl::st+Gbl::kc[i],m);

            std::cerr << "Disc I/O not enabled" << std::endl;
            exit(EXIT_FAILURE);

        }else{

            // core

            Gbl::kd[i] = Gbl::kd[i]+m;
            Gbl::kc[i] = Gbl::kd[i];
        }
    }

    void uinit(){

        //
        // purpose:
        //   re-initialise all pointers, with optional diagnostics
        //
        // parameters:
        //   argument type  i/o  dimension   description
        //     -
        //
        // globals:
        //   variable  common block  i/o  description
        //    l0     mecomp     i    diagnostic stream and flag
        //    pr          mecomc     i    read/write flags
        //
        // external calls:
        //    meinit      re-initialise all pointers
        //    udiag       user's diagnostic routine
        //
        // history:
        //   john skilling    8 nov 1985     initial release
        //
        // notes:
        //

        if(Gbl::l0>=2) std::cerr << Gbl::pr << std::endl;
        meinit();
        if(Gbl::l0>=3) udiag();
    }

    void meinit(){

        //
        // purpose:
        //   re-initialise all pointers
        //
        // parameters:
        //   argument type  i/o  dimension   description
        //     -
        //
        // globals:
        //   variable  common block  i/o  description
        //    pr          mecomc       o  read/write flags
        //
        // external calls:
        //    ureset      reset pointers for a given area
        //
        // history:
        //   john skilling    8 nov 1985     initial release
        //
        // notes:
        //

        for(int i=0; i<40; i++){
            Gbl::pr[i]= '.';
            ureset(i);
        }
        Gbl::pr[40] = 0;
    }

    void ureset(const int i){

        //
        // purpose:
        //   reset pointers for given area
        //
        // parameters:
        //   argument type  i/o  dimension   description
        //    i       i     i      -         area number
        //
        // globals:
        //   variable  common block  i/o  description
        //    ka          mecomp     i    allocation pointers
        //    kb          mecomp     i    base addresses
        //    kc          mecomp       o  core addresses
        //    kd          mecomp       o  disc or upper core addresses
        //
        // external calls:
        //    uwind       rewind disc file
        //
        // history:
        //   john skilling    8 nov 1985     initial release
        //
        // notes:
        //

        Gbl::kc[i] = Gbl::kb[i];
        if(Gbl::ka[i]>0){

            // disc

            //      uwind(i);

            std::cerr << "Disc I/O not enabled" << std::endl;
            exit(EXIT_FAILURE);

        }else{

            // core

            Gbl::kd[i] = Gbl::kb[i];
        }
    }

    void uget(const int i, float x[], const int m){

        //
        // purpose:
        //   get a buffer from disc file
        //
        // parameters:
        //   argument type  i/o  dimension   description
        //    i       i     i      -         input area number
        //    x       r       o    m         buffer's destination in core
        //    m       i     i      -         length of buffer in words
        //
        // globals:
        //   variable  common block  i/o  description
        //    ka          mecomp     i    allocation pointers to disc files
        //
        // external calls:
        //     -
        //
        // history:
        //   john skilling    8 nov 1985     initial release
        //
        // notes:
        //    (1)  this implementation uses unformatted sequential fortran i/o
        //    (2)  for direct access files, pointers kd are available for use

        //    int k = Gbl::ka[i];
        //    read(k) x;
        std::cerr << "uget: SHOULD NOT BE HERE" << std::endl;
    }

    void uput(const int i, float x[], const int m){

        //
        // purpose:
        //   put a buffer to disc file
        //
        // parameters:
        //   argument type  i/o  dimension   description
        //    i       i     i      -         output area number
        //    x       r     i      m         buffer's source in core
        //    m       i     i      -         length of buffer in words
        //
        // globals:
        //   variable  common block  i/o  description
        //    ka          mecomp     i    allocation pointers to disc files
        //
        // external calls:
        //     -
        //
        // history:
        //   john skilling    8 nov 1985     initial release
        //
        // notes:
        //    (1)  this implementation uses unformatted sequential fortran i/o
        //    (2)  for direct access files, pointers kd are available for use
        //

        //    int k = Gbl::ka[i];
        //      write(k) x
        std::cerr << "uput: SHOULD NOT BE HERE" << std::endl;
    }

    void uwind(const int i){

        //
        // purpose:
        //    rewind disc file
        //
        // parameters:
        //   argument type  i/o  dimension   description
        //    i       i     i      -         area number
        //
        // globals:
        //   variable  common block  i/o  description
        //    ka          mecomp     i    allocation pointers to disc files
        //
        // external calls:
        //     -
        //
        // history:
        //   john skilling    8 nov 1985     initial release
        //
        // notes:
        //    (1)  this implementation uses unformatted fortran i/o
        //    (2)  for direct access files, pointers kd are available for use
        //

        //    int k=Gbl::ka[i];
        //      rewind k
        std::cerr << "uwind: SHOULD NOT BE HERE" << std::endl;
    }

    void udiag(){

        //
        // purpose:
        //   user's diagnostic routine, to inspect any part of any area in core.
        //
        // parameters:
        //   argument type  i/o  dimension   description
        //    -
        //
        // globals:
        //   variable  common block  i/o  description
        //    st          mecoms     i    storage vector
        //    nj,mj       mecomp     i    sizes  (map space)
        //    nk,mk       mecomp     i    sizes (data space)
        //    ka          mecomp     i    initialised allocation pointers
        //    kb          mecomp     i    initialised base addresses
        //
        // external calls:
        //     -
        //
        // history:
        //   john skilling    8 nov 1985     initial release
        //                   10 dec 1985     modified for ibm3081
        //
        // notes:
        //    (1)  this implementation is numerical, and interrogates the user
        //         via the default fortran stream (*)
        //

        int n, m;
        for(;;){
            do{
                std::cerr << " which area would you like?  (-1 exits):  ";
                while(!(std::cin >> n));
            }while(n<-1 || n>=40);
            if(n==-1) return;
            if(n<20){
                m=Gbl::mj;
                if(Gbl::ka[n]==0) m *= Gbl::nj;
            }else{
                m=Gbl::mk;
                if(Gbl::ka[n]==0) m *= Gbl::nk;
            }
            if(Gbl::ka[n]==0){
                std::cerr << "ka = " << Gbl::ka[n] << " (core file) ";
            }else{
                std::cerr << "ka = " << Gbl::ka[n] << " (core file) ";
            }
            std::cerr << Gbl::kb[n] << " base address" << std::endl;

            std::cerr << "                    " << Gbl::kb[n]+m-1 << " last address "
                      << std::endl;
            std::cerr << " first point ?  ";
            int j, k;
            std::cin >> j;
            j = Gbl::kb[n]+j-1;
            std::cerr << "  last point ?  ";
            std::cin >> k;
            k = Gbl::kb[n]+k-1;
            for(int i=j; i<=k; i++)
                std::cerr << Gbl::st[i] << " ";
            std::cerr << std::endl;
        }
    }

    void memcore(size_t mxbuf, size_t nmod, size_t ndat){

        //  Formats memsys core buffer for images and data vectors
        //
        //  Input:
        //  mxbuf = size of buffer
        //  nmod  = size of images
        //  ndat  = size of data
        //

        int jfile[20] = {0, 1, 2, 3, 4, 5, -1, -1, -1, -1,
                         -1, -1, -1, -1, -1, -1, -1, 8, 6, 7};

        int kfile[20] = {0, 1, 2, 3, 4, 5, 6, 7, -1, -1,
                         -1, -1, -1, -1, -1, -1, -1, -1, -1, -1};

        Gbl::nj = 1;
        Gbl::mj = nmod;
        int njfile = 0;
        size_t kore = 0;

        for(int i=0; i<40; i++)
            Gbl::ka[i] = 0;

        for(int i=0; i<20; i++){
            Gbl::kb[i]=-9999999;
            if(jfile[i] >= 0){
                njfile++;
                Gbl::kb[i] = nmod*jfile[i];
                kore = std::max(kore,size_t(Gbl::kb[i]));
            }
        }
        size_t korej = kore + nmod;

        //  load data pointers
        Gbl::nk = 1;
        Gbl::mk = ndat;
        int nkfile = 0;
        for(int i=0; i<20; i++){
            Gbl::kb[20+i] = -9999999;
            if(kfile[i] >= 0){
                nkfile++;
                Gbl::kb[20+i] = korej + ndat*kfile[i];
                kore = std::max(kore,size_t(Gbl::kb[20+i]));
            }
        }
        kore += ndat;
        if(kore > mxbuf){
            std::cerr << "Buffer too small!! Currently = " << mxbuf
                      << ", but need " << kore << std::endl;
            exit(EXIT_FAILURE);
        }
        std::cerr << "Mem::memcore used " << 100.*float(kore)/mxbuf
                  << "% of " << mxbuf << " total buffer." << std::endl;
    }

    size_t memsize(size_t nmod, size_t ndat){

        //  Formats memsys core buffer for images and data vectors
        //
        //  Input:
        //	mxbuf = size of buffer
        //	nmod  = size of images
        //	ndat  = size of data
        //

        int jfile[20] = {0, 1, 2, 3, 4, 5, -1, -1, -1, -1,
                         -1, -1, -1, -1, -1, -1, -1, 8, 6, 7};

        int kfile[20] = {0, 1, 2, 3, 4, 5, 6, 7, -1, -1,
                         -1, -1, -1, -1, -1, -1, -1, -1, -1, -1};

        int njfile = 0;
        size_t kore = 0;

        for(int i=0; i<20; i++){
            Gbl::kb[i]=-9999999;
            if(jfile[i] >= 0){
                njfile++;
                kore = std::max(kore,size_t(nmod*jfile[i]));
            }
        }
        size_t korej = kore + nmod;

        // load data pointers
        int nkfile = 0;
        for(int i=0; i<20; i++){
            if(kfile[i] >= 0){
                nkfile++;
                kore = std::max(kore,korej + size_t(ndat*kfile[i]));
            }
        }
        kore += ndat;
        return kore;
    }
}
