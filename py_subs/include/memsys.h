#ifndef TRM_MEMSYS_H
#define TRM_MEMSYS_H

//! Header file for memsys
namespace Mem {

  //! Namespace of global variables
  namespace Gbl {
    int nj,mj,nk,mk,ka[40],kb[40],kc[40],kd[40];
    int l0,l1,m0,m10,m11,m20,m21,m3;
    float *st;
    char pr[41];
  }

  //! Basic mem function
  void memprm(const int method, const int level, const float aim,
              const float rmax, const float def, const float acc,
              float &c, float &test, float &cnew, float &s,
              float &rnew, float &snew, float &sumf);

  //! Standard opus declaration
  void opus(const int k, const int l);

  //! Standard tropus declaration
  void tropus(const int k, const int l);

  //! Sets up the memory
  void memcore(size_t mxbuf, size_t nmod, size_t ndat);

  //! Computes buffer size needed
  size_t memsize(size_t nmod, size_t ndat);

}

#endif
