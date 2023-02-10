#ifndef SOC_H
#define SOC_H

#include "module_base/global_function.h"
#include "module_base/global_variable.h"

class Fcoef
{

  public:
    Fcoef(){};
    ~Fcoef();

    int ind1 = 1;
    int ind4 = 1;
    int ind5 = 1;

    std::complex<double> *p = nullptr;

    inline std::complex<double> &operator()(const int &i1, const int &i2, const int &i3, const int &i4, const int &i5)
    {
        return p[ind2 * ind3 * ind4 * ind5 * i1 + ind3 * ind4 * ind5 * i2 + ind4 * ind5 * i3 + ind5 * i4 + i5];
    }

    inline const std::complex<double> &operator()(const int &i1,
                                                  const int &i2,
                                                  const int &i3,
                                                  const int &i4,
                                                  const int &i5) const
    {
        return p[ind2 * ind3 * ind4 * ind5 * i1 + ind3 * ind4 * ind5 * i2 + ind4 * ind5 * i3 + ind5 * i4 + i5];
    }

    void create(const int i1, const int i2, const int i3);

    //   void free();

  private:
    int ind2 = 2;
    int ind3 = 2;
};

//-----------------------
// spin-orbital coupling
//-----------------------
class Soc
{

  public:
    Soc(){};
    ~Soc();

    double spinor(const int l, const double j, const int m, const int spin);

    int sph_ind(const int l, const double j, const int m, const int spin);

    void rot_ylm(const int lmax);
    // std::complex<double> **rotylm;

    const std::complex<double> &rotylm(const int &i1, const int &i2) const
    {
        return p_rot[l2plus1_ * i1 + i2];
    }

    Fcoef fcoef;
    void set_fcoef(const int &l1,
                   const int &l2,
                   const int &is1,
                   const int &is2,
                   const int &m1,
                   const int &m2,
                   const double &j1,
                   const double &j2,
                   const int &it,
                   const int &ip1,
                   const int &ip2);

    // int npol;
  private:
    int l_max_ = -1; // maximum orbital index, must initialized before set_fcoef
    int l2plus1_ = -1;
    std::complex<double> *p_rot = nullptr;
};
#endif
