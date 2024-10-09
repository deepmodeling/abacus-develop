//==========================================================
// AUTHOR : Peize Lin
// DATE : 2016-01-24
//==========================================================

#ifndef CENTER2_ORB_H
#define CENTER2_ORB_H

#include "module_base/sph_bessel_recursive.h"
#include "module_basis/module_ao/ORB_atomic_lm.h"

#include <set>

class Center2_Orb
{
  public:
    class Orb11;
    class Orb21;
    class Orb22;

    //=========================================================================
    // The following functions used to be in module_ao/ORB_table_phi.h, but
    // are moved here to decouple module_ao from Center2_Orb
    //=========================================================================

    static int get_rmesh(const double& R1, const double& R2, const double dr);

    static void init_Lmax(const int orb_num,
                          const int mode,
                          int& Lmax_used,
                          int& Lmax,
                          const int& Lmax_exx,
                          const int lmax_orb,
                          const int lmax_beta);

    static void init_Table_Spherical_Bessel(const int orb_num,
                                            const int mode,
                                            int& Lmax_used,
                                            int& Lmax,
                                            const int& Lmax_exx,
                                            const int lmax_orb,
                                            const int lmax_beta,
                                            const double dr,
                                            const double dk,
                                            const int kmesh,
                                            const int Rmesh,
                                            ModuleBase::Sph_Bessel_Recursive::D2*& psb);

    static void cal_ST_Phi12_R(const int& job,
                               const int& l,
                               const Numerical_Orbital_Lm& n1,
                               const Numerical_Orbital_Lm& n2,
                               const int& rmesh,
                               double* rs,
                               double* drs,
                               const ModuleBase::Sph_Bessel_Recursive::D2* psb);

    // Peize Lin add 2017-10-13
    static void cal_ST_Phi12_R(const int& job,
                               const int& l,
                               const Numerical_Orbital_Lm& n1,
                               const Numerical_Orbital_Lm& n2,
                               const std::set<size_t>& radials, // only calculate ir in radials
                               double* rs,
                               double* drs,
                               const ModuleBase::Sph_Bessel_Recursive::D2* psb);
};

#endif // CENTER2_ORB_H
