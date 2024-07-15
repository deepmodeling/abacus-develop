//=========================================================
// AUTHOR : Peize Lin
// DATE : 2016-01-24
//=========================================================

#ifndef W_ABACUS_DEVELOP_ABACUS_DEVELOP_SOURCE_MODULE_HAMILT_LCAO_HAMILT_LCAODFT_CENTER2_ORB_ORB21_H
#define W_ABACUS_DEVELOP_ABACUS_DEVELOP_SOURCE_MODULE_HAMILT_LCAO_HAMILT_LCAODFT_CENTER2_ORB_ORB21_H

#include "center2_orb-orb11.h"
#include "center2_orb.h"
#include "module_base/vector3.h"
#include "module_basis/module_ao/ORB_atomic_lm.h"
#include "module_basis/module_ao/ORB_gaunt_table.h"

#include <map>
#include <set>
#include <vector>

class Center2_Orb::Orb21
{
  public:
    Orb21(const Numerical_Orbital_Lm& nA1_in,
          const Numerical_Orbital_Lm& nA2_in,
          const Numerical_Orbital_Lm& nB_in,
          const ModuleBase::Sph_Bessel_Recursive::D2* psb,
          const ORB_gaunt_table& MGT_in);

    void init_radial_table();
    void init_radial_table(const std::set<size_t>& radials); // unit: Bohr/MOT.dr

    double cal_overlap(const ModuleBase::Vector3<double>& RA,
                       const ModuleBase::Vector3<double>& RB, // unit: Bohr
                       const int& mA1,
                       const int& mA2,
                       const int& mB) const;

    ModuleBase::Vector3<double> cal_grad_overlap(const ModuleBase::Vector3<double>& RA,
                                                 const ModuleBase::Vector3<double>& RB, // unit: Bohr
                                                 const int& mA1,
                                                 const int& mA2,
                                                 const int& mB) const;

  private:
    const Numerical_Orbital_Lm& nA1;
    const Numerical_Orbital_Lm& nA2;
    const Numerical_Orbital_Lm& nB;

    const ModuleBase::Sph_Bessel_Recursive::D2* psb_;
    const ORB_gaunt_table& MGT;

    std::map<int, Numerical_Orbital_Lm> nA;
    std::map<int, Center2_Orb::Orb11> orb11s;
};

// this->orb11s[LA].Table_r[LAB][ir]

#endif // W_ABACUS_DEVELOP_ABACUS_DEVELOP_SOURCE_MODULE_HAMILT_LCAO_HAMILT_LCAODFT_CENTER2_ORB_ORB21_H
