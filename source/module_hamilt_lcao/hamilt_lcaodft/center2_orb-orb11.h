//=========================================================
// AUTHOR : Peize Lin
// DATE : 2016-01-24
//=========================================================

#ifndef W_ABACUS_DEVELOP_ABACUS_DEVELOP_SOURCE_MODULE_HAMILT_LCAO_HAMILT_LCAODFT_CENTER2_ORB_ORB11_H
#define W_ABACUS_DEVELOP_ABACUS_DEVELOP_SOURCE_MODULE_HAMILT_LCAO_HAMILT_LCAODFT_CENTER2_ORB_ORB11_H

#include "center2_orb.h"
#include "module_base/sph_bessel_recursive.h"
#include "module_base/vector3.h"
#include "module_basis/module_ao/ORB_atomic_lm.h"
#include "module_basis/module_ao/ORB_gaunt_table.h"

#include <map>
#include <set>
#include <vector>

class Center2_Orb::Orb11
{
  public:
    Orb11(const Numerical_Orbital_Lm& nA_in,
          const Numerical_Orbital_Lm& nB_in,
          const ModuleBase::Sph_Bessel_Recursive::D2* psb,
          const ORB_gaunt_table& MGT_in);

    void init_radial_table();
    void init_radial_table(const std::set<size_t>& radials); // unit: Bohr/MOT.dr

    double cal_overlap(const ModuleBase::Vector3<double>& RA,
                       const ModuleBase::Vector3<double>& RB, // unit: Bohr
                       const int& mA,
                       const int& mB) const;

    ModuleBase::Vector3<double> cal_grad_overlap( // caoyu add 2021-11-19
        const ModuleBase::Vector3<double>& RA,
        const ModuleBase::Vector3<double>& RB, // unit: Bohr
        const int& mA,
        const int& mB) const;

  private:
    const Numerical_Orbital_Lm& nA;
    const Numerical_Orbital_Lm& nB;

    const ModuleBase::Sph_Bessel_Recursive::D2* psb_;
    const ORB_gaunt_table& MGT;

    std::map<int, std::vector<double>> Table_r; // unit: Bohr/MOT.dr
    std::map<int, std::vector<double>> Table_dr;

    // NOTE: the following variable is introduced during refactoring in order
    // to decouple this class with ORB_table_phi. Before, MOT.dr is widely
    // used by this class. MOT.dr, however, is in fact ORB.dR, which is in
    // turn lcao_dr, one of the four parameters that define the two-center
    // integration grid. This parameter is usually not set and defaults to
    // 0.01.
    const double dr_ = 0.01;
};

// this->Table_r[LAB][ir]

#endif // W_ABACUS_DEVELOP_ABACUS_DEVELOP_SOURCE_MODULE_HAMILT_LCAO_HAMILT_LCAODFT_CENTER2_ORB_ORB11_H
