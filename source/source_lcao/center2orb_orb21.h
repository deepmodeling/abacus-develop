//=========================================================
// AUTHOR : Peize Lin
// DATE : 2016-01-24
//=========================================================

#ifndef W_ABACUS_DEVELOP_ABACUS_DEVELOP_SOURCE_MODULE_HAMILT_LCAO_HAMILT_LCAODFT_CENTER2_ORB_ORB21_H
#define W_ABACUS_DEVELOP_ABACUS_DEVELOP_SOURCE_MODULE_HAMILT_LCAO_HAMILT_LCAODFT_CENTER2_ORB_ORB21_H

#include "center2orb_orb11.h"
#include "center2orb.h"
#include "source_base/vector3.h"
#include "source_basis/module_ao/orb_atomic_lm.h"
#include "source_basis/module_ao/orb_gaunt_table.h"

#include <map>
#include <memory>
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

    static void init_radial_table(const std::vector<Orb21*>& orbitals, const bool use_gpu);
    static void init_radial_table(const std::vector<Orb21*>& orbitals, const std::set<size_t>& radials);

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

    /// All magnetic components, ordered as [mA1][mA2][mB].
    std::vector<double> cal_overlap(const ModuleBase::Vector3<double>& RA, const ModuleBase::Vector3<double>& RB) const;
    std::vector<ModuleBase::Vector3<double>> cal_grad_overlap(const ModuleBase::Vector3<double>& RA,
                                                              const ModuleBase::Vector3<double>& RB) const;

  private:
    void prepare_radial_table();
    void bind_product_orbitals();
    static std::vector<Orb11*> prepare_radial_tables(const std::vector<Orb21*>& orbitals);

    template <typename T>
    std::vector<T> cal_magnetic_components(
        const ModuleBase::Vector3<double>& RA,
        const ModuleBase::Vector3<double>& RB,
        std::vector<T> (Orb11::*evaluate)(const ModuleBase::Vector3<double>&,
                                        const ModuleBase::Vector3<double>&) const) const;

    const Numerical_Orbital_Lm& nA1;
    const Numerical_Orbital_Lm& nA2;
    const Numerical_Orbital_Lm& nB;

    const ModuleBase::Sph_Bessel_Recursive::D2* psb_ = nullptr;
    const ORB_gaunt_table& MGT;

    using ProductOrbitals = std::map<int, Numerical_Orbital_Lm>;
    std::shared_ptr<const ProductOrbitals> nA;
    std::map<int, Center2_Orb::Orb11> orb11s;
};

// this->orb11s[LA].Table_r[LAB][ir]

#endif // W_ABACUS_DEVELOP_ABACUS_DEVELOP_SOURCE_MODULE_HAMILT_LCAO_HAMILT_LCAODFT_CENTER2_ORB_ORB21_H
