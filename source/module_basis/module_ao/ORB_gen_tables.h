#ifndef W_ABACUS_DEVELOP_ABACUS_DEVELOP_SOURCE_MODULE_BASIS_MODULE_AO_ORB_GEN_TABLES_H
#define W_ABACUS_DEVELOP_ABACUS_DEVELOP_SOURCE_MODULE_BASIS_MODULE_AO_ORB_GEN_TABLES_H

#include "ORB_gaunt_table.h"
#include "ORB_read.h"
#include "ORB_table_phi.h"
#include "module_base/complexarray.h"
#include "module_base/intarray.h"
#include "module_base/matrix.h"
#include "module_base/sph_bessel_recursive.h"
#include "module_base/vector3.h"
#include "module_basis/module_nao/two_center_bundle.h"
#include "module_cell/setup_nonlocal.h"

#include <memory>

/// used to be 'Use_Overlap_Table',
/// now the name is 'ORB_gen_tables'
class ORB_gen_tables
{
  public:
    friend class ORB_control;

    ORB_gen_tables();
    ~ORB_gen_tables();

    void gen_tables(std::ofstream& ofs_in, // mohan add 2021-05-07
                    LCAO_Orbitals& orb,
                    const int& Lmax_exx,
                    const bool& deepks_setorb, ///<[in] whether to generate descriptors
                    const int& nprojmax,
                    const int* nproj,
                    const Numerical_Nonlocal* beta_);
    void set_unit(const double& v)
    {
        lat0 = v;
    }

  private:
    ModuleBase::Sph_Bessel_Recursive::D2* psb_ = nullptr;
    ORB_gaunt_table MGT;

    double get_distance(const ModuleBase::Vector3<double>& R1, const ModuleBase::Vector3<double>& R2) const;

    double lat0;
};

#endif
