#ifndef W_ABACUS_DEVELOP_ABACUS_DEVELOP_SOURCE_MODULE_BASIS_MODULE_AO_ORB_GEN_TABLES_H
#define W_ABACUS_DEVELOP_ABACUS_DEVELOP_SOURCE_MODULE_BASIS_MODULE_AO_ORB_GEN_TABLES_H

#include "ORB_gaunt_table.h"
#include "ORB_read.h"
#include "ORB_table_phi.h"
#include "module_base/complexarray.h"
#include "module_base/intarray.h"
#include "module_base/matrix.h"
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

    void snap_psipsi(const LCAO_Orbitals& orb,
                     double olm[],
                     const int& job,    ///<[in]0 for matrix element of either S or T, 1 for its derivatives
                     const char& dtype, ///<[in] derivative type, 'S' for overlap, 'T' for kinetic energy, 'D' for
                                        ///< descriptor in deepks
                     const ModuleBase::Vector3<double>& R1,
                     const int& I1,
                     const int& l1,
                     const int& m1,
                     const int& n1,
                     const ModuleBase::Vector3<double>& R2,
                     const int& I2,
                     const int& l2,
                     const int& m2,
                     const int& n2,
                     bool cal_syns = false,
                     double dmax = 0.0) const;

    /// we need to destroy the tables: SR,TR,NR
    /// after ionic optimization is done.
    ORB_table_phi MOT;

  private:
    ORB_gaunt_table MGT;

    double get_distance(const ModuleBase::Vector3<double>& R1, const ModuleBase::Vector3<double>& R2) const;

    double lat0;
};

#endif
