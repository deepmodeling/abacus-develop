#ifndef CAL_R_OVERLAP_R_H
#define CAL_R_OVERLAP_R_H

#include "module_base/abfs-vector3_order.h"
#include "module_base/sph_bessel_recursive.h"
#include "module_base/vector3.h"
#include "module_base/ylm.h"
#include "module_basis/module_ao/ORB_atomic_lm.h"
#include "module_basis/module_ao/ORB_gaunt_table.h"
#include "module_basis/module_ao/ORB_read.h"
#include "module_basis/module_ao/parallel_orbitals.h"
#include "module_hamilt_lcao/hamilt_lcaodft/center2_orb-orb11.h"
#include "module_hamilt_lcao/hamilt_lcaodft/center2_orb-orb21.h"
#include "module_hamilt_lcao/hamilt_lcaodft/center2_orb.h"
#include "single_R_io.h"
#include "module_cell/unitcell.h"
#include <map>
#include <set>
#include <vector>

// output r_R matrix, added by Jingan
class cal_r_overlap_R
{

  public:
    cal_r_overlap_R();
    ~cal_r_overlap_R();

    double kmesh_times = 4;
    double sparse_threshold = 1e-10;
    bool binary = false;

    void init(const UnitCell& ucell,const Parallel_Orbitals& pv, const LCAO_Orbitals& orb);
    void out_rR(const UnitCell& ucell, const int& istep);
    void out_rR_other(const UnitCell& ucell, const int& istep, const std::set<Abfs::Vector3_Order<int>>& output_R_coor);

  private:
    void initialize_orb_table(const UnitCell& ucell, const LCAO_Orbitals& orb);
    void construct_orbs_and_orb_r(const UnitCell& ucell,const LCAO_Orbitals& orb);

    std::vector<int> iw2ia;
    std::vector<int> iw2iL;
    std::vector<int> iw2im;
    std::vector<int> iw2iN;
    std::vector<int> iw2it;

    ModuleBase::Sph_Bessel_Recursive::D2* psb_ = nullptr;
    ORB_gaunt_table MGT;

    Numerical_Orbital_Lm orb_r;
    std::vector<std::vector<std::vector<Numerical_Orbital_Lm>>> orbs;

    std::map<
        size_t,
        std::map<size_t, std::map<size_t, std::map<size_t, std::map<size_t, std::map<size_t, Center2_Orb::Orb11>>>>>>
        center2_orb11;

    std::map<
        size_t,
        std::map<size_t, std::map<size_t, std::map<size_t, std::map<size_t, std::map<size_t, Center2_Orb::Orb21>>>>>>
        center2_orb21_r;

    const Parallel_Orbitals* ParaV;
};
#endif
