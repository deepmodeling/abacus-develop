#ifndef CAL_R_OVERLAP_R_H
#define CAL_R_OVERLAP_R_H

#include "source_lcao/module_ri/abfs-vector3_order.h"
#include "source_base/sph_bessel_recursive.h"
#include "source_base/vector3.h"
#include "source_base/ylm.h"
#include "source_basis/module_ao/ORB_atomic_lm.h"
#include "source_basis/module_ao/ORB_gaunt_table.h"
#include "source_basis/module_ao/ORB_read.h"
#include "source_basis/module_ao/parallel_orbitals.h"
#include "source_cell/module_neighbor/sltk_grid_driver.h"
#include "source_cell/unitcell.h"
#include "source_lcao/center2_orb-orb11.h"
#include "source_lcao/center2_orb-orb21.h"
#include "source_lcao/center2_orb.h"

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
    void init_nonlocal(const UnitCell& ucell,const Parallel_Orbitals& pv, const LCAO_Orbitals& orb);
    // Initialize the r-weighted two-center integral tables between numerical atomic
    // orbitals and DeePKS projected orbitals (orb.Alpha[0]). Used by the rt-TDDFT
    // current operator to evaluate the commutator -i[V_delta, r]. Added for DeePKS.
    void init_alpha(const UnitCell& ucell,const Parallel_Orbitals& pv, const LCAO_Orbitals& orb);
    ModuleBase::Vector3<double> get_psi_r_psi(
      const ModuleBase::Vector3<double>& R1,
      const int& T1,
      const int& L1,
      const int& m1,
      const int& N1,
      const ModuleBase::Vector3<double>& R2,
      const int& T2,
      const int& L2,
      const int& m2,
      const int& N2
    );
    void get_psi_r_beta(
      const UnitCell& ucell,
      std::vector<std::vector<double>>& nlm,
      const ModuleBase::Vector3<double>& R1,
      const int& T1,
      const int& L1,
      const int& m1,
      const int& N1,
      const ModuleBase::Vector3<double>& R2,
      const int& T2
    );
    // Compute <phi|alpha> and <phi|r|alpha> (4 components: nlm[0]=<phi|alpha>,
    // nlm[1..3]=<phi|x,y,z|alpha>) between an orbital phi (T1,L1,m1,N1 at R1) and all
    // DeePKS projected orbitals alpha (orb.Alpha[0]) centered at R2. The alpha index
    // ordering (L0->N0->running-m) matches TwoCenterIntegrator::snap, hence phialpha
    // and gedm. Added for DeePKS rt-TDDFT current.
    void get_psi_r_alpha(
      const LCAO_Orbitals& orb,
      std::vector<std::vector<double>>& nlm,
      const ModuleBase::Vector3<double>& R1,
      const int& T1,
      const int& L1,
      const int& m1,
      const int& N1,
      const ModuleBase::Vector3<double>& R2
    );
    void out_rR(const UnitCell& ucell, const Grid_Driver& gd, const int& istep);
    void out_rR_other(const UnitCell& ucell, const int& istep, const std::set<Abfs::Vector3_Order<int>>& output_R_coor);

  private:
    // lmax_extra (default 0) enlarges the spherical-Bessel / Gaunt tables to also
    // cover integrals against a higher-L projector set (e.g. DeePKS alpha orbitals).
    // With the default it reduces exactly to the orbital-only sizing, so existing
    // callers (init / init_nonlocal) are byte-for-byte unaffected.
    void initialize_orb_table(const UnitCell& ucell, const LCAO_Orbitals& orb, const int lmax_extra = 0);
    void construct_orbs_and_orb_r(const UnitCell& ucell,const LCAO_Orbitals& orb);
    void construct_orbs_and_nonlocal_and_orb_r(const UnitCell& ucell,const LCAO_Orbitals& orb);
    void construct_orbs_and_alpha_and_orb_r(const UnitCell& ucell,const LCAO_Orbitals& orb);

    std::vector<int> iw2ia;
    std::vector<int> iw2iL;
    std::vector<int> iw2im;
    std::vector<int> iw2iN;
    std::vector<int> iw2it;

    ModuleBase::Sph_Bessel_Recursive::D2* psb_ = nullptr;
    ORB_gaunt_table MGT;

    Numerical_Orbital_Lm orb_r;
    std::vector<std::vector<std::vector<Numerical_Orbital_Lm>>> orbs;
    std::vector<std::vector<Numerical_Orbital_Lm>> orbs_nonlocal;
    // DeePKS projected orbitals (orb.Alpha[0]), indexed [L_alpha][N_alpha].
    std::vector<std::vector<Numerical_Orbital_Lm>> orbs_alpha;

    std::map<
        size_t,
        std::map<size_t, std::map<size_t, std::map<size_t, std::map<size_t, std::map<size_t, Center2_Orb::Orb11>>>>>>
        center2_orb11;

    std::map<
        size_t,
        std::map<size_t, std::map<size_t, std::map<size_t, std::map<size_t, std::map<size_t, Center2_Orb::Orb21>>>>>>
        center2_orb21_r;

    std::map<
        size_t,
        std::map<size_t, std::map<size_t, std::map<size_t, std::map<size_t, Center2_Orb::Orb11>>>>>
        center2_orb11_nonlocal;

    std::map<
        size_t,
        std::map<size_t, std::map<size_t, std::map<size_t, std::map<size_t, Center2_Orb::Orb21>>>>>
        center2_orb21_r_nonlocal;

    // r-weighted tables between orbitals and DeePKS projectors (orb.Alpha[0]).
    // Keyed [T1][L1][N1][L_alpha] -> {N_alpha -> Orb}. Alpha is type-independent,
    // so there is no ket-type (TB) index; the ket is indexed by (L_alpha, N_alpha).
    std::map<
        size_t,
        std::map<size_t, std::map<size_t, std::map<size_t, std::map<size_t, Center2_Orb::Orb11>>>>>
        center2_orb11_alpha;

    std::map<
        size_t,
        std::map<size_t, std::map<size_t, std::map<size_t, std::map<size_t, Center2_Orb::Orb21>>>>>
        center2_orb21_r_alpha;

    const Parallel_Orbitals* ParaV = nullptr;
};
#endif
