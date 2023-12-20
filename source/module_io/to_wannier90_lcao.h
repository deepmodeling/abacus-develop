#ifndef TOWannier90_LCAO_H
#define TOWannier90_LCAO_H

#include <iostream>
#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <vector>

#include "to_wannier90.h"

#include "module_base/complexmatrix.h"
#include "module_base/global_function.h"
#include "module_base/global_variable.h"
#include "module_base/lapack_connector.h"
#include "module_base/matrix.h"
#include "module_base/matrix3.h"
#include "module_cell/klist.h"
#include "module_hamilt_lcao/hamilt_lcaodft/wavefunc_in_pw.h"
#include "module_psi/psi.h"
#include <map>
#include <set>

#include "module_hamilt_lcao/hamilt_lcaodft/center2_orb.h"
#include "module_hamilt_lcao/hamilt_lcaodft/center2_orb-orb11.h"
#include "module_hamilt_lcao/hamilt_lcaodft/center2_orb-orb21.h"

#include "module_basis/module_ao/ORB_table_phi.h"
#include "module_basis/module_ao/ORB_gaunt_table.h"
#include "module_basis/module_ao/ORB_atomic_lm.h"
#include "module_basis/module_ao/ORB_read.h"
#include "module_basis/module_ao/parallel_orbitals.h"
#include "module_base/vector3.h"
#include "module_base/abfs-vector3_order.h"
#include "module_base/ylm.h"
#include "single_R_io.h"

#include "module_base/parallel_reduce.h"
#include "module_base/timer.h"
#include "module_cell/module_neighbor/sltk_grid_driver.h"
#include "module_hamilt_pw/hamilt_pwdft/global.h"

#ifdef __LCAO
#include "module_hamilt_lcao/hamilt_lcaodft/local_orbital_wfc.h"
#include "module_hamilt_lcao/module_gint/grid_technique.h"


#include "module_base/abfs-vector3_order.h"
#include "module_base/math_lebedev_laikov.h"
#include "module_hamilt_lcao/module_hcontainer/hcontainer.h"
#include "fR_overlap.h"

class Coordinate_3D
{
public:
    double x = 0;
    double y = 0;
    double z = 0;

    Coordinate_3D(double x=0.0, double y=0.0, double z=0.0): x(x), y(y), z(z) {}

    bool operator<(const Coordinate_3D& other) const 
    {
        const double threshold = 1e-8;
        if (std::abs(x - other.x) >= threshold) return x < other.x;
        if (std::abs(y - other.y) >= threshold) return y < other.y;
        return std::abs(z - other.z) >= threshold && z < other.z;
    }
};

class toWannier90_LCAO : public toWannier90
{
  public:
    toWannier90_LCAO(
      const bool &out_wannier_mmn, 
      const bool &out_wannier_amn, 
      const bool &out_wannier_unk, 
      const bool &out_wannier_eig,
      const bool &out_wannier_wvfn_formatted, 
      const std::string &nnkpfile,
      const std::string &wannier_spin
    );
    ~toWannier90_LCAO();

    void calculate(
      const ModuleBase::matrix& ekb,
      const K_Vectors& kv,
      const psi::Psi<std::complex<double>>& psi, 
      const Parallel_Orbitals *pv
    );

    void calculate(
      const ModuleBase::matrix& ekb,
      const K_Vectors& kv,
      const psi::Psi<double>& psi, 
      const Parallel_Orbitals *pv
    )
    {
      throw std::logic_error("The wave function of toWannier90_LCAO_IN_PW is generally a std::complex<double> type.");
    }

    void cal_Amn(const K_Vectors& kv, const psi::Psi<std::complex<double>>& psi);
    void cal_Mmn(const K_Vectors& kv, const psi::Psi<std::complex<double>>& psi);
    void out_unk(const psi::Psi<std::complex<double>>& psi);

  protected:
    // Radial section of trial orbitals
    const int mesh_r = 1001; // unit is a.u.
    const double dr = 0.01;
    std::vector<std::vector<Numerical_Orbital_Lm>> A_orbs;

    // Use default element orbital information
    int orb_r_ntype = 0;

    ORB_table_phi MOT;
    ORB_gaunt_table MGT;
    double kmesh_times = 1;

    Numerical_Orbital_Lm orb_r;
    std::vector<std::vector<std::vector<Numerical_Orbital_Lm>>> orbs;
    std::map<size_t, std::map<size_t, std::map<size_t, Center2_Orb::Orb11>>> center2_orb11_A;

    std::vector<ModuleBase::Vector3<double>> R_coor_car;
    std::vector<std::vector<std::vector<double>>> psi_psiA_R;

    std::vector<int> iw2it;
    std::vector<int> iw2ia;
    std::vector<int> iw2iL;
    std::vector<int> iw2iN;
    std::vector<int> iw2im;
    std::vector<int> iw2iorb;

    const Parallel_Orbitals* ParaV;

    void initialize_orb_table();
    void produce_basis_orb();
    void set_R_coor();
    void count_delta_k(const K_Vectors& kv);

    std::vector<Coordinate_3D> delta_k_all;
    std::map<Coordinate_3D, int> delta_k_all_index;

    void unkdotkb(
      const K_Vectors& kv, 
      const psi::Psi<std::complex<double>>& psi_in, 
      const int& ik, 
      const int& ikb, 
      const ModuleBase::Vector3<double> G, 
      ModuleBase::ComplexMatrix &Mmn
    );

    void produce_trial_in_lcao();
    void construct_overlap_table_project();
    void cal_orbA_overlap_R();

    void unkdotA(
      const K_Vectors& kv, 
      const psi::Psi<std::complex<double>>& psi_in, 
      const int& ik, 
      ModuleBase::ComplexMatrix &Amn
    );

    std::vector<FR_overlap<std::complex<double>>> FR;
};
#endif
#endif
