#ifndef LCAO_DEEPKS_H
#define LCAO_DEEPKS_H

#ifdef __DEEPKS

#include "deepks_basic.h"
#include "deepks_descriptor.h"
#include "deepks_force.h"
#include "deepks_fpre.h"
#include "deepks_hmat.h"
#include "deepks_orbital.h"
#include "deepks_orbpre.h"
#include "deepks_pdm.h"
#include "deepks_phialpha.h"
#include "deepks_spre.h"
#include "deepks_vdelta.h"
#include "deepks_vdpre.h"
#include "module_base/complexmatrix.h"
#include "module_base/intarray.h"
#include "module_base/matrix.h"
#include "module_base/timer.h"
#include "module_basis/module_ao/parallel_orbitals.h"
#include "module_basis/module_nao/two_center_integrator.h"
#include "module_cell/module_neighbor/sltk_grid_driver.h"
#include "module_hamilt_lcao/module_hcontainer/hcontainer.h"
#include "module_io/winput.h"

#include <torch/script.h>
#include <torch/torch.h>

///
/// The LCAO_Deepks contains subroutines for implementation of the DeePKS method in atomic basis.
/// In essential, it is a machine-learned correction term to the XC potential
/// in the form of delta_V=|alpha> V(D) <alpha|, where D is a list of descriptors
/// The subroutines may be roughly grouped into 3 types
/// 1. generation of projected density matrices pdm=sum_i,occ <phi_i|alpha><alpha|phi_i>
///    and then descriptors D=eig(pdm)
///    as well as their gradients with regard to atomic position, gdmx = d/dX (pdm)
///    and grad_vx = d/dX (D)
/// 2. loading the model, which requires interfaces with libtorch
/// 3. applying the correction potential, delta_V, in Kohn-Sham Hamiltonian and calculation of energy, force, stress
///
/// For details of DeePKS method, you can refer to [DeePKS paper](https://pubs.acs.org/doi/10.1021/acs.jctc.0c00872).
///
///
// caoyu add 2021-03-29
// wenfei modified 2022-1-5
//
class LCAO_Deepks
{

    //-------------------
    // public variables
    //-------------------
  public:
    ///(Unit: Ry) Correction energy provided by NN
    double E_delta = 0.0;
    ///(Unit: Ry)  \f$tr(\rho H_\delta), \rho = \sum_i{c_{i, \mu}c_{i,\nu}} \f$ (for gamma_only)
    double e_delta_band = 0.0;

    /// Correction term to the Hamiltonian matrix: \f$\langle\phi|V_\delta|\phi\rangle\f$ (for gamma only)
    /// The size of first dimension is 1, which is used for the consitence with H_V_delta_k
    std::vector<std::vector<double>> H_V_delta;
    /// Correction term to Hamiltonian, for multi-k
    std::vector<std::vector<std::complex<double>>> H_V_delta_k;

    // functions for hr status: 1. get value; 2. set value;
    int get_hr_cal()
    {
        return this->hr_cal;
    }
    void set_hr_cal(bool cal)
    {
        this->hr_cal = cal;
    }

    // temporary add two getters for inl_index and gedm
    int get_inl(const int& T0, const int& I0, const int& L0, const int& N0)
    {
        return inl_index[T0](I0, L0, N0);
    }
    const double* get_gedms(const int& inl)
    {
        return gedm[inl];
    }

    int get_lmaxd()
    {
        return lmaxd;
    }
    //-------------------
    // private variables
    //-------------------
    //  private:
  public:                              // change to public to reconstuct the code, 2024-07-22 by mohan
    int lmaxd = 0;                     // max l of descirptors
    int nmaxd = 0;                     //#. descriptors per l
    int inlmax = 0;                    // tot. number {i,n,l} - atom, n, l
    int n_descriptor;                  // natoms * des_per_atom, size of descriptor(projector) basis set
    int des_per_atom;                  // \sum_L{Nchi(L)*(2L+1)}
    int* inl_l;                        // inl_l[inl_index] = l of descriptor with inl_index
    ModuleBase::IntArray* alpha_index; // seems not used in the code
    ModuleBase::IntArray* inl_index;   // caoyu add 2021-05-07

    bool init_pdm = false; // for DeePKS NSCF calculation, set init_pdm to skip the calculation of pdm in SCF iteration

    // deep neural network module that provides corrected Hamiltonian term and
    // related derivatives. Used in cal_gedm.
    torch::jit::script::Module model_deepks;

    // saves <phi(0)|alpha(R)> and its derivatives
    // index 0 for itself and index 1-3 for derivatives over x,y,z
    std::vector<hamilt::HContainer<double>*> phialpha;

    // projected density matrix
    // [tot_Inl][2l+1][2l+1], here l is corresponding to inl;
    // [nat][nlm*nlm] for equivariant version
    std::vector<torch::Tensor> pdm;

    /// dE/dD, autograd from loaded model(E: Ry)
    double** gedm; //[tot_Inl][(2l+1)*(2l+1)]

    // HR status,
    // true : HR should be calculated
    // false : HR has been calculated
    bool hr_cal = true;

    //-------------------
    // subroutines, grouped according to the file they are in:
    //-------------------

    //-------------------
    // LCAO_deepks.cpp
    //-------------------

    // This file contains constructor and destructor of the class LCAO_deepks,
    // as well as subroutines for initializing and releasing relevant data structures

    // Other than the constructor and the destructor, it contains 3 types of subroutines:
    // 1. subroutines that are related to calculating descriptors:
    //   - init : allocates some arrays
    //   - init_index : records the index (inl)
    // 2. subroutines that are related to V_delta:
    //   - allocate_V_delta : allocates H_V_delta; if calculating force, it also allocates F_delta

  public:
    explicit LCAO_Deepks();
    ~LCAO_Deepks();

    /// Allocate memory and calculate the index of descriptor in all atoms.
    ///(only for descriptor part, not including scf)
    void init(const LCAO_Orbitals& orb,
              const int nat,
              const int ntype,
              const int nks,
              const Parallel_Orbitals& pv_in,
              std::vector<int> na);

    /// Allocate memory for correction to Hamiltonian
    void allocate_V_delta(const int nat, const int nks = 1);

    //! a temporary interface for cal_e_delta_band
    template <typename TK>
    void dpks_cal_e_delta_band(const std::vector<std::vector<TK>>& dm, const int nks);

  private:
    // arrange index of descriptor in all atoms
    void init_index(const int ntype, const int nat, std::vector<int> na, const int tot_inl, const LCAO_Orbitals& orb);

    const Parallel_Orbitals* pv;
};

namespace GlobalC
{
extern LCAO_Deepks ld;
}

#endif
#endif
