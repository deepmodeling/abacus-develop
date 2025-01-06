#ifndef LCAO_DEEPKS_H
#define LCAO_DEEPKS_H

#ifdef __DEEPKS

#include "deepks_descriptor.h"
#include "deepks_force.h"
#include "deepks_hmat.h"
#include "deepks_orbital.h"
#include "deepks_orbpre.h"
#include "deepks_vdpre.h"
#include "module_base/complexmatrix.h"
#include "module_base/intarray.h"
#include "module_base/matrix.h"
#include "module_base/timer.h"
#include "module_basis/module_ao/parallel_orbitals.h"
#include "module_basis/module_nao/two_center_integrator.h"
#include "module_cell/module_neighbor/sltk_grid_driver.h"
#include "module_elecstate/module_dm/density_matrix.h"
#include "module_hamilt_lcao/module_hcontainer/hcontainer.h"
#include "module_io/winput.h"

#include <torch/script.h>
#include <torch/torch.h>
#include <unordered_map>

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
  public:           // change to public to reconstuct the code, 2024-07-22 by mohan
    int lmaxd = 0;  // max l of descirptors
    int nmaxd = 0;  //#. descriptors per l
    int inlmax = 0; // tot. number {i,n,l} - atom, n, l
    int nat_gdm = 0;
    int nks_V_delta = 0;

    bool init_pdm = false; // for DeePKS NSCF calculation

    // deep neural network module that provides corrected Hamiltonian term and
    // related derivatives.
    torch::jit::script::Module module;

    // saves <phi(0)|alpha(R)> and its derivatives
    // index 0 for itself and index 1-3 for derivatives over x,y,z
    std::vector<hamilt::HContainer<double>*> phialpha;

    // projected density matrix
    // [tot_Inl][2l+1][2l+1], here l is corresponding to inl;
    // [nat][nlm*nlm] for equivariant version
    std::vector<torch::Tensor> pdm;

    // gedm:dE/dD, [tot_Inl][2l+1][2l+1]	(E: Hartree)
    std::vector<torch::Tensor> gedm_tensor;

    /// dE/dD, autograd from loaded model(E: Ry)
    double** gedm; //[tot_Inl][2l+1][2l+1]

    /// size of descriptor(projector) basis set
    int n_descriptor;

    // \sum_L{Nchi(L)*(2L+1)}
    int des_per_atom;

    ModuleBase::IntArray* alpha_index;
    ModuleBase::IntArray* inl_index; // caoyu add 2021-05-07
    int* inl_l;                      // inl_l[inl_index] = l of descriptor with inl_index

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

  private:
    // arrange index of descriptor in all atoms
    void init_index(const int ntype, const int nat, std::vector<int> na, const int tot_inl, const LCAO_Orbitals& orb);

    //-------------------
    // LCAO_deepks_phialpha.cpp
    //-------------------

    // E.Wu 2024-12-24
    // This file contains 3 subroutines:
    // 1. allocate_phialpha, which allocates memory for phialpha
    // 2. build_phialpha, which calculates the overlap
    // between atomic basis and projector alpha : <phi_mu|alpha>
    // which will be used in calculating pdm, gdmx, H_V_delta, F_delta;
    // 3. check_phialpha, which prints the results into .dat files
    // for checking

  public:
    // calculates <chi|alpha>
    void allocate_phialpha(const bool& cal_deri,
                           const UnitCell& ucell,
                           const LCAO_Orbitals& orb,
                           const Grid_Driver& GridD);

    void build_phialpha(const bool& cal_deri /**< [in] 0 for 2-center intergration, 1 for its derivation*/,
                        const UnitCell& ucell,
                        const LCAO_Orbitals& orb,
                        const Grid_Driver& GridD,
                        const TwoCenterIntegrator& overlap_orb_alpha);

    void check_phialpha(const bool& cal_deri /**< [in] 0 for 2-center intergration, 1 for its derivation*/,
                        const UnitCell& ucell,
                        const LCAO_Orbitals& orb,
                        const Grid_Driver& GridD);

    //-------------------
    // LCAO_deepks_pdm.cpp
    //-------------------

    // This file contains subroutines for calculating pdm,
    // which is defind as sum_mu,nu rho_mu,nu (<chi_mu|alpha><alpha|chi_nu>);
    // as well as gdmx, which is the gradient of pdm, defined as
    // sum_mu,nu rho_mu,nu d/dX(<chi_mu|alpha><alpha|chi_nu>)

    // It also contains subroutines for printing pdm and gdmx
    // for checking purpose

    // There are 4 subroutines in this file:
    // 1. cal_projected_DM, which is used for calculating pdm
    // 2. check_projected_dm, which prints pdm to descriptor.dat

    // 3. cal_gdmx, calculating gdmx (and optionally gdmepsl for stress)
    // 4. check_gdmx, which prints gdmx to a series of .dat files

  public:
    /**
     * @brief calculate projected density matrix:
     * pdm = sum_i,occ <phi_i|alpha1><alpha2|phi_k>
     * 3 cases to skip calculation of pdm:
     *    1. NSCF calculation of DeePKS, init_chg = file and pdm has been read
     *    2. SCF calculation of DeePKS with init_chg = file and pdm has been read for restarting SCF
     *    3. Relax/Cell-Relax/MD calculation, non-first step will use the convergence pdm from the last step as initial
     * pdm
     */
    template <typename TK>
    void cal_projected_DM(const elecstate::DensityMatrix<TK, double>* dm,
                          const UnitCell& ucell,
                          const LCAO_Orbitals& orb,
                          const Grid_Driver& GridD);

    void check_projected_dm();

    // calculate the gradient of pdm with regard to atomic positions
    // d/dX D_{Inl,mm'}
    template <typename TK>
    void cal_gdmx( // const ModuleBase::matrix& dm,
        const std::vector<std::vector<TK>>& dm,
        const UnitCell& ucell,
        const LCAO_Orbitals& orb,
        const Grid_Driver& GridD,
        const int nks,
        const std::vector<ModuleBase::Vector3<double>>& kvec_d,
        std::vector<hamilt::HContainer<double>*> phialpha,
        torch::Tensor& gdmx);

    void check_gdmx(const int nat, const torch::Tensor& gdmx);

    template <typename TK>
    void cal_gdmepsl( // const ModuleBase::matrix& dm,
        const std::vector<std::vector<TK>>& dm,
        const UnitCell& ucell,
        const LCAO_Orbitals& orb,
        const Grid_Driver& GridD,
        const int nks,
        const std::vector<ModuleBase::Vector3<double>>& kvec_d,
        std::vector<hamilt::HContainer<double>*> phialpha,
        torch::Tensor& gdmepsl);

    void check_gdmepsl(const torch::Tensor& gdmepsl);

    /**
     * @brief set init_pdm to skip the calculation of pdm in SCF iteration
     */
    void set_init_pdm(bool ipdm)
    {
        this->init_pdm = ipdm;
    }
    /**
     * @brief read pdm from file, do it only once in whole calculation
     */
    void read_projected_DM(bool read_pdm_file, bool is_equiv, const Numerical_Orbital& alpha);

    //-------------------
    // LCAO_deepks_vdelta.cpp
    //-------------------

    // This file contains subroutines related to V_delta, which is the deepks contribution to Hamiltonian
    // defined as |alpha>V(D)<alpha|
    // as well as subroutines for printing them for checking
    // It also contains subroutine related to calculating e_delta_bands, which is basically
    // tr (rho * V_delta)

    // Four subroutines are contained in the file:
    // 5. cal_e_delta_band : calculates e_delta_bands

  public:
    /// calculate tr(\rho V_delta)
    // void cal_e_delta_band(const std::vector<ModuleBase::matrix>& dm/**<[in] density matrix*/);
    template <typename TK>
    void cal_e_delta_band(const std::vector<std::vector<TK>>& dm /**<[in] density matrix*/, const int nks);

    //! a temporary interface for cal_e_delta_band and cal_e_delta_band_k
    template <typename TK>
    void dpks_cal_e_delta_band(const std::vector<std::vector<TK>>& dm, const int nks);

  public:
    //-------------------
    // LCAO_deepks_torch.cpp
    //-------------------

    // This file contains interfaces with libtorch,
    // including loading of model and calculating gradients
    // as well as subroutines that prints the results for checking

    // The file contains 8 subroutines:
    // 3. cal_gvx : gvx is used for training with force label, which is gradient of descriptors,
    //       calculated by d(des)/dX = d(pdm)/dX * d(des)/d(pdm) = gdmx * gvdm
    //       using einsum
    // 4. check_gvx : prints gvx into gvx.dat for checking
    // 5. cal_gvepsl : gvepsl is used for training with stress label, which is derivative of
    //       descriptors wrt strain tensor, calculated by
    //       d(des)/d\epsilon_{ab} = d(pdm)/d\epsilon_{ab} * d(des)/d(pdm) = gdmepsl * gvdm
    //       using einsum
    // 6. cal_gevdm : d(des)/d(pdm)
    //       calculated using torch::autograd::grad
    // 7. load_model : loads model for applying V_delta
    // 8. cal_gedm : calculates d(E_delta)/d(pdm)
    //       this is the term V(D) that enters the expression H_V_delta = |alpha>V(D)<alpha|
    //       caculated using torch::autograd::grad
    // 9. check_gedm : prints gedm for checking

  public:
    /// calculates gradient of descriptors w.r.t atomic positions
    ///----------------------------------------------------
    /// m, n: 2*l+1
    /// v: eigenvalues of dm , 2*l+1
    /// a,b: natom
    ///  - (a: the center of descriptor orbitals
    ///  - b: the atoms whose force being calculated)
    /// gvdm*gdmx->gvx
    ///----------------------------------------------------
    void cal_gvx(const int nat, const std::vector<torch::Tensor>& gevdm, const torch::Tensor& gdmx, torch::Tensor& gvx);
    void check_gvx(const int nat, const torch::Tensor& gvx);

    // for stress
    void cal_gvepsl(const int nat,
                    const std::vector<torch::Tensor>& gevdm,
                    const torch::Tensor& gdmepsl,
                    torch::Tensor& gvepsl);

    // load the trained neural network model
    void load_model(const std::string& model_file);

    /// calculate partial of energy correction to descriptors
    void cal_gedm(const int nat, const std::vector<torch::Tensor>& descriptor);
    void check_gedm();
    void cal_gedm_equiv(const int nat, const std::vector<torch::Tensor>& descriptor);

    // calculate gevdm
    void cal_gevdm(const int nat, std::vector<torch::Tensor>& gevdm);

  private:
    const Parallel_Orbitals* pv;
};

namespace GlobalC
{
extern LCAO_Deepks ld;
}

#endif
#endif
