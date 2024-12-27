#ifndef LCAO_DEEPKS_H
#define LCAO_DEEPKS_H

#ifdef __DEEPKS

#include "deepks_force.h"
#include "deepks_hmat.h"
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

    ///(Unit: Ry)  \f$tr(\rho_{HL} H_\delta),
    ///\rho_{HL} = c_{L, \mu}c_{L,\nu} - c_{H, \mu}c_{H,\nu} \f$ (for gamma_only)
    ModuleBase::matrix o_delta;

    /// Correction term to the Hamiltonian matrix: \f$\langle\phi|V_\delta|\phi\rangle\f$ (for gamma only)
    /// The size of first dimension is 1, which is used for the consitence with H_V_delta_k
    std::vector<std::vector<double>> H_V_delta;
    /// Correction term to Hamiltonian, for multi-k
    std::vector<std::vector<std::complex<double>>> H_V_delta_k;

    // F_delta will be deleted soon, mohan 2024-07-25
    ///(Unit: Ry/Bohr) Total Force due to the DeePKS correction term \f$E_{\delta}\f$
    ModuleBase::matrix F_delta;

    // k index of HOMO for multi-k bandgap label. QO added 2022-01-24
    int h_ind = 0;

    // k index of LUMO for multi-k bandgap label. QO added 2022-01-24
    int l_ind = 0;

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
    double** pdm; //[tot_Inl][2l+1][2l+1]	caoyu modified 2021-05-07; if equivariant version: [nat][nlm*nlm]
    std::vector<torch::Tensor> pdm_tensor;

    // descriptors
    std::vector<torch::Tensor> d_tensor;

    // gedm:dE/dD, [tot_Inl][2l+1][2l+1]	(E: Hartree)
    std::vector<torch::Tensor> gedm_tensor;

    // gdmx: dD/dX		\sum_{mu,nu} 2*c_mu*c_nu * <dphi_mu/dx|alpha_m><alpha_m'|phi_nu>
    double*** gdmx; //[natom][tot_Inl][2l+1][2l+1]
    double*** gdmy;
    double*** gdmz;

    // gdm_epsl: dD/d\epsilon_{\alpha\beta}
    double*** gdm_epsl; //[6][tot_Inl][2l+1][2l+1]

    // dD/d\epsilon_{\alpha\beta}, tensor form of gdm_epsl
    std::vector<torch::Tensor> gdmepsl_vector;

    // gv_epsl:d(d)/d\epsilon_{\alpha\beta}, [natom][6][des_per_atom]
    torch::Tensor gvepsl_tensor;

    /// dE/dD, autograd from loaded model(E: Ry)
    double** gedm; //[tot_Inl][2l+1][2l+1]

    // gvx:d(d)/dX, [natom][3][natom][des_per_atom]
    torch::Tensor gvx_tensor;

    // d(d)/dD, autograd from torch::linalg::eigh
    std::vector<torch::Tensor> gevdm_vector;

    // dD/dX, tensor form of gdmx
    std::vector<torch::Tensor> gdmr_vector;

    // orbital_pdm_shell:[1,Inl,nm*nm]; \langle \phi_\mu|\alpha\rangle\langle\alpha|\phi_\nu\ranlge
    double**** orbital_pdm_shell;
    // orbital_precalc:[1,NAt,NDscrpt]; gvdm*orbital_pdm_shell
    torch::Tensor orbital_precalc_tensor;

    // v_delta_pdm_shell[nks,nlocal,nlocal,Inl,nm*nm] = overlap * overlap
    double***** v_delta_pdm_shell;
    std::complex<double>***** v_delta_pdm_shell_complex; // for multi-k
    // v_delta_precalc[nks,nlocal,nlocal,NAt,NDscrpt] = gvdm * v_delta_pdm_shell;
    torch::Tensor v_delta_precalc_tensor;
    // for v_delta==2 , new v_delta_precalc storage method
    torch::Tensor phialpha_tensor;
    torch::Tensor gevdm_tensor;

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
    // 2. subroutines that are related to calculating force label:
    //   - init_gdmx : allocates gdmx; it is a private subroutine
    //   - del_gdmx : releases gdmx
    // 3. subroutines that are related to V_delta:
    //   - allocate_V_delta : allocates H_V_delta; if calculating force, it also calls
    //       init_gdmx, as well as allocating F_delta

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

    // array for storing gdmx, used for calculating gvx
    void init_gdmx(const int nat);
    // void del_gdmx(const int nat);
    void del_gdmx();

    // array for storing gdm_epsl, used for calculating gvx
    void init_gdmepsl();
    void del_gdmepsl();

  private:
    // arrange index of descriptor in all atoms
    void init_index(const int ntype, const int nat, std::vector<int> na, const int tot_inl, const LCAO_Orbitals& orb);
    // data structure that saves <phi|alpha>
    void allocate_nlm(const int nat);

    // for bandgap label calculation; QO added on 2022-1-7
    void init_orbital_pdm_shell(const int nks);
    void del_orbital_pdm_shell(const int nks);

    // for v_delta label calculation; xinyuan added on 2023-2-22
    void init_v_delta_pdm_shell(const int nks, const int nlocal);
    void del_v_delta_pdm_shell(const int nks, const int nlocal);

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

    // 3. cal_gdmx, calculating gdmx (and optionally gdm_epsl for stress)
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
        const bool isstress);

    void check_gdmx(const int nat);

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

    //-------------------
    // LCAO_deepks_odelta.cpp
    //-------------------

    // This file contains subroutines for calculating O_delta,
    // which corresponds to the correction of the band gap.

  public:
    template <typename TK, typename TH>
    void cal_o_delta(
        const std::vector<std::vector<TH>>& dm_hl /**<[in] modified density matrix that contains HOMO and LUMO only*/,
        const int nks);

    //-------------------
    // LCAO_deepks_torch.cpp
    //-------------------

    // This file contains interfaces with libtorch,
    // including loading of model and calculating gradients
    // as well as subroutines that prints the results for checking

    // The file contains 8 subroutines:
    // 1. cal_descriptor : obtains descriptors which are eigenvalues of pdm
    //       by calling torch::linalg::eigh
    // 2. check_descriptor : prints descriptor for checking
    // 3. cal_gvx : gvx is used for training with force label, which is gradient of descriptors,
    //       calculated by d(des)/dX = d(pdm)/dX * d(des)/d(pdm) = gdmx * gvdm
    //       using einsum
    // 4. check_gvx : prints gvx into gvx.dat for checking
    // 5. cal_gvepsl : gvepsl is used for training with stress label, which is derivative of
    //       descriptors wrt strain tensor, calculated by
    //       d(des)/d\epsilon_{ab} = d(pdm)/d\epsilon_{ab} * d(des)/d(pdm) = gdm_epsl * gvdm
    //       using einsum
    // 6. cal_gvdm : d(des)/d(pdm)
    //       calculated using torch::autograd::grad
    // 7. load_model : loads model for applying V_delta
    // 8. cal_gedm : calculates d(E_delta)/d(pdm)
    //       this is the term V(D) that enters the expression H_V_delta = |alpha>V(D)<alpha|
    //       caculated using torch::autograd::grad
    // 9. check_gedm : prints gedm for checking
    // 10. cal_orbital_precalc : orbital_precalc is usted for training with orbital label,
    //                          which equals gvdm * orbital_pdm_shell,
    //                          orbital_pdm_shell[1,Inl,nm*nm] = dm_hl * overlap * overlap
    // 11. cal_v_delta_precalc : v_delta_precalc is used for training with v_delta label,
    //                         which equals gvdm * v_delta_pdm_shell,
    //                         v_delta_pdm_shell = overlap * overlap
    // 12. check_v_delta_precalc : check v_delta_precalc
    // 13. prepare_phialpha : prepare phialpha for outputting npy file
    // 14. prepare_gevdm : prepare gevdm for outputting npy file

  public:
    /// Calculates descriptors
    /// which are eigenvalues of pdm in blocks of I_n_l
    void cal_descriptor(const int nat);
    /// print descriptors based on LCAO basis
    void check_descriptor(const UnitCell& ucell, const std::string& out_dir);

    void cal_descriptor_equiv(const int nat);

    /// calculates gradient of descriptors w.r.t atomic positions
    ///----------------------------------------------------
    /// m, n: 2*l+1
    /// v: eigenvalues of dm , 2*l+1
    /// a,b: natom
    ///  - (a: the center of descriptor orbitals
    ///  - b: the atoms whose force being calculated)
    /// gvdm*gdmx->gvx
    ///----------------------------------------------------
    void cal_gvx(const int nat);
    void check_gvx(const int nat);

    // for stress
    void cal_gvepsl(const int nat);

    // load the trained neural network model
    void load_model(const std::string& model_file);

    /// calculate partial of energy correction to descriptors
    void cal_gedm(const int nat);
    void check_gedm();
    void cal_gedm_equiv(const int nat);

    // calculates orbital_precalc
    template <typename TK, typename TH>
    void cal_orbital_precalc(const std::vector<std::vector<TH>>& dm_hl /**<[in] density matrix*/,
                             const int nat,
                             const int nks,
                             const std::vector<ModuleBase::Vector3<double>>& kvec_d,
                             const UnitCell& ucell,
                             const LCAO_Orbitals& orb,
                             const Grid_Driver& GridD);

    // calculates v_delta_precalc
    template <typename TK>
    void cal_v_delta_precalc(const int nlocal,
                             const int nat,
                             const int nks,
                             const std::vector<ModuleBase::Vector3<double>>& kvec_d,
                             const UnitCell& ucell,
                             const LCAO_Orbitals& orb,
                             const Grid_Driver& GridD);

    template <typename TK>
    void check_v_delta_precalc(const int nat, const int nks, const int nlocal);

    // prepare phialpha for outputting npy file
    template <typename TK>
    void prepare_phialpha(const int nlocal,
                          const int nat,
                          const int nks,
                          const std::vector<ModuleBase::Vector3<double>>& kvec_d,
                          const UnitCell& ucell,
                          const LCAO_Orbitals& orb,
                          const Grid_Driver& GridD);

    template <typename TK>
    void check_vdp_phialpha(const int nat, const int nks, const int nlocal);

    // prepare gevdm for outputting npy file
    void prepare_gevdm(const int nat, const LCAO_Orbitals& orb);
    void check_vdp_gevdm(const int nat);

  private:
    const Parallel_Orbitals* pv;
    void cal_gvdm(const int nat);

#ifdef __MPI

  public:
    // reduces a dim 2 array
    void allsum_deepks(int inlmax,    // first dimension
                       int ndim,      // second dimension
                       double** mat); // the array being reduced
#endif
};

namespace GlobalC
{
extern LCAO_Deepks ld;
}

#endif
#endif
