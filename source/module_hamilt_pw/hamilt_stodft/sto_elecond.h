#ifndef STOELECOND_H
#define STOELECOND_H

#include "module_hamilt_general/hamilt.h"
#include "module_hamilt_pw/hamilt_pwdft/elecond.h"
#include "module_hamilt_pw/hamilt_stodft/sto_wf.h"
#include "source_hsolver/hsolver_pw_sdft.h"

template <typename FPTYPE, typename Device> 
class Sto_EleCond : protected EleCond<FPTYPE, Device>
{
  public:
#ifdef __ENABLE_FLOAT_FFTW
    using lowTYPE = float; // Here we use float to accelerate the calculation, which is enough for the accuracy
#else
    using lowTYPE = double;
#endif
    using lcomplex = std::complex<lowTYPE>;
    using resmem_lcomplex_op = base_device::memory::resize_memory_op<std::complex<lowTYPE>, Device>; 
    using delmem_lcomplex_op = base_device::memory::delete_memory_op<std::complex<lowTYPE>, Device>;
    using cpymem_lcomplex_op = base_device::memory::synchronize_memory_op<std::complex<lowTYPE>, Device, Device>;
    using castmem_lcomplex_op = base_device::memory::cast_memory_op<std::complex<lowTYPE>, std::complex<FPTYPE>, Device, Device>;
    using cpymem_complex_op = base_device::memory::synchronize_memory_op<std::complex<FPTYPE>, Device, Device>;
  public:
    Sto_EleCond(UnitCell* p_ucell_in,
                K_Vectors* p_kv_in,
                elecstate::ElecState* p_elec_in,
                ModulePW::PW_Basis_K* p_wfcpw_in,
                psi::Psi<std::complex<FPTYPE>, Device>* p_psi_in,
                pseudopot_cell_vnl* p_ppcell_in,
                hamilt::Hamilt<std::complex<FPTYPE>, Device>* p_hamilt_in,
                StoChe<FPTYPE, Device>& stoche,
                Stochastic_WF<std::complex<FPTYPE>, Device>* p_stowf_in);
    ~Sto_EleCond(){
        delete hamilt_sto_;
    };
    /**
     * @brief Set the N order of Chebyshev expansion for conductivities
     *        It will change class member : fd_nche, cond_nche
     *
     * @param dt t step
     * @param nbatch number of t batch
     * @param cond_thr threshold of errors for conductivities
     * @param fd_nche  N order of Chebyshev for Fermi-Dirac function
     * @param try_emin trial Emin
     * @param try_emax trial Emax
     *
     */
    void decide_nche(const FPTYPE dt, const FPTYPE cond_thr, const int& fd_nche, FPTYPE try_emin, FPTYPE try_emax);
    /**
     * @brief calculate Stochastic Kubo-Greenwood
     *
     * @param fwhmin FWHM
     * @param smear_type 1: Gaussian, 2: Lorentzian
     * @param wcut cutoff omega
     * @param dw_in omega step
     * @param dt_in t step
     * @param nonlocal whether to include the nonlocal potential corrections for velocity operator
     * @param nbatch t step batch
     * @param npart_sto number stochastic wavefunctions parts to evalution simultaneously
     */
    void sKG(const int& smear_type,
             const double& fwhmin,
             const double& wcut,
             const double& dw_in,
             const double& dt_in,
             const bool& nonlocal,
             const int& npart_sto);

  protected:
    int nbands_ks = 0;    ///< number of KS bands
    int nbands_sto = 0;   ///< number of stochastic bands
    int cond_nche = 0;    ///< number of Chebyshev orders for conductivities
    int fd_nche = 0;      ///< number of Chebyshev orders for Fermi-Dirac function
    int cond_dtbatch = 0; ///< number of time steps in a batch
    hamilt::Hamilt<std::complex<FPTYPE>, Device>* p_hamilt = nullptr; ///< pointer to the Hamiltonian
    Stochastic_WF<std::complex<FPTYPE>, Device>* p_stowf = nullptr;   ///< pointer to the stochastic wavefunctions
    Sto_Func<FPTYPE> stofunc;                                         ///< functions

    hamilt::HamiltSdftPW<std::complex<FPTYPE>, Device>* p_hamilt_sto = nullptr; ///< pointer to the Hamiltonian for sDFT
    hamilt::HamiltSdftPW<std::complex<lowTYPE>, Device>* hamilt_sto_ = nullptr; ///< pointer to the Hamiltonian for sDFT
    lowTYPE low_emin_ = 0;                                                      ///< Emin of the Hamiltonian for sDFT
    lowTYPE low_emax_ = 0;                                                      ///< Emax of the Hamiltonian for sDFT
  protected:
    /**
     * @brief calculate Jmatrix  <leftv|J|rightv>
     *
     */
    void cal_jmatrix(hamilt::HamiltSdftPW<std::complex<lowTYPE>, Device>* hamilt,
                     const psi::Psi<std::complex<lowTYPE>, Device>& kspsi_all,
                     const psi::Psi<std::complex<lowTYPE>, Device>& vkspsi,
                     const double* en,
                     const double* en_all,
                     std::complex<FPTYPE>* leftfact,
                     std::complex<FPTYPE>* rightfact,
                     psi::Psi<std::complex<lowTYPE>, Device>& leftchi,
                     psi::Psi<std::complex<lowTYPE>, Device>& rightchi,
                     psi::Psi<std::complex<lowTYPE>, Device>& left_hchi,
                     psi::Psi<std::complex<lowTYPE>, Device>& right_hchi,
                     psi::Psi<std::complex<lowTYPE>, Device>& batch_vchi,
                     psi::Psi<std::complex<lowTYPE>, Device>& batch_vhchi,
#ifdef __MPI
                     psi::Psi<std::complex<lowTYPE>, Device>& chi_all,
                     psi::Psi<std::complex<lowTYPE>, Device>& hchi_all,
                     void* gatherinfo_ks,
                     void* gatherinfo_sto,
#endif
                     const int& bsize_psi,
                     std::complex<lowTYPE>* j1,
                     std::complex<lowTYPE>* j2,
                     std::complex<lowTYPE>* tmpj,
                     hamilt::Velocity<lowTYPE, Device>& velop,
                     const int& ik,
                     const std::complex<lowTYPE>& factor,
                     const int bandinfo[6]);
};
#endif // ELECOND_H