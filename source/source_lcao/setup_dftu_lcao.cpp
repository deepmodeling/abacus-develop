#include "setup_dftu_lcao.h"
#include "source_lcao/module_dftu/dftu_lcao.h"
#include "source_lcao/module_dftu/dftu_occup.h"
#include "source_pw/module_pwdft/dftu_output.h" // mohan add 2025-11-08
#include "source_estate/module_dm/density_matrix.h"
#include "source_lcao/hamilt_lcao.h"

namespace ModuleESolver
{

template <typename TK>
void init_dftu_lcao(const int istep,
                     const int iter,
                     int dft_plus_u,
                     void* dftu,
                     void* dm,
                     const UnitCell& ucell,
                     double** rho,
                     const int nrxx)
{
    if (!dft_plus_u)
    {
        return;
    }
    
    auto* dftu_ptr = static_cast<Plus_U*>(dftu);
    auto* dm_ptr = static_cast<elecstate::DensityMatrix<TK, double>*>(dm);
    
    if (istep != 0 || iter != 1)
    {
        dftu_ptr->set_dmr(dm_ptr);
    }
    
    /// Calculate U and J if Yukawa potential is used
    dftu_ptr->cal_slater_UJ(ucell, rho, nrxx);
}

template <typename TK>
void finish_dftu_lcao(const int iter,
                       const bool conv_esolver,
                       int dft_plus_u,
                       bool out_chg,
                       void* dftu,
                       const UnitCell& ucell,
                       const std::vector<std::vector<TK>>& dm_vec,
                       const K_Vectors& kv,
                       const double mixing_beta,
                       void* hamilt_lcao,
                       const std::string& global_out_dir,
                       int nspin,
                       int npol)
{
    if (!dft_plus_u)
    {
        return;
    }
    
    auto* dftu_ptr = static_cast<Plus_U*>(dftu);
    auto* hamilt_lcao_ptr = static_cast<hamilt::HamiltLCAO<TK, double>*>(hamilt_lcao);
    
    /// old DFT+U method calculates energy correction in esolver,
    /// new DFT+U method calculates energy in Hamiltonian
    if (dft_plus_u == 2)
    {
        if (dftu_ptr->get_occ_mat_ctrl() != 2)
        {
            DFTU_LCAO::cal_occ_mat(iter, ucell, dm_vec, kv, mixing_beta,
                                   static_cast<hamilt::Hamilt<TK>*>(hamilt_lcao_ptr), *dftu_ptr);
        }
        dftu_ptr->cal_energy_correction(ucell, iter);
    }
    dftu_io::output(*dftu_ptr, ucell, out_chg, global_out_dir, nspin, npol);
    
    /// use the converged occupation matrix for next MD/Relax SCF calculation
    if (conv_esolver)
    {
        dftu_ptr->mark_occ_mat_initialized();
    }
}

/// Template instantiation
template void init_dftu_lcao<double>(const int istep,
                                      const int iter,
                                      int dft_plus_u,
                                      void* dftu,
                                      void* dm,
                                      const UnitCell& ucell,
                                      double** rho,
                                      const int nrxx);

template void init_dftu_lcao<std::complex<double>>(const int istep,
                                                    const int iter,
                                                    int dft_plus_u,
                                                    void* dftu,
                                                    void* dm,
                                                    const UnitCell& ucell,
                                                    double** rho,
                                                    const int nrxx);

template void finish_dftu_lcao<double>(const int iter,
                                        const bool conv_esolver,
                                        int dft_plus_u,
                                        bool out_chg,
                                        void* dftu,
                                        const UnitCell& ucell,
                                        const std::vector<std::vector<double>>& dm_vec,
                                        const K_Vectors& kv,
                                        const double mixing_beta,
                                        void* hamilt_lcao,
                                        const std::string& global_out_dir,
                                        int nspin,
                                        int npol);

template void finish_dftu_lcao<std::complex<double>>(const int iter,
                                                      const bool conv_esolver,
                                                      int dft_plus_u,
                                                      bool out_chg,
                                                      void* dftu,
                                                      const UnitCell& ucell,
                                                      const std::vector<std::vector<std::complex<double>>>& dm_vec,
                                                      const K_Vectors& kv,
                                                      const double mixing_beta,
                                                      void* hamilt_lcao,
                                                      const std::string& global_out_dir,
                                                      int nspin,
                                                      int npol);

} // namespace ModuleESolver
