#include "dftu_lcao.h"
#include "source_lcao/module_dftu/dftu.h"
#include "source_estate/module_dm/density_matrix.h"
#include "source_lcao/hamilt_lcao.h"

namespace ModuleESolver
{

template <typename TK>
void init_dftu_lcao(const int istep,
                     const int iter,
                     const Input_para& inp,
                     void* dftu,
                     void* dm,
                     const UnitCell& ucell,
                     double** rho,
                     const int nrxx)
{
    if (!inp.dft_plus_u)
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
                       const Input_para& inp,
                       void* dftu,
                       const UnitCell& ucell,
                       const std::vector<std::vector<TK>>& dm_vec,
                       const K_Vectors& kv,
                       const double mixing_beta,
                       void* hamilt_lcao)
{
    if (!inp.dft_plus_u)
    {
        return;
    }
    
    auto* dftu_ptr = static_cast<Plus_U*>(dftu);
    auto* hamilt_lcao_ptr = static_cast<hamilt::HamiltLCAO<TK, double>*>(hamilt_lcao);
    
    /// old DFT+U method calculates energy correction in esolver,
    /// new DFT+U method calculates energy in Hamiltonian
    if (inp.dft_plus_u == 2)
    {
        if (dftu_ptr->omc != 2)
        {
            dftu_cal_occup_m(iter, ucell, dm_vec, kv, mixing_beta, 
                             static_cast<hamilt::Hamilt<TK>*>(hamilt_lcao_ptr), *dftu_ptr);
        }
        dftu_ptr->cal_energy_correction(ucell, iter);
    }
    dftu_ptr->output(ucell);
    
    /// use the converged occupation matrix for next MD/Relax SCF calculation
    if (conv_esolver)
    {
        dftu_ptr->initialed_locale = true;
    }
}

/// Template instantiation
template void init_dftu_lcao<double>(const int istep,
                                      const int iter,
                                      const Input_para& inp,
                                      void* dftu,
                                      void* dm,
                                      const UnitCell& ucell,
                                      double** rho,
                                      const int nrxx);
template void init_dftu_lcao<std::complex<double>>(const int istep,
                                                    const int iter,
                                                    const Input_para& inp,
                                                    void* dftu,
                                                    void* dm,
                                                    const UnitCell& ucell,
                                                    double** rho,
                                                    const int nrxx);

template void finish_dftu_lcao<double>(const int iter,
                                        const bool conv_esolver,
                                        const Input_para& inp,
                                        void* dftu,
                                        const UnitCell& ucell,
                                        const std::vector<std::vector<double>>& dm_vec,
                                        const K_Vectors& kv,
                                        const double mixing_beta,
                                        void* hamilt_lcao);
template void finish_dftu_lcao<std::complex<double>>(const int iter,
                                                      const bool conv_esolver,
                                                      const Input_para& inp,
                                                      void* dftu,
                                                      const UnitCell& ucell,
                                                      const std::vector<std::vector<std::complex<double>>>& dm_vec,
                                                      const K_Vectors& kv,
                                                      const double mixing_beta,
                                                      void* hamilt_lcao);

} // namespace ModuleESolver
