#include "hsolver_lcaopw.h"

#include "source_base/global_variable.h"
#include "source_base/parallel_global.h" // for MPI
#include "source_base/timer.h"
#include "source_base/tool_quit.h"
#include "source_estate/elecstate_pw.h"
#include "source_pw/module_pwdft/global.h"
#include "source_pw/module_pwdft/hamilt_pw.h"
#include "source_hsolver/diago_iter_assist.h"
#include "source_io/module_parameter/parameter.h"
#include "source_estate/elecstate_tools.h"


#ifdef __EXX
#include "source_pw/module_pwdft/hamilt_lcaopw.h"
#endif

namespace hsolver
{

/*
    lcao_in_pw
*/
template <typename T>
void HSolverLIP<T>::solve(hamilt::Hamilt<T>* pHamilt, // ESolver_KS_PW::p_hamilt
                          psi::Psi<T>& psi,           // ESolver_KS_PW::kspw_psi
                          elecstate::ElecState* pes,  // ESolver_KS_PW::pes
                          psi::Psi<T>& transform,
                          const bool skip_charge,
                          const double tpiba,
                          const int nat)
{
    ModuleBase::TITLE("HSolverLIP", "solve");
    ModuleBase::timer::tick("HSolverLIP", "solve");
    std::vector<Real> eigenvalues(pes->ekb.nr * pes->ekb.nc, 0);
    for (int ik = 0; ik < this->wfc_basis->nks; ++ik)
    {
        /// update H(k) for each k point
        pHamilt->updateHk(ik);

        psi.fix_k(ik);
        transform.fix_k(ik);

#ifdef __EXX
        auto& exx_lip = dynamic_cast<hamilt::HamiltLIP<T>*>(pHamilt)->exx_lip;
        auto add_exx_to_subspace_hamilt = [&ik, &exx_lip](T* hcc, const int naos) -> void {
            if (GlobalC::exx_info.info_global.cal_exx)
            {
                for (int n = 0; n < naos; ++n)
                {
                    for (int m = 0; m < naos; ++m)
                    {
                        hcc[n * naos + m]
                            += (T)GlobalC::exx_info.info_global.hybrid_alpha * exx_lip.get_exx_matrix()[ik][m][n];
                    }
                }
            }
        };
        auto set_exxlip_lcaowfc = [&ik, &exx_lip](const T* const vcc, const int naos, const int nbands) -> void {
            if (GlobalC::exx_info.info_global.cal_exx)
            {
                exx_lip.set_hvec(ik, vcc, naos, nbands);
            }
        };
#endif
        /// solve eigenvector and eigenvalue for H(k)
        hsolver::DiagoIterAssist<T>::diagH_subspace_init(
            pHamilt,                 // interface to hamilt
            transform.get_pointer(), // transform matrix between lcao and pw
            transform.get_nbands(),
            transform.get_nbasis(),
            psi,                                  // psi in pw basis
            eigenvalues.data() + ik * pes->ekb.nc // eigenvalues
#ifdef __EXX
            ,
            add_exx_to_subspace_hamilt,
            set_exxlip_lcaowfc
#endif
        );

        if (skip_charge)
        {
            GlobalV::ofs_running << "Average iterative diagonalization steps for k-points " << ik
                                 << " is: " << DiagoIterAssist<T>::avg_iter
                                 << " ; where current threshold is: " << DiagoIterAssist<T>::PW_DIAG_THR << " . "
                                 << std::endl;
            DiagoIterAssist<T>::avg_iter = 0.0;
        }
        /// calculate the contribution of Psi for charge density rho
    }
    base_device::memory::cast_memory_op<double, Real, base_device::DEVICE_CPU, base_device::DEVICE_CPU>()(
        pes->ekb.c,
        eigenvalues.data(),
        pes->ekb.nr * pes->ekb.nc);

    elecstate::calculate_weights(pes->ekb,
                                 pes->wg,
                                 pes->klist,
                                 pes->eferm,
                                 pes->f_en,
                                 pes->nelec_spin,
                                 pes->skip_weights);
    elecstate::calEBand(pes->ekb,pes->wg,pes->f_en);
    if (skip_charge)
    {
        if (PARAM.globalv.use_uspp)
        {
            reinterpret_cast<elecstate::ElecStatePW<T>*>(pes)->cal_becsum(psi);
        }
        ModuleBase::timer::tick("HSolverLIP", "solve");
        return;
    }
    reinterpret_cast<elecstate::ElecStatePW<T>*>(pes)->psiToRho(psi);

    ModuleBase::timer::tick("HSolverLIP", "solve");
    return;
}

template class HSolverLIP<std::complex<float>>;
template class HSolverLIP<std::complex<double>>;

} // namespace hsolver
