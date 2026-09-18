#include "hsolver_lcaopw.h"

#include "source_base/parallel_global.h" // for MPI
#include "source_base/timer.h"
#include "source_base/tool_quit.h"
#include "source_estate/elecstate_pw.h"
#include "source_estate/elecstate_tools.h"
#include "source_hsolver/diag_comm_info.h"
#include "source_hsolver/diago_iter_assist.h"

#include <ostream>

namespace hsolver
{

/*
    lcao_in_pw
*/
template <typename T>
void HSolverLIP<T>::solve(HSOperator<T>& op,       // ESolver_KS_PW::p_hamilt behind the operator interface
                          psi::Psi<T>& psi,           // ESolver_KS_PW::kspw_psi
                          elecstate::ElecState* pes,  // ESolver_KS_PW::pes
                          psi::Psi<T>& transform,
                          const diag_comm_info& diag_comm,
                          std::ostream& log,
                          const bool skip_charge)
{
    ModuleBase::TITLE("HSolverLIP", "solve");
    ModuleBase::timer::start("HSolverLIP", "solve");
    std::vector<Real> eigenvalues(pes->ekb.nr * pes->ekb.nc, 0);
    for (int ik = 0; ik < this->wfc_basis->nks; ++ik)
    {
        /// update H(k) for each k point
        op.update_k(ik);

        psi.fix_k(ik);
        transform.fix_k(ik);

        /// solve eigenvector and eigenvalue for H(k)
        hsolver::DiagoIterAssist<T>::diag_subspace_init(op,
                                                        transform.get_pointer(), // transform matrix between lcao and pw
                                                        transform.get_nbands(),
                                                        transform.get_nbasis(),
                                                        psi,                                   // psi in pw basis
                                                        eigenvalues.data() + ik * pes->ekb.nc, // eigenvalues
                                                        this->basis_type,
                                                        this->calculation,
                                                        diag_comm);

        if (skip_charge)
        {
            log << "Average iterative diagonalization steps for k-points " << ik
                << " is: " << DiagoIterAssist<T>::avg_iter
                << " ; where current threshold is: " << DiagoIterAssist<T>::PW_DIAG_THR << " . " << std::endl;
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
                                 this->global_nbands,
                                 pes->skip_weights);
    elecstate::calEBand(pes->ekb,pes->wg,pes->f_en);
    if (skip_charge)
    {
        if (this->use_uspp)
        {
            reinterpret_cast<elecstate::ElecStatePW<T>*>(pes)->cal_becsum(psi);
        }
        ModuleBase::timer::end("HSolverLIP", "solve");
        return;
    }
    reinterpret_cast<elecstate::ElecStatePW<T>*>(pes)->psiToRho(psi);

    ModuleBase::timer::end("HSolverLIP", "solve");
    return;
}

template class HSolverLIP<std::complex<float>>;
template class HSolverLIP<std::complex<double>>;

} // namespace hsolver
