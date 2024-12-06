#include "hsolver_pw_sdft.h"

#include "module_base/global_function.h"
#include "module_base/parallel_device.h"
#include "module_base/timer.h"
#include "module_base/tool_title.h"
#include "module_elecstate/module_charge/symmetry_rho.h"

#include <algorithm>

namespace hsolver
{
template <typename T, typename Device>
void HSolverPW_SDFT<T, Device>::solve(const UnitCell& ucell,
                                      hamilt::Hamilt<T, Device>* pHamilt,
                                      psi::Psi<T, Device>& psi,
                                      psi::Psi<T>& psi_cpu,
                                      elecstate::ElecState* pes,
                                      ModulePW::PW_Basis_K* wfc_basis,
                                      Stochastic_WF<T, Device>& stowf,
                                      const int istep,
                                      const int iter,
                                      const bool skip_charge)
{
    ModuleBase::TITLE("HSolverPW_SDFT", "solve");
    ModuleBase::timer::tick("HSolverPW_SDFT", "solve");

    const int npwx = psi.get_nbasis();
    const int nbands = psi.get_nbands();
    const int nks = psi.get_nk();

    //---------------------------------------------------------------------------------------------------------------
    //---------------------------------for psi init guess!!!!--------------------------------------------------------
    //---------------------------------------------------------------------------------------------------------------
    // if (!PARAM.inp.psi_initializer && !this->initialed_psi && this->basis_type == "pw")
    // {
    //     for (int ik = 0; ik < nks; ++ik)
    //     {
    //         /// update H(k) for each k point
    //         pHamilt->updateHk(ik);

    //         if (nbands > 0 && GlobalV::MY_STOGROUP == 0)
    //         {
    //             /// update psi pointer for each k point
    //             psi.fix_k(ik);

    //             /// for psi init guess!!!!
    //             hamilt::diago_PAO_in_pw_k2(this->ctx, ik, psi, this->wfc_basis, this->pwf, pHamilt);
    //         }
    //     }
    // }
    //---------------------------------------------------------------------------------------------------------------
    //---------------------------------------------------------------------------------------------------------------
    //---------------------------------------------------------------------------------------------------------------

    // prepare for the precondition of diagonalization
    std::vector<double> precondition(psi.get_nbasis(), 0.0);

    // report if the specified diagonalization method is not supported
    const std::initializer_list<std::string> _methods = {"cg", "dav", "dav_subspace", "bpcg"};
    if (std::find(std::begin(_methods), std::end(_methods), this->method) == std::end(_methods))
    {
        ModuleBase::WARNING_QUIT("HSolverPW::solve", "This type of eigensolver is not supported!");
    }

    // part of KSDFT to get KS orbitals
    for (int ik = 0; ik < nks; ++ik)
    {
        ModuleBase::timer::tick("HSolverPW_SDFT", "solve_KS");
        pHamilt->updateHk(ik);
        if (nbands > 0 && GlobalV::MY_STOGROUP == 0)
        {
            /// update psi pointer for each k point
            psi.fix_k(ik);
            /// template add precondition calculating here
            this->update_precondition(precondition, ik, this->wfc_basis->npwk[ik], pes->pot->get_vl_of_0());
            /// solve eigenvector and eigenvalue for H(k)
            double* p_eigenvalues = &(pes->ekb(ik, 0));
            this->hamiltSolvePsiK(pHamilt, psi, precondition, p_eigenvalues);
        }

#ifdef __MPI
        if (nbands > 0 && PARAM.inp.bndpar > 1)
        {
            Parallel_Common::bcast_dev(this->ctx, &psi(ik, 0, 0), npwx * nbands, PARAPW_WORLD, &psi_cpu(ik, 0, 0));
            MPI_Bcast(&pes->ekb(ik, 0), nbands, MPI_DOUBLE, 0, PARAPW_WORLD);
        }
#endif
        ModuleBase::timer::tick("HSolverPW_SDFT", "solve_KS");
        stoiter.orthog(ik, psi, stowf);
        stoiter.checkemm(ik, istep, iter, stowf); // check and reset emax & emin
    }

    this->output_iterInfo();

    for (int ik = 0; ik < nks; ik++)
    {
        // init k
        if (nks > 1)
        {
            pHamilt->updateHk(ik); // necessary , because emax and emin should be decided first
        }
        stoiter.calPn(ik, stowf);
    }

    // iterate to get mu
    stoiter.itermu(iter, pes);

    // prepare sqrt{f(\hat{H})}|\chi> to calculate density, force and stress
    stoiter.calHsqrtchi(stowf);

    // calculate eband = \sum_{ik,ib} w(ik)f(ik,ib)e_{ikib}, demet = -TS
    elecstate::ElecStatePW<T, Device>* pes_pw = static_cast<elecstate::ElecStatePW<T, Device>*>(pes);
    if (GlobalV::MY_STOGROUP == 0)
    {
        pes_pw->calEBand();
    }
    if (nbands > 0)
    {
#ifdef __MPI
        pes->f_en.eband /= GlobalV::NPROC_IN_POOL;
        MPI_Allreduce(MPI_IN_PLACE, &pes->f_en.eband, 1, MPI_DOUBLE, MPI_SUM, STO_WORLD);
        MPI_Bcast(&pes->f_en.eband, 1, MPI_DOUBLE, 0, PARAPW_WORLD);
#endif
    }
    stoiter.sum_stoeband(stowf, pes_pw, pHamilt, wfc_basis);
    
    

    // for nscf, skip charge
    if (skip_charge)
    {
        ModuleBase::timer::tick("HSolverPW_SDFT", "solve");
        return;
    }

    //(5) calculate new charge density
    // calculate KS rho.
    pes_pw->init_rho_data();
    if (nbands > 0)
    {
        pes_pw->psiToRho(psi);
    }
    // calculate stochastic rho
    stoiter.cal_storho(ucell, stowf, pes_pw,wfc_basis);

    // will do rho symmetry and energy calculation in esolver
    ModuleBase::timer::tick("HSolverPW_SDFT", "solve");
    return;
}

// template class HSolverPW_SDFT<std::complex<float>, base_device::DEVICE_CPU>;
template class HSolverPW_SDFT<std::complex<double>, base_device::DEVICE_CPU>;
#if ((defined __CUDA) || (defined __ROCM))
// template class HSolverPW_SDFT<std::complex<float>, base_device::DEVICE_GPU>;
template class HSolverPW_SDFT<std::complex<double>, base_device::DEVICE_GPU>;
#endif
} // namespace hsolver