#include "esolver_ks_lcaopw.h"

#include "module_hamilt_pw/hamilt_pwdft/elecond.h"
#include "module_io/input_conv.h"
#include "module_io/nscf_band.h"
#include "module_io/output_log.h"
#include "module_io/write_dos_pw.h"
#include "module_io/write_istate_info.h"
#include "module_io/write_wfc_pw.h"

#include <iostream>

//--------------temporary----------------------------
#include "module_elecstate/module_charge/symmetry_rho.h"
#include "module_elecstate/occupy.h"
#include "module_hamilt_general/module_ewald/H_Ewald_pw.h"
#include "module_hamilt_pw/hamilt_pwdft/global.h"
#include "module_io/print_info.h"
//-----force-------------------
#include "module_hamilt_pw/hamilt_pwdft/forces.h"
//-----stress------------------
#include "module_hamilt_pw/hamilt_pwdft/stress_pw.h"
//---------------------------------------------------
#include "module_base/memory.h"
#include "module_elecstate/elecstate_pw.h"
#include "module_hamilt_pw/hamilt_pwdft/hamilt_pw.h"
#include "module_hsolver/diago_iter_assist.h"
#include "module_hsolver/hsolver_lcaopw.h"
#include "module_hsolver/kernels/dngvd_op.h"
#include "module_hsolver/kernels/math_kernel_op.h"
#include "module_io/berryphase.h"
#include "module_io/numerical_basis.h"
#include "module_io/numerical_descriptor.h"
#include "module_io/rho_io.h"
#include "module_io/to_wannier90_pw.h"
#include "module_io/winput.h"
#include "module_io/write_pot.h"
#include "module_io/write_wfc_r.h"

#include <ATen/kernels/blas.h>
#include <ATen/kernels/lapack.h>

namespace ModuleESolver
{

    template <typename T>
    ESolver_KS_LIP<T>::ESolver_KS_LIP()
    {
        this->classname = "ESolver_KS_LIP";
        this->basisname = "LIP";
    }

    template <typename T>
    void ESolver_KS_LIP<T>::allocate_hsolver()
    {
        this->phsol = new hsolver::HSolverLIP<T>(this->pw_wfc);
    }
    template <typename T>
    void ESolver_KS_LIP<T>::deallocate_hsolver()
    {
        if (this->phsol != nullptr)
        {
            delete reinterpret_cast<hsolver::HSolverLIP<T>*>(this->phsol);
            this->phsol = nullptr;
        }
    }
    template <typename T>
    void ESolver_KS_LIP<T>::hamilt2density(const int istep, const int iter, const double ethr)
    {
        ModuleBase::TITLE("ESolver_KS_LIP", "hamilt2density");
        ModuleBase::timer::tick("ESolver_KS_LIP", "hamilt2density");

        if (this->phsol != nullptr)
        {
            // reset energy
            this->pelec->f_en.eband = 0.0;
            this->pelec->f_en.demet = 0.0;
            // choose if psi should be diag in subspace
            // be careful that istep start from 0 and iter start from 1
            // if (iter == 1)
            hsolver::DiagoIterAssist<T>::need_subspace = ((istep == 0 || istep == 1) && iter == 1) ? false : true;
            hsolver::DiagoIterAssist<T>::SCF_ITER = iter;
            hsolver::DiagoIterAssist<T>::PW_DIAG_THR = ethr;
            hsolver::DiagoIterAssist<T>::PW_DIAG_NMAX = GlobalV::PW_DIAG_NMAX;

            // It is not a good choice to overload another solve function here, this will spoil the concept of
            // multiple inheritance and polymorphism. But for now, we just do it in this way.
            // In the future, there will be a series of class ESolver_KS_LCAO_PW, HSolver_LCAO_PW and so on.
            std::weak_ptr<psi::Psi<T>> psig = this->p_wf_init->get_psig();

            if (psig.expired())
            {
                ModuleBase::WARNING_QUIT("ESolver_KS_PW::hamilt2density", "psig lifetime is expired");
            }

            // from HSolverLIP
            this->phsol->solve(this->p_hamilt,        // hamilt::Hamilt<T>* pHamilt,
                this->kspw_psi[0],     // psi::Psi<T>& psi,
                this->pelec,           // elecstate::ElecState<T>* pelec,
                psig.lock().get()[0]); // psi::Psi<T>& transform,

            if (GlobalV::out_bandgap)
            {
                if (!GlobalV::TWO_EFERMI)
                {
                    this->pelec->cal_bandgap();
                }
                else
                {
                    this->pelec->cal_bandgap_updw();
                }
            }
        }
        else
        {
            ModuleBase::WARNING_QUIT("ESolver_KS_LIP", "HSolver has not been allocated.");
        }
        // add exx
#ifdef __LCAO
#ifdef __EXX
        this->pelec->set_exx(GlobalC::exx_lip.get_exx_energy()); // Peize Lin add 2019-03-09
#endif
#endif

        // calculate the delta_harris energy
        // according to new charge density.
        // mohan add 2009-01-23
        this->pelec->cal_energies(1);

        Symmetry_rho srho;
        for (int is = 0; is < GlobalV::NSPIN; is++)
        {
            srho.begin(is, *(this->pelec->charge), this->pw_rhod, GlobalC::Pgrid, GlobalC::ucell.symm);
        }

        // compute magnetization, only for LSDA(spin==2)
        GlobalC::ucell.magnet.compute_magnetization(this->pelec->charge->nrxx,
            this->pelec->charge->nxyz,
            this->pelec->charge->rho,
            this->pelec->nelec_spin.data());

        // deband is calculated from "output" charge density calculated
        // in sum_band
        // need 'rho(out)' and 'vr (v_h(in) and v_xc(in))'
        this->pelec->f_en.deband = this->pelec->cal_delta_eband();

        ModuleBase::timer::tick("ESolver_KS_LIP", "hamilt2density");
    }

    template class ESolver_KS_LIP<std::complex<float>>;
    template class ESolver_KS_LIP<std::complex<double>>;
    // LIP is not supported on GPU yet.
} // namespace ModuleESolver
