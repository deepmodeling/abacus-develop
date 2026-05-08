#include "source_base/timer.h"
#include "source_base/tool_title.h"
#include "source_hsolver/diago_iter_assist.h"
#include "source_io/module_parameter/parameter.h"
#include "spin_constrain.h"
#include "source_pw/module_pwdft/onsite_proj.h"
#include "source_base/parallel_reduce.h"
#include "source_base/kernels/math_kernel_op.h"
#include "source_hsolver/hsolver_lcao.h"
#include "source_hsolver/hsolver_pw.h"
#include "source_estate/elecstate_pw.h"
#include "source_estate/elecstate_tools.h"

#ifdef __LCAO
#include "source_estate/elecstate_lcao.h"
#include "source_estate/elecstate_tools.h"
#include "source_estate/module_dm/cal_dm_psi.h"
#include "source_lcao/module_operator_lcao/dspin_lcao.h"
#endif

template <>
void spinconstrain::SpinConstrain<std::complex<double>>::calculate_delta_hcc(std::complex<double>* h_tmp, const std::complex<double>* becp_k, const ModuleBase::Vector3<double>* delta_lambda, const int nbands, const int nkb, const int* nh_iat, const int ik)
{
    ModuleBase::TITLE("spinconstrain::SpinConstrain", "calculate_delta_hcc");
    ModuleBase::timer::start("spinconstrain::SpinConstrain", "calculate_delta_hcc");
    
    int sum = 0;
    int size_ps = nkb * this->npol_ * nbands;
    std::complex<double>* becp_cpu = nullptr;
    if(PARAM.inp.device == "gpu")
    {
#if ((defined __CUDA) || (defined __ROCM))
        base_device::memory::resize_memory_op<std::complex<double>, base_device::DEVICE_CPU>()(becp_cpu, size_ps);
        base_device::memory::synchronize_memory_op<std::complex<double>, base_device::DEVICE_CPU, base_device::DEVICE_GPU>()(becp_cpu, becp_k, size_ps);   
#endif
    }
    else if (PARAM.inp.device == "cpu")
    {
        becp_cpu = const_cast<std::complex<double>*>(becp_k);
    }

    // Compute modified projector coefficients
    std::vector<std::complex<double>> ps(size_ps, 0.0);
    if(this->npol_ == 2)
    {
        // nspin=4: full Pauli matrix treatment
        for (int iat = 0; iat < this->Mi_.size(); iat++)
        {
            const int nproj = nh_iat[iat];
            const std::complex<double> coefficients0(delta_lambda[iat][2], 0.0);
            const std::complex<double> coefficients1(delta_lambda[iat][0] , delta_lambda[iat][1]);
            const std::complex<double> coefficients2(delta_lambda[iat][0] , -1 * delta_lambda[iat][1]);
            const std::complex<double> coefficients3(-1 * delta_lambda[iat][2], 0.0);
            for (int ib = 0; ib < nbands * this->npol_; ib += this->npol_)
            {
                for (int ip = 0; ip < nproj; ip++)
                {
                    const int becpind = ib * nkb + sum + ip;
                    const std::complex<double> becp1 = becp_cpu[becpind];
                    const std::complex<double> becp2 = becp_cpu[becpind + nkb];
                    ps[becpind] += coefficients0 * becp1
                                    + coefficients2 * becp2;
                    ps[becpind + nkb] += coefficients1 * becp1
                                        + coefficients3 * becp2;
                }
            }
            sum += nproj;
        }
    }
    else if(this->npol_ == 1)
    {
        // nspin=2: only z-component (spin collinear)
        for (int iat = 0; iat < this->Mi_.size(); iat++)
        {
            const int nproj = nh_iat[iat];
            double coefficients0 = delta_lambda[iat][2] * this->get_spin_sign(ik);
            for (int ib = 0; ib < nbands; ib++)
            {
                for (int ip = 0; ip < nproj; ip++)
                {
                    const int becpind = ib * nkb + sum + ip;
                    const std::complex<double> becp1 = becp_cpu[becpind];
                    ps[becpind] += coefficients0 * becp1;
                }
            }
            sum += nproj;
        }
    }
    
    std::complex<double>* ps_pointer = nullptr;
    if(PARAM.inp.device == "gpu")
    {
#if ((defined __CUDA) || (defined __ROCM))
        base_device::memory::resize_memory_op<std::complex<double>, base_device::DEVICE_GPU>()(ps_pointer, size_ps);
        base_device::memory::synchronize_memory_op<std::complex<double>, base_device::DEVICE_GPU, base_device::DEVICE_CPU>()(ps_pointer, ps.data(), size_ps);   
#endif
    }
    else if (PARAM.inp.device == "cpu")
    {
        ps_pointer = ps.data();
    }
    
    // update h_tmp by becp_k * ps: H += becp^† * ps
    char transa = 'C';
    char transb = 'N';
    const int npm = nkb * this->npol_;
    if (PARAM.inp.device == "gpu")
    {
#if ((defined __CUDA) || (defined __ROCM))
        ModuleBase::gemm_op<std::complex<double>, base_device::DEVICE_GPU>()(
            transa,
            transb,
            nbands,
            nbands,
            npm,
            &ModuleBase::ONE,
            becp_k,
            npm,
            ps_pointer,
            npm,
            &ModuleBase::ONE,
            h_tmp,
            nbands
        );
        base_device::memory::delete_memory_op<std::complex<double>, base_device::DEVICE_GPU>()(ps_pointer);
        base_device::memory::delete_memory_op<std::complex<double>, base_device::DEVICE_CPU>()(becp_cpu);
#endif

    }
    else if (PARAM.inp.device == "cpu")
    {
        ModuleBase::gemm_op<std::complex<double>, base_device::DEVICE_CPU>()(
            transa,
            transb,
            nbands,
            nbands,
            npm,
            &ModuleBase::ONE,
            becp_k,
            npm,
            ps_pointer,
            npm,
            &ModuleBase::ONE,
            h_tmp,
            nbands
        );
    }
    ModuleBase::timer::end("spinconstrain::SpinConstrain", "calculate_delta_hcc");
}

template <>
void spinconstrain::SpinConstrain<std::complex<double>>::update_psi_charge_pw_cpu(const ModuleBase::Vector3<double>* delta_lambda, bool pw_solve)
{
    ModuleBase::TITLE("spinconstrain::SpinConstrain", "update_psi_charge_pw_cpu");
    ModuleBase::timer::start("spinconstrain::SpinConstrain", "update_psi_charge_pw_cpu");
    
    psi::Psi<std::complex<double>>* psi_t = static_cast<psi::Psi<std::complex<double>>*>(this->psi);
    hamilt::Hamilt<std::complex<double>, base_device::DEVICE_CPU>* hamilt_t = static_cast<hamilt::Hamilt<std::complex<double>, base_device::DEVICE_CPU>*>(this->p_hamilt);
    auto* onsite_p = projectors::OnsiteProjector<double, base_device::DEVICE_CPU>::get_instance();
    
    int nbands = psi_t->get_nbands();
    int npol = psi_t->get_npol();
    int nkb = onsite_p->get_tot_nproj();
    int nk = psi_t->get_nk();
    int size_becp = nbands * nkb * npol;
    const int* nh_iat = &onsite_p->get_nh(0);
    
    std::vector<std::complex<double>> h_tmp(nbands * nbands), s_tmp(nbands * nbands);
    
    assert(this->sub_h_save != nullptr);
    assert(this->sub_s_save != nullptr);
    assert(this->becp_save != nullptr);
    
    for (int ik = 0; ik < nk; ++ik)
    {
        std::complex<double>* h_k = this->sub_h_save + ik * nbands * nbands;
        std::complex<double>* s_k = this->sub_s_save + ik * nbands * nbands;
        std::complex<double>* becp_k = this->becp_save + ik * size_becp;

        psi_t->fix_k(ik);
        
        memcpy(h_tmp.data(), h_k, sizeof(std::complex<double>) * nbands * nbands);
        memcpy(s_tmp.data(), s_k, sizeof(std::complex<double>) * nbands * nbands);
        
        // Apply DeltaSpin correction: H' = H_k + delta_H(lambda)
        this->calculate_delta_hcc(h_tmp.data(), becp_k, delta_lambda, nbands, nkb, nh_iat, ik);
        
        // Diagonalize in subspace to update wavefunction
        hsolver::DiagoIterAssist<std::complex<double>>::diag_subspace_psi(h_tmp.data(),
                                                                        s_tmp.data(),
                                                                        nbands,
                                                                        psi_t[0],
                                                                        &this->pelec->ekb(ik, 0));
    }

    // Clean up saved subspace data
    delete[] this->sub_h_save;
    delete[] this->sub_s_save;
    delete[] this->becp_save;
    this->sub_h_save = nullptr;
    this->sub_s_save = nullptr;
    this->becp_save = nullptr;

    // Subspace diagonalization already includes DeltaSpin correction via calculate_delta_hcc.
    // For the PW case, the full-space HSolverPW does NOT include the DeltaSpin correction
    // (it only exists in the subspace), so calling HSolverPW::solve would overwrite the
    // corrected psi with an uncorrected one, causing density explosion. Always use psiToRho.
    reinterpret_cast<elecstate::ElecStatePW<std::complex<double>, base_device::DEVICE_CPU>*>(this->pelec)->psiToRho(*psi_t);
    ModuleBase::timer::end("spinconstrain::SpinConstrain", "update_psi_charge_pw_cpu");
}

#if ((defined __CUDA) || (defined __ROCM))
template <>
void spinconstrain::SpinConstrain<std::complex<double>>::update_psi_charge_pw_gpu(const ModuleBase::Vector3<double>* delta_lambda, bool pw_solve)
{
    ModuleBase::TITLE("spinconstrain::SpinConstrain", "update_psi_charge_pw_gpu");
    ModuleBase::timer::start("spinconstrain::SpinConstrain", "update_psi_charge_pw_gpu");
    
    psi::Psi<std::complex<double>, base_device::DEVICE_GPU>* psi_t = static_cast<psi::Psi<std::complex<double>, base_device::DEVICE_GPU>*>(this->psi);
    hamilt::Hamilt<std::complex<double>, base_device::DEVICE_GPU>* hamilt_t = static_cast<hamilt::Hamilt<std::complex<double>, base_device::DEVICE_GPU>*>(this->p_hamilt);
    auto* onsite_p = projectors::OnsiteProjector<double, base_device::DEVICE_GPU>::get_instance();
    
    int nbands = psi_t->get_nbands();
    int npol = psi_t->get_npol();
    int nkb = onsite_p->get_tot_nproj();
    int nk = psi_t->get_nk();
    int size_becp = nbands * nkb * npol;
    const int* nh_iat = &onsite_p->get_nh(0);
    
    std::complex<double>* h_tmp = nullptr;
    std::complex<double>* s_tmp = nullptr;
    base_device::memory::resize_memory_op<std::complex<double>, base_device::DEVICE_GPU>()(h_tmp, nbands * nbands);
    base_device::memory::resize_memory_op<std::complex<double>, base_device::DEVICE_GPU>()(s_tmp, nbands * nbands);
    
    assert(this->sub_h_save != nullptr);
    assert(this->sub_s_save != nullptr);
    assert(this->becp_save != nullptr);
    
    for (int ik = 0; ik < nk; ++ik)
    {
        std::complex<double>* h_k = this->sub_h_save + ik * nbands * nbands;
        std::complex<double>* s_k = this->sub_s_save + ik * nbands * nbands;
        std::complex<double>* becp_k = this->becp_save + ik * size_becp;

        psi_t->fix_k(ik);
        
        base_device::memory::synchronize_memory_op<std::complex<double>, base_device::DEVICE_GPU, base_device::DEVICE_GPU>()(h_tmp, h_k, nbands * nbands);
        base_device::memory::synchronize_memory_op<std::complex<double>, base_device::DEVICE_GPU, base_device::DEVICE_GPU>()(s_tmp, s_k, nbands * nbands);
        
        // Apply DeltaSpin correction: H' = H_k + delta_H(lambda)
        this->calculate_delta_hcc(h_tmp, becp_k, delta_lambda, nbands, nkb, nh_iat, ik);
        
        // Diagonalize in subspace to update wavefunction
        hsolver::DiagoIterAssist<std::complex<double>, base_device::DEVICE_GPU>::diag_subspace_psi(h_tmp,
                                                                                s_tmp,
                                                                                nbands,
                                                                                psi_t[0],
                                                                                &this->pelec->ekb(ik, 0));
    }

    // Clean up saved subspace data
    base_device::memory::delete_memory_op<std::complex<double>, base_device::DEVICE_GPU>()(sub_h_save);
    base_device::memory::delete_memory_op<std::complex<double>, base_device::DEVICE_GPU>()(sub_s_save);
    base_device::memory::delete_memory_op<std::complex<double>, base_device::DEVICE_GPU>()(becp_save);
    this->sub_h_save = nullptr;
    this->sub_s_save = nullptr;
    this->becp_save = nullptr;

    // Subspace diagonalization already includes DeltaSpin correction via calculate_delta_hcc.
    // For the PW case, the full-space HSolverPW does NOT include the DeltaSpin correction,
    // so calling HSolverPW::solve would overwrite the corrected psi. Always use psiToRho.
    reinterpret_cast<elecstate::ElecStatePW<std::complex<double>, base_device::DEVICE_GPU>*>(this->pelec)->psiToRho(*psi_t);
    ModuleBase::timer::end("spinconstrain::SpinConstrain", "update_psi_charge_pw_gpu");
}
#endif

template <>
void spinconstrain::SpinConstrain<std::complex<double>>::cal_mw_from_lambda(
		int i_step, 
		const ModuleBase::Vector3<double>* delta_lambda)
{
    ModuleBase::TITLE("spinconstrain::SpinConstrain", "cal_mw_from_lambda");
    ModuleBase::timer::start("spinconstrain::SpinConstrain", "cal_mw_from_lambda");
    
#ifdef __LCAO
    if (PARAM.inp.basis_type == "lcao")
    {
        psi::Psi<std::complex<double>>* psi_t = static_cast<psi::Psi<std::complex<double>>*>(this->psi);
        hamilt::Hamilt<std::complex<double>>* hamilt_t = static_cast<hamilt::Hamilt<std::complex<double>>*>(this->p_hamilt);
        hsolver::HSolverLCAO<std::complex<double>> hsolver_t(this->ParaV, PARAM.inp.ks_solver);
        if (this->nspin_ == 2)
        {
            dynamic_cast<hamilt::DeltaSpin<hamilt::OperatorLCAO<std::complex<double>, double>>*>(this->p_operator)
                ->update_lambda();
        }
        else if (this->nspin_ == 4)
        {
            dynamic_cast<hamilt::DeltaSpin<hamilt::OperatorLCAO<std::complex<double>, std::complex<double>>>*>(
                this->p_operator)
                ->update_lambda();
        }
        // diagonalization without update charge
        hsolver_t.solve(hamilt_t, psi_t[0], this->pelec, *this->dm_, *this->pelec->charge, this->nspin_, true);
        elecstate::calculate_weights(this->pelec->ekb,
                                     this->pelec->wg,
                                     this->pelec->klist,
                                     this->pelec->eferm,
                                     this->pelec->f_en,
                                     this->pelec->nelec_spin,
                                     this->pelec->skip_weights);
        elecstate::calEBand(this->pelec->ekb,this->pelec->wg,this->pelec->f_en);

        elecstate::cal_dm_psi(this->ParaV, this->pelec->wg, *psi_t, *this->dm_);

        this->dm_->cal_DMR();

        this->cal_mi_lcao(i_step);
    }
    else
#endif
    {
        {
            this->zero_Mi();
            int size_becp = 0;
            std::vector<std::complex<double>> becp_tmp;
            int nk = 0;
            int nkb = 0;
            int nbands = 0;
            int npol = 0;
            const int* nh_iat = nullptr;
            if (PARAM.inp.device == "cpu")
            {
                psi::Psi<std::complex<double>>* psi_t = static_cast<psi::Psi<std::complex<double>>*>(this->psi);
                hamilt::Hamilt<std::complex<double>, base_device::DEVICE_CPU>* hamilt_t = static_cast<hamilt::Hamilt<std::complex<double>, base_device::DEVICE_CPU>*>(this->p_hamilt);
                auto* onsite_p = projectors::OnsiteProjector<double, base_device::DEVICE_CPU>::get_instance();
                nbands = psi_t->get_nbands();
                npol = psi_t->get_npol();
                nkb = onsite_p->get_tot_nproj();
                nk = psi_t->get_nk();
                nh_iat = &onsite_p->get_nh(0);
                size_becp = nbands * nkb * npol;
                becp_tmp.resize(size_becp * nk);
                std::vector<std::complex<double>> h_tmp(nbands * nbands), s_tmp(nbands * nbands);
                int initial_hs = 0;
                if(this->sub_h_save == nullptr)
                {
                    initial_hs = 1;
                    this->sub_h_save = new std::complex<double>[nbands * nbands * nk];
                    this->sub_s_save = new std::complex<double>[nbands * nbands * nk];
                    this->becp_save = new std::complex<double>[size_becp * nk];
                }
                for (int ik = 0; ik < nk; ++ik)
                {

                    psi_t->fix_k(ik);

                    std::complex<double>* h_k = this->sub_h_save + ik * nbands * nbands;
                    std::complex<double>* s_k = this->sub_s_save + ik * nbands * nbands;
                    std::complex<double>* becp_k = this->becp_save + ik * size_becp;
                    if(initial_hs)
                    {
                        /// update H(k) for each k point
                        hamilt_t->updateHk(ik);
                        hsolver::DiagoIterAssist<std::complex<double>>::cal_hs_subspace(hamilt_t, psi_t[0], h_k, s_k);
                        memcpy(becp_k, onsite_p->get_becp(), sizeof(std::complex<double>) * size_becp);
                    }
                    memcpy(h_tmp.data(), h_k, sizeof(std::complex<double>) * nbands * nbands);
                    memcpy(s_tmp.data(), s_k, sizeof(std::complex<double>) * nbands * nbands);
                    // update h_tmp by delta_lambda
                    if (i_step != -1) this->calculate_delta_hcc(h_tmp.data(), becp_k, delta_lambda, nbands, nkb, nh_iat, ik);

                    hsolver::DiagoIterAssist<std::complex<double>>::diag_responce(h_tmp.data(),
                                                                                  s_tmp.data(),
                                                                                  nbands,
                                                                                  becp_k,
                                                                                  &becp_tmp[ik * size_becp],
                                                                                  nkb * npol,
                                                                                  &this->pelec->ekb(ik, 0));
                }
            }
#if ((defined __CUDA) || (defined __ROCM))
            else
            {
                psi::Psi<std::complex<double>, base_device::DEVICE_GPU>* psi_t = static_cast<psi::Psi<std::complex<double>, base_device::DEVICE_GPU>*>(this->psi);
                hamilt::Hamilt<std::complex<double>, base_device::DEVICE_GPU>* hamilt_t = static_cast<hamilt::Hamilt<std::complex<double>, base_device::DEVICE_GPU>*>(this->p_hamilt);
                auto* onsite_p = projectors::OnsiteProjector<double, base_device::DEVICE_GPU>::get_instance();
                nbands = psi_t->get_nbands();
                npol = psi_t->get_npol();
                nkb = onsite_p->get_tot_nproj();
                nk = psi_t->get_nk();
                nh_iat = &onsite_p->get_nh(0);
                size_becp = nbands * nkb * npol;
                becp_tmp.resize(size_becp * nk);
                std::complex<double>* h_tmp = nullptr;
                std::complex<double>* s_tmp = nullptr;
                base_device::memory::resize_memory_op<std::complex<double>, base_device::DEVICE_GPU>()(h_tmp, nbands * nbands);
                base_device::memory::resize_memory_op<std::complex<double>, base_device::DEVICE_GPU>()(s_tmp, nbands * nbands);
                int initial_hs = 0;
                if(this->sub_h_save == nullptr)
                {
                    initial_hs = 1;
                    base_device::memory::resize_memory_op<std::complex<double>, base_device::DEVICE_GPU>()(this->sub_h_save, nbands * nbands * nk);
                    base_device::memory::resize_memory_op<std::complex<double>, base_device::DEVICE_GPU>()(this->sub_s_save, nbands * nbands * nk);
                    base_device::memory::resize_memory_op<std::complex<double>, base_device::DEVICE_GPU>()(this->becp_save, size_becp * nk);
                }
                std::complex<double>* becp_pointer = nullptr;
                base_device::memory::resize_memory_op<std::complex<double>, base_device::DEVICE_GPU>()(becp_pointer, size_becp);
                for (int ik = 0; ik < nk; ++ik)
                {
                    psi_t->fix_k(ik);

                    std::complex<double>* h_k = this->sub_h_save + ik * nbands * nbands;
                    std::complex<double>* s_k = this->sub_s_save + ik * nbands * nbands;
                    std::complex<double>* becp_k = this->becp_save + ik * size_becp;
                    if(initial_hs)
                    {
                        hamilt_t->updateHk(ik);
                        hsolver::DiagoIterAssist<std::complex<double>, base_device::DEVICE_GPU>::cal_hs_subspace(hamilt_t, psi_t[0], h_k, s_k);
                        base_device::memory::synchronize_memory_op<std::complex<double>, base_device::DEVICE_GPU, base_device::DEVICE_GPU>()(becp_k, onsite_p->get_becp(), size_becp);
                    }
                    base_device::memory::synchronize_memory_op<std::complex<double>, base_device::DEVICE_GPU, base_device::DEVICE_GPU>()(h_tmp, h_k, nbands * nbands);
                    base_device::memory::synchronize_memory_op<std::complex<double>, base_device::DEVICE_GPU, base_device::DEVICE_GPU>()(s_tmp, s_k, nbands * nbands);
                    if (i_step != -1) this->calculate_delta_hcc(h_tmp, becp_k, delta_lambda, nbands, nkb, nh_iat, ik);

                    hsolver::DiagoIterAssist<std::complex<double>, base_device::DEVICE_GPU>::diag_responce(h_tmp,
                                                                                  s_tmp,
                                                                                  nbands,
                                                                                  becp_k,
                                                                                  becp_pointer,
                                                                                  nkb * npol,
                                                                                  &this->pelec->ekb(ik, 0));
                    base_device::memory::synchronize_memory_op<std::complex<double>, base_device::DEVICE_CPU, base_device::DEVICE_GPU>()(&becp_tmp[ik * size_becp], becp_pointer, size_becp);   
                }

                base_device::memory::delete_memory_op<std::complex<double>, base_device::DEVICE_GPU>()(becp_pointer);
            }
#endif

            // calculate weights from ekb to update wg
            elecstate::calculate_weights(this->pelec->ekb,
                                         this->pelec->wg,
                                         this->pelec->klist,
                                         this->pelec->eferm,
                                         this->pelec->f_en,
                                         this->pelec->nelec_spin,
                                         this->pelec->skip_weights);
            // calculate Mi from existed becp
            for (int ik = 0; ik < nk; ik++)
            {
                const std::complex<double>* becp = &becp_tmp[ik * size_becp];
                this->accumulate_Mi_from_becp(becp, nkb, nbands, this->npol_, ik,
                    &this->pelec->wg(ik, 0), nh_iat);
            }
            Parallel_Reduce::reduce_double_allpool(PARAM.inp.kpar,
                                                    GlobalV::NPROC_IN_POOL,
                                                    &(this->Mi_[0][0]),
                                                    3 * this->Mi_.size());
        }
    }
    ModuleBase::timer::end("spinconstrain::SpinConstrain", "cal_mw_from_lambda");
}

template <>
void spinconstrain::SpinConstrain<std::complex<double>>::update_psi_charge(const ModuleBase::Vector3<double>* delta_lambda, bool pw_solve)
{
    ModuleBase::TITLE("spinconstrain::SpinConstrain", "update_psi_charge");
    ModuleBase::timer::start("spinconstrain::SpinConstrain", "update_psi_charge");
#ifdef __LCAO
    if (PARAM.inp.basis_type == "lcao")
    {
        psi::Psi<std::complex<double>>* psi_t = static_cast<psi::Psi<std::complex<double>>*>(this->psi);
        this->pelec->psiToRho(*psi_t);
    }
    else
#endif
    {
        if (PARAM.inp.device == "cpu")
        {
            this->update_psi_charge_pw_cpu(delta_lambda, pw_solve);
        }
#if ((defined __CUDA) || (defined __ROCM))
        else
        {
            this->update_psi_charge_pw_gpu(delta_lambda, pw_solve);
        }
#endif
    }
    ModuleBase::timer::end("spinconstrain::SpinConstrain", "update_psi_charge");
}
