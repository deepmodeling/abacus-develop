#include "module_base/timer.h"
#include "module_base/tool_title.h"
#include "module_hsolver/diago_iter_assist.h"
#include "module_parameter/parameter.h"
#include "spin_constrain.h"
#include "module_hamilt_pw/hamilt_pwdft/onsite_projector.h"
#include "module_base/parallel_reduce.h"
#include "module_base/kernels/math_kernel_op.h"
#include "module_hsolver/hsolver_lcao.h"
#include "module_hsolver/hsolver_pw.h"
#include "module_elecstate/elecstate_pw.h"

#ifdef __LCAO
#include "module_elecstate/elecstate_lcao.h"
#include "module_elecstate/module_dm/cal_dm_psi.h"
#include "module_hamilt_lcao/hamilt_lcaodft/operator_lcao/dspin_lcao.h"
#endif

template <>
void spinconstrain::SpinConstrain<std::complex<double>>::calculate_delta_hcc(std::complex<double>* h_tmp, const std::complex<double>* becp_k, const ModuleBase::Vector3<double>* delta_lambda, const int nbands, const int nkb, const int* nh_iat)
{
    int sum = 0;
    int size_ps = nkb * 2 * nbands;
    std::complex<double>* becp_cpu = nullptr;
    if(PARAM.inp.device == "gpu")
    {
#if ((defined __CUDA) || (defined __ROCM))
        base_device::DEVICE_GPU* ctx = {};
        base_device::DEVICE_CPU* cpu_ctx = {};
        base_device::memory::resize_memory_op<std::complex<double>, base_device::DEVICE_CPU>()(becp_cpu, size_ps);
        base_device::memory::synchronize_memory_op<std::complex<double>, base_device::DEVICE_CPU, base_device::DEVICE_GPU>()(becp_cpu, becp_k, size_ps);   
#endif
    }
    else if (PARAM.inp.device == "cpu")
    {
        becp_cpu = const_cast<std::complex<double>*>(becp_k);
    }

    std::vector<std::complex<double>> ps(size_ps, 0.0);
    for (int iat = 0; iat < this->Mi_.size(); iat++)
    {
        const int nproj = nh_iat[iat];
        const std::complex<double> coefficients0(delta_lambda[iat][2], 0.0);
        const std::complex<double> coefficients1(delta_lambda[iat][0] , delta_lambda[iat][1]);
        const std::complex<double> coefficients2(delta_lambda[iat][0] , -1 * delta_lambda[iat][1]);
        const std::complex<double> coefficients3(-1 * delta_lambda[iat][2], 0.0);
        // each atom has nproj, means this is with structure factor;
        // each projector (each atom) must multiply coefficient
        // with all the other projectors.
        for (int ib = 0; ib < nbands * 2; ib+=2)
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
            } // end ip
        } // end ib
        sum += nproj;
    } // end iat
    std::complex<double>* ps_pointer = nullptr;
    if(PARAM.inp.device == "gpu")
    {
#if ((defined __CUDA) || (defined __ROCM))
        base_device::DEVICE_GPU* ctx = {};
        base_device::DEVICE_CPU* cpu_ctx = {};
        base_device::memory::resize_memory_op<std::complex<double>, base_device::DEVICE_GPU>()(ps_pointer, size_ps);
        base_device::memory::synchronize_memory_op<std::complex<double>, base_device::DEVICE_GPU, base_device::DEVICE_CPU>()(ps_pointer, ps.data(), size_ps);   
#endif
    }
    else if (PARAM.inp.device == "cpu")
    {
        ps_pointer = ps.data();
    }
    // update h_tmp by becp_k * ps
    char transa = 'C';
    char transb = 'N';
    const int npm = nkb * 2;
    if (PARAM.inp.device == "gpu")
    {
#if ((defined __CUDA) || (defined __ROCM))
        base_device::DEVICE_GPU* ctx = {};
        ModuleBase::gemm_op<std::complex<double>, base_device::DEVICE_GPU>()(
            ctx,
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
        delete[] becp_cpu;
#endif

    }
    else if (PARAM.inp.device == "cpu")
    {
        base_device::DEVICE_CPU* ctx = {};
        ModuleBase::gemm_op<std::complex<double>, base_device::DEVICE_CPU>()(
            ctx,
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
}

template <>
void spinconstrain::SpinConstrain<std::complex<double>>::cal_mw_from_lambda(int i_step, const ModuleBase::Vector3<double>* delta_lambda)
{
    ModuleBase::TITLE("spinconstrain::SpinConstrain", "cal_mw_from_lambda");
    ModuleBase::timer::tick("spinconstrain::SpinConstrain", "cal_mw_from_lambda");
    // lambda has been updated in the lambda loop
#ifdef __LCAO
    if (PARAM.inp.basis_type == "lcao")
    {
        psi::Psi<std::complex<double>>* psi_t = static_cast<psi::Psi<std::complex<double>>*>(this->psi);
        hamilt::Hamilt<std::complex<double>>* hamilt_t = static_cast<hamilt::Hamilt<std::complex<double>>*>(this->p_hamilt);
        hsolver::HSolverLCAO<std::complex<double>> hsolver_t(this->ParaV, PARAM.inp.ks_solver);
        if (PARAM.inp.nspin == 2)
        {
            dynamic_cast<hamilt::DeltaSpin<hamilt::OperatorLCAO<std::complex<double>, double>>*>(this->p_operator)
                ->update_lambda();
        }
        else if (PARAM.inp.nspin == 4)
        {
            dynamic_cast<hamilt::DeltaSpin<hamilt::OperatorLCAO<std::complex<double>, std::complex<double>>>*>(
                this->p_operator)
                ->update_lambda();
        }
        // diagonalization without update charge
        hsolver_t.solve(hamilt_t, psi_t[0], this->pelec, true);
        this->pelec->calculate_weights();
        this->pelec->calEBand();
        elecstate::ElecStateLCAO<std::complex<double>>* pelec_lcao
            = dynamic_cast<elecstate::ElecStateLCAO<std::complex<double>>*>(this->pelec);
        elecstate::cal_dm_psi(this->ParaV, pelec_lcao->wg, *psi_t, *(pelec_lcao->get_DM()));
        pelec_lcao->get_DM()->cal_DMR();
        this->cal_mi_lcao(i_step);
    }
    else
#endif
    {
        /*if (i_step == -1 && this->higher_mag_prec)
        {
            // std::cout<<__FILE__<<__LINE__<<"istep == 0"<<std::endl;
            if (PARAM.inp.device == "cpu")
            {
                psi::Psi<std::complex<double>>* psi_t = static_cast<psi::Psi<std::complex<double>>*>(this->psi);
                hamilt::Hamilt<std::complex<double>>* hamilt_t = static_cast<hamilt::Hamilt<std::complex<double>>*>(this->p_hamilt);
                hsolver::HSolver<std::complex<double>, base_device::DEVICE_CPU>* hsolver_t = static_cast<hsolver::HSolver<std::complex<double>, base_device::DEVICE_CPU>*>(this->phsol);
                hsolver_t->solve(hamilt_t, psi_t[0], this->pelec, this->KS_SOLVER, true);
            }
            else
            {
                psi::Psi<std::complex<double>, base_device::DEVICE_GPU>* psi_t = static_cast<psi::Psi<std::complex<double>, base_device::DEVICE_GPU>*>(this->psi);
                hamilt::Hamilt<std::complex<double>, base_device::DEVICE_GPU>* hamilt_t = static_cast<hamilt::Hamilt<std::complex<double>, base_device::DEVICE_GPU>*>(this->p_hamilt);
                hsolver::HSolver<std::complex<double>, base_device::DEVICE_GPU>* hsolver_t = static_cast<hsolver::HSolver<std::complex<double>, base_device::DEVICE_GPU>*>(this->phsol);
                hsolver_t->solve(hamilt_t, psi_t[0], this->pelec, this->KS_SOLVER, true);
            }
            this->pelec->calculate_weights();
            this->cal_Mi_pw();
        }
        else*/
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
                    if (i_step != -1) this->calculate_delta_hcc(h_tmp.data(), becp_k, delta_lambda, nbands, nkb, nh_iat);

                    hsolver::DiagoIterAssist<std::complex<double>>::diag_responce(h_tmp.data(),
                                                                                  s_tmp.data(),
                                                                                  nbands,
                                                                                  becp_k,
                                                                                  &becp_tmp[ik * size_becp],
                                                                                  nkb * 2,
                                                                                  &this->pelec->ekb(ik, 0));
                }
            }
#if ((defined __CUDA) || (defined __ROCM))
            else
            {
                base_device::DEVICE_GPU* ctx = {};
                base_device::DEVICE_CPU* cpu_ctx = {};
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
                // allocate memory for becp_pointer in GPU device
                base_device::memory::resize_memory_op<std::complex<double>, base_device::DEVICE_GPU>()(becp_pointer, size_becp);
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
                        hsolver::DiagoIterAssist<std::complex<double>, base_device::DEVICE_GPU>::cal_hs_subspace(hamilt_t, psi_t[0], h_k, s_k);
                        base_device::memory::synchronize_memory_op<std::complex<double>, base_device::DEVICE_GPU, base_device::DEVICE_GPU>()(becp_k, onsite_p->get_becp(), size_becp);
                    }
                    base_device::memory::synchronize_memory_op<std::complex<double>, base_device::DEVICE_GPU, base_device::DEVICE_GPU>()(h_tmp, h_k, nbands * nbands);
                    base_device::memory::synchronize_memory_op<std::complex<double>, base_device::DEVICE_GPU, base_device::DEVICE_GPU>()(s_tmp, s_k, nbands * nbands);
                    // update h_tmp by delta_lambda
                    if (i_step != -1) this->calculate_delta_hcc(h_tmp, becp_k, delta_lambda, nbands, nkb, nh_iat);

                    hsolver::DiagoIterAssist<std::complex<double>, base_device::DEVICE_GPU>::diag_responce(h_tmp,
                                                                                  s_tmp,
                                                                                  nbands,
                                                                                  becp_k,
                                                                                  becp_pointer,
                                                                                  nkb * npol,
                                                                                  &this->pelec->ekb(ik, 0));
                    // copy becp_pointer from GPU to CPU
                    base_device::memory::synchronize_memory_op<std::complex<double>, base_device::DEVICE_CPU, base_device::DEVICE_GPU>()(&becp_tmp[ik * size_becp], becp_pointer, size_becp);   
                }

                // free memory for becp_pointer in GPU device
                base_device::memory::delete_memory_op<std::complex<double>, base_device::DEVICE_GPU>()(becp_pointer);
            }
#endif
            // calculate weights from ekb to update wg
            this->pelec->calculate_weights();
            // calculate Mi from existed becp
            for (int ik = 0; ik < nk; ik++)
            {
                const std::complex<double>* becp = &becp_tmp[ik * size_becp];
                // becp(nbands*npol , nkb)
                // mag = wg * \sum_{nh}becp * becp
                for (int ib = 0; ib < nbands; ib++)
                {
                    const double weight = this->pelec->wg(ik, ib);
                    int begin_ih = 0;
                    for (int iat = 0; iat < this->Mi_.size(); iat++)
                    {
                        const int nh = nh_iat[iat];
                        std::complex<double> occ[4]
                            = {ModuleBase::ZERO, ModuleBase::ZERO, ModuleBase::ZERO, ModuleBase::ZERO};
                        for (int ih = 0; ih < nh; ih++)
                        {
                            const int index = ib * npol * nkb + begin_ih + ih;
                            occ[0] += conj(becp[index]) * becp[index];
                            occ[1] += conj(becp[index]) * becp[index + nkb];
                            occ[2] += conj(becp[index + nkb]) * becp[index];
                            occ[3] += conj(becp[index + nkb]) * becp[index + nkb];
                        }
                        // occ has been reduced and calculate mag
                        this->Mi_[iat].x += weight * (occ[1] + occ[2]).real();
                        this->Mi_[iat].y += weight * (occ[1] - occ[2]).imag();
                        this->Mi_[iat].z += weight * (occ[0] - occ[3]).real();
                        begin_ih += nh;
                    }
                }
            }
            Parallel_Reduce::reduce_double_allpool(GlobalV::KPAR,
                                                   GlobalV::NPROC_IN_POOL,
                                                   &(this->Mi_[0][0]),
                                                   3 * this->Mi_.size());
            // for(int i = 0; i < this->Mi_.size(); i++)
            //{
            //     std::cout<<"atom"<<i<<": "<<" mag: "<<this->Mi_[i].x<<" "<<this->Mi_[i].y<<" "<<this->Mi_[i].z<<"
            //     "<<this->lambda_[i].x<<" "<<this->lambda_[i].y<<" "<<this->lambda_[i].z<<std::endl;
            // }
        }
    }
    ModuleBase::timer::tick("SpinConstrain", "cal_mw_from_lambda");
}

template <>
void spinconstrain::SpinConstrain<std::complex<double>>::update_psi_charge(const ModuleBase::Vector3<double>* delta_lambda, bool pw_solve)
{
    ModuleBase::TITLE("spinconstrain::SpinConstrain", "update_psi_charge");
    ModuleBase::timer::tick("spinconstrain::SpinConstrain", "update_psi_charge");
#ifdef __LCAO
    if (PARAM.inp.basis_type == "lcao")
    {
        psi::Psi<std::complex<double>>* psi_t = static_cast<psi::Psi<std::complex<double>>*>(this->psi);
        this->pelec->psiToRho(*psi_t);
    }
    else
#endif
    {
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
                this->calculate_delta_hcc(h_tmp.data(), becp_k, delta_lambda, nbands, nkb, nh_iat);
                hsolver::DiagoIterAssist<std::complex<double>>::diag_subspace_psi(h_tmp.data(),
                                                                                s_tmp.data(),
                                                                                nbands,
                                                                                psi_t[0],
                                                                                &this->pelec->ekb(ik, 0));
            }

            delete[] this->sub_h_save;
            delete[] this->sub_s_save;
            delete[] this->becp_save;
            this->sub_h_save = nullptr;
            this->sub_s_save = nullptr;
            this->becp_save = nullptr;

            if(pw_solve)
            {
                hsolver::HSolverPW<std::complex<double>, base_device::DEVICE_CPU> hsolver_pw_obj(this->pw_wfc_,
                                                 PARAM.inp.calculation,
                                                 PARAM.inp.basis_type,
                                                 PARAM.inp.ks_solver,
                                                 PARAM.inp.use_paw,
                                                 PARAM.globalv.use_uspp,
                                                 PARAM.inp.nspin,
                                                 hsolver::DiagoIterAssist<std::complex<double>, base_device::DEVICE_CPU>::SCF_ITER,
                                                 hsolver::DiagoIterAssist<std::complex<double>, base_device::DEVICE_CPU>::PW_DIAG_NMAX,
                                                 hsolver::DiagoIterAssist<std::complex<double>, base_device::DEVICE_CPU>::PW_DIAG_THR,
                                                 hsolver::DiagoIterAssist<std::complex<double>, base_device::DEVICE_CPU>::need_subspace);

                hsolver_pw_obj.solve(hamilt_t,
                         psi_t[0],
                         this->pelec,
                         this->pelec->ekb.c,
                         GlobalV::RANK_IN_POOL,
                         GlobalV::NPROC_IN_POOL,
                         false,
                         this->tpiba,
                         this->get_nat());
            }
            else
            {// update charge density only
                this->pelec->psiToRho(*psi_t);
            }
        }
#if ((defined __CUDA) || (defined __ROCM))
        else
        {
            base_device::DEVICE_GPU* ctx = {};
            base_device::DEVICE_CPU* cpu_ctx = {};
            psi::Psi<std::complex<double>, base_device::DEVICE_GPU>* psi_t = static_cast<psi::Psi<std::complex<double>, base_device::DEVICE_GPU>*>(this->psi);
            hamilt::Hamilt<std::complex<double>, base_device::DEVICE_GPU>* hamilt_t = static_cast<hamilt::Hamilt<std::complex<double>, base_device::DEVICE_GPU>*>(this->p_hamilt);
            auto* onsite_p = projectors::OnsiteProjector<double, base_device::DEVICE_GPU>::get_instance();
            nbands = psi_t->get_nbands();
            npol = psi_t->get_npol();
            nkb = onsite_p->get_tot_nproj();
            nk = psi_t->get_nk();
            nh_iat = &onsite_p->get_nh(0);
            size_becp = nbands * nkb * npol;

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
                this->calculate_delta_hcc(h_tmp, becp_k, delta_lambda, nbands, nkb, nh_iat);
                hsolver::DiagoIterAssist<std::complex<double>, base_device::DEVICE_GPU>::diag_subspace_psi(h_tmp,
                                                                                s_tmp,
                                                                                nbands,
                                                                                psi_t[0],
                                                                                &this->pelec->ekb(ik, 0));
            }

            base_device::memory::delete_memory_op<std::complex<double>, base_device::DEVICE_GPU>()(sub_h_save);
            base_device::memory::delete_memory_op<std::complex<double>, base_device::DEVICE_GPU>()(sub_s_save);
            base_device::memory::delete_memory_op<std::complex<double>, base_device::DEVICE_GPU>()(becp_save);
            this->sub_h_save = nullptr;
            this->sub_s_save = nullptr;
            this->becp_save = nullptr;

            if(pw_solve)
            {
                hsolver::HSolverPW<std::complex<double>, base_device::DEVICE_GPU> hsolver_pw_obj(this->pw_wfc_,
                                                 PARAM.inp.calculation,
                                                 PARAM.inp.basis_type,
                                                 PARAM.inp.ks_solver,
                                                 PARAM.inp.use_paw,
                                                 PARAM.globalv.use_uspp,
                                                 PARAM.inp.nspin,
                                                 hsolver::DiagoIterAssist<std::complex<double>, base_device::DEVICE_GPU>::SCF_ITER,
                                                 hsolver::DiagoIterAssist<std::complex<double>, base_device::DEVICE_GPU>::PW_DIAG_NMAX,
                                                 hsolver::DiagoIterAssist<std::complex<double>, base_device::DEVICE_GPU>::PW_DIAG_THR,
                                                 hsolver::DiagoIterAssist<std::complex<double>, base_device::DEVICE_GPU>::need_subspace);

                hsolver_pw_obj.solve(hamilt_t,
                         psi_t[0],
                         this->pelec,
                         this->pelec->ekb.c,
                         GlobalV::RANK_IN_POOL,
                         GlobalV::NPROC_IN_POOL,
                         false,
                         this->tpiba,
                         this->get_nat());
            }
            else
            {// update charge density only
                reinterpret_cast<elecstate::ElecStatePW<std::complex<double>, base_device::DEVICE_GPU>*>(this->pelec)->psiToRho(*psi_t);
            }
            
        }
#endif       
    }
    ModuleBase::timer::tick("spinconstrain::SpinConstrain", "update_psi_charge");
}
