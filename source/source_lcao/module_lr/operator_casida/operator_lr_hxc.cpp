#include "operator_lr_hxc.h"
#include <cstdlib>
#include <vector>
#include "source_io/module_parameter/parameter.h"
#include "source_base/timer.h"
#include "source_lcao/module_lr/utils/lr_util.h"
#include "source_lcao/module_lr/utils/lr_util_hcontainer.h"
#include "source_lcao/module_lr/utils/lr_util_print.h"
// #include "source_lcao/DM_gamma_2d_to_grid.h"
#include "source_hamilt/module_hcontainer/hcontainer_funcs.h"
#include "source_lcao/module_lr/ao_to_mo_transformer/ao_to_mo.h"
#include "source_hamilt/module_gint/gint_interface.h"
#include "source_lcao/module_lr/ao_to_mo_transformer/cvcx.h"

inline double conj(double a) { return a; }
inline std::complex<double> conj(std::complex<double> a) { return std::conj(a); }

namespace LR
{
    template<typename T, typename Device>
    void OperatorLRHxc<T, Device>::act(const int nbands, const int nbasis, const int npol, const T* psi_in, T* hpsi, const int ngk_ik, const bool is_first_node)const
    {
        TransitionDensityCache density_cache;
        this->act_with_shared_density(psi_in, hpsi, density_cache);
    }

    template<typename T, typename Device>
    void OperatorLRHxc<T, Device>::act_with_shared_density(
        const T* psi_in, T* hpsi, TransitionDensityCache& density_cache) const
    {
        ModuleBase::TITLE("OperatorLRHxc", "act_with_shared_density");
        ModuleBase::timer::start("OperatorLRHxc", "act_with_shared_density");

        const int& sl = ispin_ks[0];
        const auto psil_ks = LR_Util::get_psi_spin(psi_ks, sl, nk);

        if (density_cache.count(&this->DM_trans) == 0)
        {
            this->DM_trans.cal_dmr(-1);  // DM_trans.get_dmr_vec() is 2D-block parallelized.
            // cal_dmr preserves AO indices while converting DMK's storage layout to DMR.
        }
        // ========================= begin grid calculation =========================
        this->grid_calculation(density_cache);   // DM(R) -> rho(r) -> V_Hxc(r) -> H(R)
        // ========================= end grid calculation ===========================

        // V(R)->V(k) 
        std::vector<ct::Tensor> v_hxc_2d(nk, LR_Util::newTensor<T>({ pmat.get_col_size(), pmat.get_row_size() }));
        for (auto& v : v_hxc_2d) v.zero();
        int nrow = ModuleBase::GlobalFunc::IS_COLUMN_MAJOR_KS_SOLVER(PARAM.inp.ks_solver) ? this->pmat.get_row_size() : this->pmat.get_col_size();
        for (int ik = 0;ik < nk;++ik) { folding_HR(*this->hR, v_hxc_2d[ik].data<T>(), this->kv.kvec_d[ik], nrow, 1); }  // V(R) -> V(k)
        // LR_Util::print_HR(*this->hR, this->ucell.nat, "4.VR");
        // if (this->first_print)
        // for (int ik = 0;ik < nk;++ik)
        //     LR_Util::print_tensor<T>(v_hxc_2d[ik], "4.V(k)[ik=" + std::to_string(ik) + "]", &this->pmat);

        // 5. AO to MO transformation
        switch (this->dm_pq_)
        {
        case MO_TO_AO_TYPE::CC_vo:  //[AX]^{Hxc}_{ai}=\sum_{\mu,\nu}c^*_{\mu,a}V^{Hxc}_{\mu,\nu}c_{\nu,i}
#ifdef __MPI
            ao_to_mo_pblas(v_hxc_2d, this->pmat, psil_ks, this->pc, this->naos, this->nocc[sl], this->nvirt[sl], this->pX[sl], hpsi, /*add_on=*/true, LR_Util::MO_TYPE::VO, this->factor_);
#else
            ao_to_mo_blas(v_hxc_2d, psil_ks, this->nocc[sl], this->nvirt[sl], hpsi, /*add_on=*/true, LR_Util::MO_TYPE::VO, this->factor_);
#endif
            break;
        case MO_TO_AO_TYPE::CC_oo:  //[AX]^{Hxc}_{ij}=\sum_{\mu,\nu}c^*_{\mu,i}V^{Hxc}_{\mu,\nu}c_{\nu,j}
#ifdef __MPI
            ao_to_mo_pblas(v_hxc_2d, this->pmat, psil_ks, this->pc, this->naos, this->nocc[sl], this->nvirt[sl], this->pX[sl], hpsi, /*add_on=*/true, LR_Util::MO_TYPE::OO, this->factor_);
#else
            ao_to_mo_blas(v_hxc_2d, psil_ks, this->nocc[sl], this->nvirt[sl], hpsi, /*add_on=*/true, LR_Util::MO_TYPE::OO, this->factor_);
#endif
            break;
        case MO_TO_AO_TYPE::CXC:
        {
#ifdef __MPI
            CVCX_virt_pblas(v_hxc_2d, this->pmat, psil_ks, this->pc, psi_in, this->pX[sl],
                this->naos, this->nocc[sl], this->nvirt[sl], hpsi, /*add_on=*/true, this->factor_);
            CVCX_occ_pblas(v_hxc_2d, this->pmat, psil_ks, this->pc, psi_in, this->pX[sl],
                this->naos, this->nocc[sl], this->nvirt[sl], hpsi, /*add_on=*/true, -this->factor_);
#else
            CVCX_virt_blas(v_hxc_2d, psil_ks, psi_in, this->naos, this->nocc[sl], this->nvirt[sl], hpsi, /*add_on=*/true, this->factor_);
            CVCX_occ_blas(v_hxc_2d, psil_ks, psi_in, this->naos, this->nocc[sl], this->nvirt[sl], hpsi, /*add_on=*/true, -this->factor_);
#endif
            break;
        }
        case MO_TO_AO_TYPE::CXC_o:
#ifdef __MPI
            CVCX_occ_pblas(v_hxc_2d, this->pmat, psil_ks, this->pc, psi_in, this->pX[sl],
                this->naos, this->nocc[sl], this->nvirt[sl], hpsi, /*add_on=*/true, this->factor_);
#else
            CVCX_occ_blas(v_hxc_2d, psil_ks, psi_in, this->naos, this->nocc[sl], this->nvirt[sl], hpsi, /*add_on=*/true, this->factor_);
#endif
            break;
        default:
            throw std::runtime_error("Unknown DM_TYPE");
            break;
        }
        // for debug
        //std::cout << "After Hxc, hpsi: [nvirt= " << nvirt[sl] << " nocc= " << nocc[sl] << " nk= " << nk << " ]" << std::endl;
        //LR_Util::print_value(hpsi, nk, nocc[sl], nvirt[sl]);

        ModuleBase::timer::end("OperatorLRHxc", "act_with_shared_density");
    }


    template<>
    void OperatorLRHxc<double, base_device::DEVICE_CPU>::grid_calculation(TransitionDensityCache& density_cache) const
    {
        ModuleBase::TITLE("OperatorLRHxc", "grid_calculation(real)");
        ModuleBase::timer::start("OperatorLRHxc", "grid_calculation");

        // 2. transition electron density
        // \f[ \tilde{\rho}(r)=\sum_{\mu_j, \mu_b}\tilde{\rho}_{\mu_j,\mu_b}\phi_{\mu_b}(r)\phi_{\mu_j}(r) \f]
        // For a fixed input spin, rho is shared by both output-spin kernels.
        const int nrxx = this->pot.lock()->nrxx;
        auto& density = density_cache[&this->DM_trans];
        if (density.empty())
        {
            const std::vector<double> zeros(nrxx, 0.0);
            density.assign(1, zeros);  // nspin=1 for transition density
            double* rho = density[0].data();
            const auto& dmr = this->DM_trans.get_dmr_vec();
            ModuleGint::cal_gint_rho(dmr, 1, &rho, false);
        }
        // 3. v_hxc = f_hxc * rho_trans, evaluated separately for each output spin
        double* rho = density[0].data();
        ModuleBase::matrix vr_hxc(1, nrxx);   // grid
        this->pot.lock()->cal_v_eff(&rho, ucell, vr_hxc, ispin_ks);

        // 4. V^{Hxc}_{\mu,\nu}=\int{dr} \phi_\mu(r) v_{Hxc}(r) \phi_\nu(r)
        this->hR->set_zero();   // clear hR for each bands
        ModuleGint::cal_gint_vl(vr_hxc.c, &*this->hR);
        ModuleBase::timer::end("OperatorLRHxc", "grid_calculation");
    }

    template<>
    void OperatorLRHxc<std::complex<double>, base_device::DEVICE_CPU>::grid_calculation(TransitionDensityCache& density_cache) const
    {
        ModuleBase::TITLE("OperatorLRHxc", "grid_calculation(complex)");
        ModuleBase::timer::start("OperatorLRHxc", "grid_calculation");

        // 2. transition electron density, evaluated for each real/imaginary part
        // \f[ \tilde{\rho}(r)=\sum_{\mu_j, \mu_b}\tilde{\rho}_{\mu_j,\mu_b}\phi_{\mu_b}(r)\phi_{\mu_j}(r) \f]
        const int nrxx = this->pot.lock()->nrxx;
        const int nparts = nk > 1 ? 2 : 1;  // real; also imaginary for multi-k
        auto& density = density_cache[&this->DM_trans];
        if (density.empty())
        {
            module_dm::DensityMatrix<std::complex<double>, double> dm_real_imag(&pmat, 1, kv.kvec_d, nk);
            dm_real_imag.init_dmr(*this->hR);
            const std::vector<double> zeros(nrxx, 0.0);
            density.assign(nparts, zeros);
            for (int ipart = 0; ipart < nparts; ++ipart)
            {
                const char part = ipart == 0 ? 'R' : 'I';
                LR_Util::get_DMR_real_imag_part(this->DM_trans, dm_real_imag, part);
                // nspin=1 for each real/imaginary transition-density component
                double* rho = density[ipart].data();
                const auto& dmr = dm_real_imag.get_dmr_vec();
                ModuleGint::cal_gint_rho(dmr, 1, &rho, false);
            }
        }

        hamilt::HContainer<double> hr_real_imag(ucell, &this->pmat);
        LR_Util::initialize_HR<std::complex<double>, double>(hr_real_imag, ucell, gd, orb_cutoff_);
        this->hR->set_zero();   // clear hR for each band
        for (int ipart = 0; ipart < nparts; ++ipart)
        {
            const char part = ipart == 0 ? 'R' : 'I';
            // 3. v_hxc = f_hxc * rho_trans, evaluated separately for each output spin
            double* rho = density[ipart].data();
            ModuleBase::matrix vr_hxc(1, nrxx);   // grid
            this->pot.lock()->cal_v_eff(&rho, ucell, vr_hxc, ispin_ks);

            // 4. V^{Hxc}_{\mu,\nu}=\int{dr} \phi_\mu(r) v_{Hxc}(r) \phi_\nu(r)
            hr_real_imag.set_zero();
            ModuleGint::cal_gint_vl(vr_hxc.c, &hr_real_imag);
            LR_Util::set_HR_real_imag_part(hr_real_imag, *this->hR, part);
        }
        ModuleBase::timer::end("OperatorLRHxc", "grid_calculation");
    }

    template class OperatorLRHxc<double>;
    template class OperatorLRHxc<std::complex<double>>;
}
