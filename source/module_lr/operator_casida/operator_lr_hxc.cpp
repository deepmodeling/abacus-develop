#include "operator_lr_hxc.h"
#include <vector>
#include "module_parameter/parameter.h"
#include "module_base/blas_connector.h"
#include "module_base/timer.h"
#include "module_lr/utils/lr_util.h"
#include "module_lr/utils/lr_util_hcontainer.h"
#include "module_lr/utils/lr_util_print.h"
// #include "module_hamilt_lcao/hamilt_lcaodft/DM_gamma_2d_to_grid.h"
#include "module_hamilt_lcao/module_hcontainer/hcontainer_funcs.h"
#include "module_lr/AX/AX.h"
#include "module_hamilt_pw/hamilt_pwdft/global.h"

inline double conj(double a) { return a; }
inline std::complex<double> conj(std::complex<double> a) { return std::conj(a); }

namespace LR
{
    template<typename T, typename Device>
    void OperatorLRHxc<T, Device>::act(const int nbands, const int nbasis, const int npol, const T* psi_in, T* hpsi, const int ngk_ik, const bool is_first_node)const
    {
        ModuleBase::TITLE("OperatorLRHxc", "act");
        const int& sl = ispin_ks[0];
        const auto psil_ks = LR_Util::get_psi_spin(psi_ks, sl, nk);
        const int& lgd = gint->gridt->lgd;

        this->DM_trans->cal_DMR();  //DM_trans->get_DMR_vector() is 2d-block parallized
        // LR_Util::print_DMR(*DM_trans, ucell.nat, "DMR");

        // ========================= begin grid calculation=========================
        this->grid_calculation(nbands);   //DM(R) to H(R)
        // ========================= end grid calculation =========================

        // V(R)->V(k) 
        std::vector<ct::Tensor> v_hxc_2d(nk, LR_Util::newTensor<T>({ pmat.get_col_size(), pmat.get_row_size() }));
        for (auto& v : v_hxc_2d) v.zero();
        int nrow = ModuleBase::GlobalFunc::IS_COLUMN_MAJOR_KS_SOLVER(PARAM.inp.ks_solver) ? this->pmat.get_row_size() : this->pmat.get_col_size();
        for (int ik = 0;ik < nk;++ik) { folding_HR(*this->hR, v_hxc_2d[ik].data<T>(), this->kv.kvec_d[ik], nrow, 1); }  // V(R) -> V(k)
        // LR_Util::print_HR(*this->hR, this->ucell.nat, "4.VR");
        // if (this->first_print)
        // for (int ik = 0;ik < nk;++ik)
        //     LR_Util::print_tensor<T>(v_hxc_2d[ik], "4.V(k)[ik=" + std::to_string(ik) + "]", &this->pmat);

        // 5. [AX]^{Hxc}_{ai}=\sum_{\mu,\nu}c^*_{a,\mu,}V^{Hxc}_{\mu,\nu}c_{\nu,i}
#ifdef __MPI
        cal_AX_pblas(v_hxc_2d, this->pmat, psil_ks, this->pc, naos, nocc[sl], nvirt[sl], this->pX[sl], hpsi);
#else
        cal_AX_blas(v_hxc_2d, psil_ks, nocc[sl], nvirt[sl], hpsi);
#endif
    }


    template<>
    void OperatorLRHxc<double, base_device::DEVICE_CPU>::grid_calculation(const int& nbands) const
    {
        ModuleBase::TITLE("OperatorLRHxc", "grid_calculation(real)");
        ModuleBase::timer::tick("OperatorLRHxc", "grid_calculation");
        this->gint->transfer_DM2DtoGrid(this->DM_trans->get_DMR_vector());     // 2d block to grid

        // 2. transition electron density
        // \f[ \tilde{\rho}(r)=\sum_{\mu_j, \mu_b}\tilde{\rho}_{\mu_j,\mu_b}\phi_{\mu_b}(r)\phi_{\mu_j}(r) \f]
        double** rho_trans;
        const int& nrxx = this->pot.lock()->nrxx;
        LR_Util::_allocate_2order_nested_ptr(rho_trans, 1, nrxx); // currently gint_kernel_rho uses PARAM.inp.nspin, it needs refactor
        ModuleBase::GlobalFunc::ZEROS(rho_trans[0], nrxx);
        Gint_inout inout_rho(rho_trans, Gint_Tools::job_type::rho, 1, false);
        this->gint->cal_gint(&inout_rho);

        // 3. v_hxc = f_hxc * rho_trans
        ModuleBase::matrix vr_hxc(1, nrxx);   //grid
        this->pot.lock()->cal_v_eff(rho_trans, GlobalC::ucell, vr_hxc, ispin_ks);
        LR_Util::_deallocate_2order_nested_ptr(rho_trans, 1);

        // 4. V^{Hxc}_{\mu,\nu}=\int{dr} \phi_\mu(r) v_{Hxc}(r) \phi_\mu(r)
        Gint_inout inout_vlocal(vr_hxc.c, 0, Gint_Tools::job_type::vlocal);
        this->gint->get_hRGint()->set_zero();
        this->gint->cal_gint(&inout_vlocal);
        this->hR->set_zero();   // clear hR for each bands
        this->gint->transfer_pvpR(&*this->hR, &GlobalC::ucell);    //grid to 2d block
        ModuleBase::timer::tick("OperatorLRHxc", "grid_calculation");
    }

    template<>
    void OperatorLRHxc<std::complex<double>, base_device::DEVICE_CPU>::grid_calculation(const int& nbands) const
    {
        ModuleBase::TITLE("OperatorLRHxc", "grid_calculation(complex)");
        ModuleBase::timer::tick("OperatorLRHxc", "grid_calculation");

        elecstate::DensityMatrix<std::complex<double>, double> DM_trans_real_imag(&pmat, 1, kv.kvec_d, kv.get_nks() / nspin);
        DM_trans_real_imag.init_DMR(*this->hR);
        hamilt::HContainer<double> HR_real_imag(GlobalC::ucell, &this->pmat);
        LR_Util::initialize_HR<std::complex<double>, double>(HR_real_imag, ucell, gd, orb_cutoff_);

        auto dmR_to_hR = [&, this](const char& type) -> void
            {
                LR_Util::get_DMR_real_imag_part(*this->DM_trans, DM_trans_real_imag, ucell.nat, type);
                // if (this->first_print)LR_Util::print_DMR(DM_trans_real_imag, ucell.nat, "DMR(2d, real)");

                this->gint->transfer_DM2DtoGrid(DM_trans_real_imag.get_DMR_vector());
                // LR_Util::print_HR(*this->gint->get_DMRGint()[0], this->ucell.nat, "DMR(grid, real)");

                // 2. transition electron density
                double** rho_trans;
                const int& nrxx = this->pot.lock()->nrxx;

                LR_Util::_allocate_2order_nested_ptr(rho_trans, 1, nrxx); // nspin=1 for transition density
                ModuleBase::GlobalFunc::ZEROS(rho_trans[0], nrxx);
                Gint_inout inout_rho(rho_trans, Gint_Tools::job_type::rho, 1, false);
                this->gint->cal_gint(&inout_rho);
                // print_grid_nonzero(rho_trans[0], nrxx, 10, "rho_trans");

                // 3. v_hxc = f_hxc * rho_trans
                ModuleBase::matrix vr_hxc(1, nrxx);   //grid
                this->pot.lock()->cal_v_eff(rho_trans, GlobalC::ucell, vr_hxc, ispin_ks);
                // print_grid_nonzero(vr_hxc.c, this->poticab->nrxx, 10, "vr_hxc");

                LR_Util::_deallocate_2order_nested_ptr(rho_trans, 1);

                // 4. V^{Hxc}_{\mu,\nu}=\int{dr} \phi_\mu(r) v_{Hxc}(r) \phi_\mu(r)
                Gint_inout inout_vlocal(vr_hxc.c, 0, Gint_Tools::job_type::vlocal);
                this->gint->get_hRGint()->set_zero();
                this->gint->cal_gint(&inout_vlocal);

                // LR_Util::print_HR(*this->gint->get_hRGint(), this->ucell.nat, "VR(grid)");
                HR_real_imag.set_zero();
                this->gint->transfer_pvpR(&HR_real_imag, &GlobalC::ucell, &GlobalC::GridD);
                // LR_Util::print_HR(HR_real_imag, this->ucell.nat, "VR(real, 2d)");
                LR_Util::set_HR_real_imag_part(HR_real_imag, *this->hR, GlobalC::ucell.nat, type);
            };
        this->hR->set_zero();
        dmR_to_hR('R');   //real
        if (kv.get_nks() / this->nspin > 1) { dmR_to_hR('I'); }   //imag for multi-k
        ModuleBase::timer::tick("OperatorLRHxc", "grid_calculation");
    }

    template class OperatorLRHxc<double>;
    template class OperatorLRHxc<std::complex<double>>;
}