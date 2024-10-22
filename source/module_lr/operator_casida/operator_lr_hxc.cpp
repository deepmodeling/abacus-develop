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
#include "module_lr/dm_trans/dm_trans.h"
#include "module_lr/AX/AX.h"
#include "module_hamilt_pw/hamilt_pwdft/global.h"

inline double conj(double a) { return a; }
inline std::complex<double> conj(std::complex<double> a) { return std::conj(a); }

namespace LR
{
    template<typename T, typename Device>
    void OperatorLRHxc<T, Device>::act(const psi::Psi<T>& psi_in, psi::Psi<T>& psi_out, const int nbands) const
    {
        ModuleBase::TITLE("OperatorLRHxc", "act");
        assert(nbands <= psi_in.get_nbands());
        const int& nk = this->kv.get_nks() / this->nspin;

        //print 
        // if (this->first_print) LR_Util::print_psi_kfirst(*psi_ks, "psi_ks");

        this->init_DM_trans(nbands, this->DM_trans);    // initialize transion density matrix

        psi::Psi<T> psi_in_bfirst = LR_Util::k1_to_bfirst_wrapper(psi_in, nk, this->pX->get_local_size());
        psi::Psi<T> psi_out_bfirst = LR_Util::k1_to_bfirst_wrapper(psi_out, nk, this->pX->get_local_size());

        const int& lgd = gint->gridt->lgd;
        for (int ib = 0;ib < nbands;++ib)
        {
            // if (this->first_print) LR_Util::print_psi_bandfirst(psi_in_bfirst, "psi_in_bfirst", ib);

            // if Hxc-only, the memory of single-band DM_trans is enough.
            // if followed by EXX, we need to allocate memory for all bands.
            int ib_dm = (this->next_op == nullptr) ? 0 : ib;
            psi_in_bfirst.fix_b(ib);
            psi_out_bfirst.fix_b(ib);

            // 1. transition density matrix
#ifdef __MPI
            std::vector<container::Tensor>  dm_trans_2d = cal_dm_trans_pblas(psi_in_bfirst, *pX, *psi_ks, *pc, naos, nocc, nvirt, *pmat);
            if (this->tdm_sym) for (auto& t : dm_trans_2d) LR_Util::matsym(t.data<T>(), naos, *pmat);
#else
            std::vector<container::Tensor>  dm_trans_2d = cal_dm_trans_blas(psi_in_bfirst, *psi_ks, nocc, nvirt);
            if (this->tdm_sym) for (auto& t : dm_trans_2d) LR_Util::matsym(t.data<T>(), naos);
#endif
            // tensor to vector, then set DMK
            for (int ik = 0;ik < nk;++ik) { this->DM_trans[ib_dm]->set_DMK_pointer(ik, dm_trans_2d[ik].data<T>()); }

            // if (this->first_print)
            //     for (int ik = 0;ik < nk;++ik)
            //         LR_Util::print_tensor<std::complex<double>>(dm_trans_2d[ik], "1.DMK[ik=" + std::to_string(ik) + "]", this->pmat);

            // use cal_DMR to get DMR form DMK by FT
            this->DM_trans[ib_dm]->cal_DMR();  //DM_trans->get_DMR_vector() is 2d-block parallized
            // LR_Util::print_DMR(*this->DM_trans[0], ucell.nat, "DM(R) (complex)");

            // ========================= begin grid calculation=========================
            this->grid_calculation(nbands, ib_dm);   //DM(R) to H(R)
            // ========================= end grid calculation =========================

            // V(R)->V(k)
            std::vector<ct::Tensor> v_hxc_2d(this->kv.get_nks(),
                ct::Tensor(ct::DataTypeToEnum<T>::value, ct::DeviceTypeToEnum<base_device::DEVICE_CPU>::value,
                    { pmat->get_col_size(), pmat->get_row_size() }));
            for (auto& v : v_hxc_2d) v.zero();
            int nrow = ModuleBase::GlobalFunc::IS_COLUMN_MAJOR_KS_SOLVER(PARAM.inp.ks_solver) ? this->pmat->get_row_size() : this->pmat->get_col_size();
            for (int ik = 0;ik < nk;++ik) { folding_HR(*this->hR, v_hxc_2d[ik].data<T>(), this->kv.kvec_d[ik], nrow, 1); }  // V(R) -> V(k)
            // LR_Util::print_HR(*this->hR, this->ucell.nat, "4.VR");
            // if (this->first_print)
            //     for (int ik = 0;ik < nk;++ik)
            //         LR_Util::print_tensor<T>(v_hxc_2d[ik], "4.V(k)[ik=" + std::to_string(ik) + "]", this->pmat);

            // 5. [AX]^{Hxc}_{ai}=\sum_{\mu,\nu}c^*_{a,\mu,}V^{Hxc}_{\mu,\nu}c_{\nu,i}
#ifdef __MPI
            cal_AX_pblas(v_hxc_2d, *this->pmat, *this->psi_ks, *this->pc, naos, nocc, nvirt, *this->pX, psi_out_bfirst);
#else
            cal_AX_blas(v_hxc_2d, *this->psi_ks, nocc, nvirt, psi_out_bfirst);
#endif
            // if (this->first_print) LR_Util::print_psi_bandfirst(psi_out_bfirst, "5.AX", ib);
        }
    }


    template<>
    void OperatorLRHxc<double, base_device::DEVICE_CPU>::grid_calculation(const int& nbands, const int& iband_dm) const
    {
        ModuleBase::TITLE("OperatorLRHxc", "grid_calculation(real)");
        ModuleBase::timer::tick("OperatorLRHxc", "grid_calculation");
        this->gint->transfer_DM2DtoGrid(this->DM_trans[iband_dm]->get_DMR_vector());     // 2d block to grid

        // 2. transition electron density
        // \f[ \tilde{\rho}(r)=\sum_{\mu_j, \mu_b}\tilde{\rho}_{\mu_j,\mu_b}\phi_{\mu_b}(r)\phi_{\mu_j}(r) \f]
        double** rho_trans;
        const int& nrxx = this->pot.lock()->nrxx;
        // LR_Util::new_p2(rho_trans, nspin_solve, nrxx);
        LR_Util::new_p2(rho_trans, 1, nrxx); // currently gint_kernel_rho uses PARAM.inp.nspin, it needs refactor
        for (int is = 0;is < nspin_solve;++is)ModuleBase::GlobalFunc::ZEROS(rho_trans[is], nrxx);
        Gint_inout inout_rho(rho_trans, Gint_Tools::job_type::rho, 1, false);
        this->gint->cal_gint(&inout_rho);

        // 3. v_hxc = f_hxc * rho_trans
        ModuleBase::matrix vr_hxc(nspin_solve, nrxx);   //grid
        this->pot.lock()->cal_v_eff(rho_trans, &GlobalC::ucell, vr_hxc);
        LR_Util::delete_p2(rho_trans, nspin_solve);

        // 4. V^{Hxc}_{\mu,\nu}=\int{dr} \phi_\mu(r) v_{Hxc}(r) \phi_\mu(r)
        // V(R) for each spin
        for (int is = 0;is < nspin_solve;++is)
        {
            double* vr_hxc_is = &vr_hxc.c[is * nrxx];   //v(r) at current spin
            Gint_inout inout_vlocal(vr_hxc_is, is, Gint_Tools::job_type::vlocal);
            this->gint->get_hRGint()->set_zero();
            this->gint->cal_gint(&inout_vlocal);
        }
        this->hR->set_zero();   // clear hR for each bands
        this->gint->transfer_pvpR(&*this->hR, &GlobalC::ucell);    //grid to 2d block
        ModuleBase::timer::tick("OperatorLRHxc", "grid_calculation");
    }

    template<>
    void OperatorLRHxc<std::complex<double>, base_device::DEVICE_CPU>::grid_calculation(const int& nbands, const int& iband_dm) const
    {
        ModuleBase::TITLE("OperatorLRHxc", "grid_calculation(complex)");
        ModuleBase::timer::tick("OperatorLRHxc", "grid_calculation");

        elecstate::DensityMatrix<std::complex<double>, double> DM_trans_real_imag(pmat, 1, kv.kvec_d, kv.get_nks() / nspin);
        DM_trans_real_imag.init_DMR(*this->hR);
        hamilt::HContainer<double> HR_real_imag(GlobalC::ucell, this->pmat);
        this->initialize_HR(HR_real_imag, ucell, gd, this->pmat);

        auto dmR_to_hR = [&, this](const int& iband_dm, const char& type) -> void
            {
                LR_Util::get_DMR_real_imag_part(*this->DM_trans[iband_dm], DM_trans_real_imag, ucell.nat, type);
                // if (this->first_print)LR_Util::print_DMR(DM_trans_real_imag, ucell.nat, "DMR(2d, real)");

                this->gint->transfer_DM2DtoGrid(DM_trans_real_imag.get_DMR_vector());
                // LR_Util::print_HR(*this->gint->get_DMRGint()[0], this->ucell.nat, "DMR(grid, real)");

                // 2. transition electron density
                double** rho_trans;
                const int& nrxx = this->pot.lock()->nrxx;
                // LR_Util::new_p2(rho_trans, nspin_solve, nrxx);
                LR_Util::new_p2(rho_trans, nspin, nrxx); // currently gint_kernel_rho uses PARAM.inp.nspin, it needs refactor
                for (int is = 0;is < nspin_solve;++is)ModuleBase::GlobalFunc::ZEROS(rho_trans[is], nrxx);
                Gint_inout inout_rho(rho_trans, Gint_Tools::job_type::rho, 1, false);
                this->gint->cal_gint(&inout_rho);
                // print_grid_nonzero(rho_trans[0], nrxx, 10, "rho_trans");

                // 3. v_hxc = f_hxc * rho_trans
                ModuleBase::matrix vr_hxc(nspin_solve, nrxx);   //grid
                this->pot.lock()->cal_v_eff(rho_trans, &GlobalC::ucell, vr_hxc);
                // print_grid_nonzero(vr_hxc.c, this->poticab->nrxx, 10, "vr_hxc");

                LR_Util::delete_p2(rho_trans, nspin_solve);

                // 4. V^{Hxc}_{\mu,\nu}=\int{dr} \phi_\mu(r) v_{Hxc}(r) \phi_\mu(r)
                for (int is = 0;is < nspin_solve;++is)
                {
                    double* vr_hxc_is = &vr_hxc.c[is * nrxx];   //v(r) at current spin
                    Gint_inout inout_vlocal(vr_hxc_is, is, Gint_Tools::job_type::vlocal);
                    this->gint->get_hRGint()->set_zero();
                    this->gint->cal_gint(&inout_vlocal);
                }
                // LR_Util::print_HR(*this->gint->get_hRGint(), this->ucell.nat, "VR(grid)");
                HR_real_imag.set_zero();
                this->gint->transfer_pvpR(&HR_real_imag, &GlobalC::ucell, &GlobalC::GridD);
                // LR_Util::print_HR(HR_real_imag, this->ucell.nat, "VR(real, 2d)");
                LR_Util::set_HR_real_imag_part(HR_real_imag, *this->hR, GlobalC::ucell.nat, type);
            };
        this->hR->set_zero();
        dmR_to_hR(iband_dm, 'R');   //real
        if (kv.get_nks() / this->nspin > 1) { dmR_to_hR(iband_dm, 'I'); }   //imag for multi-k
        ModuleBase::timer::tick("OperatorLRHxc", "grid_calculation");
    }

    template class OperatorLRHxc<double>;
    template class OperatorLRHxc<std::complex<double>>;
}