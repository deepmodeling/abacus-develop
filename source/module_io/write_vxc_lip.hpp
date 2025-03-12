#ifndef __WRITE_VXC_LIP_H_
#define __WRITE_VXC_LIP_H_
#include "module_parameter/parameter.h"
#include "module_base/parallel_reduce.h"
#include "module_base/module_container/base/third_party/blas.h"
#include "module_hamilt_pw/hamilt_pwdft/operator_pw/veff_pw.h"
#include "module_psi/psi.h"
#include "module_cell/unitcell.h"
#include "module_cell/klist.h"
#include "module_elecstate/potentials/potential_new.h"
#include "module_io/write_HS.h"
#include <type_traits>

namespace ModuleIO
{
    template<typename T>
    using Real = typename GetTypeReal<T>::type;

    template<typename FPTYPE>
    inline FPTYPE get_real(const std::complex<FPTYPE>& c)
    {
        return c.real();
    }
    template<typename FPTYPE>
    inline FPTYPE get_real(const FPTYPE& d)
    {
        return d;
    }

    template <typename T>
    std::vector<T> cVc(T* const V, T* const c, int nbasis, int nbands)
    {
        std::vector<T> Vc(nbasis * nbands, 0.0);
        char transa = 'N';
        char transb = 'N';
        const T alpha(1.0, 0.0);
        const T beta(0.0, 0.0);
        container::BlasConnector::gemm(transa, transb, nbasis, nbands, nbasis, alpha, V, nbasis, c, nbasis, beta, Vc.data(), nbasis);

        std::vector<T> cVc(nbands * nbands, 0.0);
        transa = ((std::is_same<T, double>::value || std::is_same<T, float>::value) ? 'T' : 'C');
        container::BlasConnector::gemm(transa, transb, nbands, nbands, nbasis, alpha, c, nbasis, Vc.data(), nbasis, beta, cVc.data(), nbands);
        return cVc;
    }
    template <typename FPTYPE>
    std::vector<std::complex<FPTYPE>> psi_Hpsi(std::complex<FPTYPE>* const psi, std::complex<FPTYPE>* const hpsi, const int nbasis, const int nbands)
    {
        using T = std::complex<FPTYPE>;
        std::vector<T> cVc(nbands * nbands, (T)0.0);
        const T alpha(1.0, 0.0);
        const T beta(0.0, 0.0);
        container::BlasConnector::gemm('C', 'N', nbands, nbands, nbasis, alpha, psi, nbasis, hpsi, nbasis, beta, cVc.data(), nbands);
        return cVc;
    }
    template <typename FPTYPE>
    std::vector<FPTYPE> orbital_energy(const int ik, const int nbands, const std::vector<std::complex<FPTYPE>>& mat_mo)
    {
#ifdef __DEBUG
        assert(nbands >= 0);
#endif
        std::vector<FPTYPE> e(nbands, 0.0);
        for (int i = 0; i < nbands; ++i) {
            e[i] = get_real(mat_mo[i * nbands + i]);
}
        return e;
    }
    template <typename FPTYPE>
    FPTYPE all_band_energy(const int ik, const int nbands, const std::vector<std::complex<FPTYPE>>& mat_mo, const ModuleBase::matrix& wg)
    {
        FPTYPE e = 0.0;
        for (int i = 0; i < nbands; ++i) {
            e += get_real(mat_mo[i * nbands + i]) * (FPTYPE)wg(ik, i);
}
        return e;
    }
    template <typename FPTYPE>
    FPTYPE all_band_energy(const int ik, const std::vector<FPTYPE>& orbital_energy, const ModuleBase::matrix& wg)
    {
        FPTYPE e = 0.0;
        for (int i = 0; i < orbital_energy.size(); ++i) {
            e += orbital_energy[i] * (FPTYPE)wg(ik, i);
}
        return e;
    }

    /// @brief  write the Vxc matrix in KS orbital representation, usefull for GW calculation
    /// including terms: local/semi-local XC and EXX
    template <typename FPTYPE>
    void write_Vxc(int nspin,
                   int naos,
                   int drank,
                   const psi::Psi<std::complex<FPTYPE>>& psi_pw,
                   // const psi::Psi<T>& psi_lcao,
                   const UnitCell& ucell,
                   Structure_Factor& sf,
                   surchem& solvent,
                   const ModulePW::PW_Basis_K& wfc_basis,
                   const ModulePW::PW_Basis& rho_basis,
                   const ModulePW::PW_Basis& rhod_basis,
                   const ModuleBase::matrix& vloc,
                   const Charge& chg,
                   const K_Vectors& kv,
                   const ModuleBase::matrix& wg
#ifdef __EXX
                   ,
                   const Exx_Lip<std::complex<FPTYPE>>& exx_lip
#endif
    )
    {
        using T = std::complex<FPTYPE>;
        ModuleBase::TITLE("ModuleIO", "write_Vxc_LIP");
        int nbands = wg.nc;
        // 1. real-space xc potential
        // ModuleBase::matrix vr_xc(nspin, chg.nrxx);
        double etxc = 0.0;
        double vtxc = 0.0;
        // elecstate::PotXC* potxc(&rho_basis, &etxc, vtxc, nullptr);
        // potxc.cal_v_eff(&chg, &ucell, vr_xc);
        elecstate::Potential* potxc
            = new elecstate::Potential(&rhod_basis, &rho_basis, &ucell, &vloc, &sf, &solvent, &etxc, &vtxc);
        std::vector<std::string> compnents_list = { "xc" };

        potxc->pot_register(compnents_list);
        potxc->update_from_charge(&chg, &ucell);
        // const ModuleBase::matrix vr_localxc = potxc->get_veff_smooth();

        // 2. allocate xc operator
        psi::Psi<T> hpsi_localxc(psi_pw.get_nk(), psi_pw.get_nbands(), psi_pw.get_nbasis(), kv.ngk, true);
        hpsi_localxc.zero_out();
        // std::cout << "hpsi.nk=" << hpsi_localxc.get_nk() << std::endl;
        // std::cout << "hpsi.nbands=" << hpsi_localxc.get_nbands() << std::endl;
        // std::cout << "hpsi.nbasis=" << hpsi_localxc.get_nbasis() << std::endl;
        hamilt::Veff<hamilt::OperatorPW<T>>* vxcs_op_pw;


        std::vector<std::vector<FPTYPE>> e_orb_locxc; // orbital energy (local XC)
        std::vector<std::vector<FPTYPE>> e_orb_tot;   // orbital energy (total)
        std::vector<std::vector<FPTYPE>> e_orb_exx; // orbital energy (EXX)
        Parallel_2D p2d_serial;
        p2d_serial.set_serial(nbands, nbands);
        // ======test=======
        std::vector<FPTYPE> exx_energy(kv.get_nks());
        // ======test=======s
        for (int ik = 0; ik < kv.get_nks(); ++ik)
        {
            // 2.1 local xc
            vxcs_op_pw = new hamilt::Veff<hamilt::OperatorPW<T>>(kv.isk.data(),
                potxc->get_veff_smooth_data<FPTYPE>(), potxc->get_veff_smooth().nr, potxc->get_veff_smooth().nc, &wfc_basis);
            vxcs_op_pw->init(ik);   // set k-point index
            psi_pw.fix_k(ik);
            hpsi_localxc.fix_k(ik);
#ifdef __DEBUG
            assert(hpsi_localxc.get_current_nbas() == psi_pw.get_current_nbas());
            assert(hpsi_localxc.get_current_nbas() == hpsi_localxc.get_ngk(ik));
#endif
            /// wrap psi and act band-by-band (the same result as act all bands at once)
            // for (int ib = 0;ib < psi_pw.get_nbands();++ib)
            // {
            //     std::cout<<"ib="<<ib<<std::endl;
            //     psi::Psi<T> psi_single_band(&psi_pw(ik, ib, 0), 1, 1, psi_pw.get_current_nbas());
            //     psi::Psi<T> hpsi_single_band(&hpsi_localxc(ik, ib, 0), 1, 1, hpsi_localxc.get_current_nbas());
            //     vxcs_op_pw->act(1, psi_pw.get_current_nbas(), psi_pw.npol, psi_single_band.get_pointer(), hpsi_single_band.get_pointer(), psi_pw.get_ngk(ik));
            // }
            vxcs_op_pw->act(psi_pw.get_nbands(), psi_pw.get_nbasis(), psi_pw.get_npol(), &psi_pw(ik, 0, 0), &hpsi_localxc(ik, 0, 0), psi_pw.get_ngk(ik));
            delete vxcs_op_pw;
            std::vector<T> vxc_local_k_mo = psi_Hpsi(&psi_pw(ik, 0, 0), &hpsi_localxc(ik, 0, 0), psi_pw.get_nbasis(), psi_pw.get_nbands());
            Parallel_Reduce::reduce_pool(vxc_local_k_mo.data(), nbands * nbands);
            e_orb_locxc.emplace_back(orbital_energy(ik, nbands, vxc_local_k_mo));

            // 2.2 exx
            std::vector<T> vxc_tot_k_mo(std::move(vxc_local_k_mo));
            std::vector<T> vexx_k_ao(naos * naos);
#if((defined __LCAO)&&(defined __EXX) && !(defined __CUDA)&& !(defined __ROCM))
            if (GlobalC::exx_info.info_global.cal_exx)
            {
                for (int n = 0; n < naos; ++n) {
                    for (int m = 0; m < naos; ++m) {
                        vexx_k_ao[n * naos + m] += (T)GlobalC::exx_info.info_global.hybrid_alpha * exx_lip.get_exx_matrix()[ik][m][n];
}
}
                std::vector<T> vexx_k_mo = cVc(vexx_k_ao.data(), &(exx_lip.get_hvec()(ik, 0, 0)), naos, nbands);
                Parallel_Reduce::reduce_pool(vexx_k_mo.data(), nbands * nbands);
                e_orb_exx.emplace_back(orbital_energy(ik, nbands, vexx_k_mo));
                // ======test=======    
                // std::cout << "exx_energy from matrix:" << all_band_energy(ik, nbands, vexx_k_mo, wg) << std::endl;
                // std::cout << "exx_energy from orbitals: " << all_band_energy(ik, e_orb_exx.at(ik), wg) << std::endl;
                // std::cout << "exx_energy from exx_lip: " << GlobalC::exx_info.info_global.hybrid_alpha * exx_lip.get_exx_energy() << std::endl;
                // ======test=======
                container::BlasConnector::axpy(nbands * nbands, 1.0, vexx_k_mo.data(), 1, vxc_tot_k_mo.data(), 1);
            }
#endif
            /// add-up and write
            ModuleIO::save_mat(-1, vxc_tot_k_mo.data(), nbands, false, PARAM.inp.out_ndigits, true, false, "Vxc", "k-" + std::to_string(ik), p2d_serial, drank, false);
            e_orb_tot.emplace_back(orbital_energy(ik, nbands, vxc_tot_k_mo));
        }
        //===== test total xc energy =======
        // std::cout << "xc energy =" << etxc << std::endl;
        // std::cout << "vtxc=" << vtxc << std::endl;
        // FPTYPE exc_by_orb = 0.0;
        // for (int ik = 0;ik < e_orb_locxc.size();++ik)
        //     exc_by_orb += all_band_energy(ik, e_orb_locxc[ik], wg);
        // std::cout << "xc all-bands energy by orbital =" << exc_by_orb << std::endl;
        // /// calculate orbital energy by grid integration of vtxc*rho
        // FPTYPE exc_by_rho = 0.0;
        // for (int ir = 0;ir < potxc->get_veff_smooth().nc;++ir)
        //     exc_by_rho += potxc->get_veff_smooth()(0, ir) * chg.rho[0][ir];
        // Parallel_Reduce::reduce_all(exc_by_rho);
        // exc_by_rho *= ((FPTYPE)ucell.omega * (FPTYPE)GlobalV::NPROC / (FPTYPE)potxc->get_veff_smooth().nc);
        // std::cout << "xc all-bands energy by rho =" << exc_by_rho << std::endl;
        //===== test total xc energy =======
        //===== test total exx energy =======
// #if((defined __LCAO)&&(defined __EXX) && !(defined __CUDA)&& !(defined __ROCM))
//         if (GlobalC::exx_info.info_global.cal_exx)
//         {
//             FPTYPE exx_by_orb = 0.0;
//             for (int ik = 0;ik < e_orb_exx.size();++ik)
//                 exx_by_orb += all_band_energy(ik, e_orb_exx[ik], wg);
//             exx_by_orb /= 2;
//             std::cout << "exx all-bands energy by orbital =" << exx_by_orb << std::endl;
//             FPTYPE exx_from_lip = GlobalC::exx_info.info_global.hybrid_alpha * exx_lip.get_exx_energy();
//             std::cout << "exx all-bands energy from exx_lip =" << exx_from_lip << std::endl;
//         }
// #endif
        //===== test total exx energy =======
        // write the orbital energy for xc and exx in LibRPA format
        const int nspin0 = (nspin == 2) ? 2 : 1;
        auto write_orb_energy = [&kv, &nspin0, &nbands](const std::vector<std::vector<FPTYPE>>& e_orb,
            const std::string& label,
            const bool app = false) {
                assert(e_orb.size() == kv.get_nks());
                const int nk = kv.get_nks() / nspin0;
                std::ofstream ofs;
                ofs.open(PARAM.globalv.global_out_dir + "vxc_" + (label == "" ? "out.dat" : label + "_out.dat"),
                    app ? std::ios::app : std::ios::out);
                ofs << nk << "\n" << nspin0 << "\n" << nbands << "\n";
                ofs << std::scientific << std::setprecision(16);
                for (int ik = 0; ik < nk; ++ik)
                {
                    for (int is = 0; is < nspin0; ++is)
                    {
                        for (auto e : e_orb[is * nk + ik])
                        { // Hartree and eV
                            ofs << e / 2. << "\t" << e * (FPTYPE)ModuleBase::Ry_to_eV << "\n";
                        }
                    }
                }
            };
        if (GlobalV::MY_RANK == 0)
        {
            write_orb_energy(e_orb_tot, "");
#if((defined __LCAO)&&(defined __EXX) && !(defined __CUDA)&& !(defined __ROCM))
            if (GlobalC::exx_info.info_global.cal_exx)
            {
                write_orb_energy(e_orb_locxc, "local");
                write_orb_energy(e_orb_exx, "exx");
            }
#endif
        }
    }
} // namespace ModuleIO
#endif