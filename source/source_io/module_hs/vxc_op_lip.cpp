#include "vxc_op_lip.h"

#include "source_base/module_out/filename.h"
#include "source_base/parallel_reduce.h"
#include "source_base/module_container/base/third_party/blas.h"
#include "source_io/module_hs/hs_dense_io.h"

#include <fstream>
#include <iomanip>
#include <memory>

namespace ModuleIO
{

namespace
{

template <typename FPTYPE>
FPTYPE get_real_lip(const std::complex<FPTYPE>& c)
{
    return c.real();
}

template <typename FPTYPE>
FPTYPE get_real_lip(const FPTYPE& d)
{
    return d;
}

template <typename T>
std::vector<T> cVc_lip(const T* V, const T* c, int nbasis, int nbands)
{
    std::vector<T> Vc(nbasis * nbands, 0.0);
    char transa = 'N';
    char transb = 'N';
    const T alpha(1.0, 0.0);
    const T beta(0.0, 0.0);
    container::BlasConnector::gemm(transa, transb, nbasis, nbands, nbasis,
        alpha, V, nbasis, c, nbasis, beta, Vc.data(), nbasis);

    std::vector<T> cVc(nbands * nbands, 0.0);
    transa = ((std::is_same<T, double>::value || std::is_same<T, float>::value) ? 'T' : 'C');
    container::BlasConnector::gemm(transa, transb, nbands, nbands, nbasis,
        alpha, c, nbasis, Vc.data(), nbasis, beta, cVc.data(), nbands);
    return cVc;
}

template <typename FPTYPE>
std::vector<std::complex<FPTYPE>> psi_Hpsi(const std::complex<FPTYPE>* psi,
                                           const std::complex<FPTYPE>* hpsi,
                                           const int nbasis,
                                           const int nbands)
{
    using T = std::complex<FPTYPE>;
    std::vector<T> cVc(nbands * nbands, static_cast<T>(0.0));
    const T alpha(1.0, 0.0);
    const T beta(0.0, 0.0);
    container::BlasConnector::gemm('C', 'N', nbands, nbands, nbasis, alpha,
        psi, nbasis, hpsi, nbasis, beta, cVc.data(), nbands);
    return cVc;
}

template <typename FPTYPE>
std::vector<FPTYPE> orbital_energy_lip(const int ik, const int nbands,
                                       const std::vector<std::complex<FPTYPE>>& mat_mo)
{
#ifdef __DEBUG
    assert(nbands >= 0);
#endif
    std::vector<FPTYPE> e(nbands, 0.0);
    for (int i = 0; i < nbands; ++i)
    {
        e[i] = get_real_lip(mat_mo[i * nbands + i]);
    }
    return e;
}

template <typename FPTYPE>
FPTYPE all_band_energy_lip(const int ik, const int nbands,
                           const std::vector<std::complex<FPTYPE>>& mat_mo,
                           const ModuleBase::matrix& wg)
{
    FPTYPE e = 0.0;
    for (int i = 0; i < nbands; ++i)
    {
        e += get_real_lip(mat_mo[i * nbands + i]) * static_cast<FPTYPE>(wg(ik, i));
    }
    return e;
}

template <typename FPTYPE>
FPTYPE all_band_energy_lip(const int ik, const std::vector<FPTYPE>& orbital_energy,
                           const ModuleBase::matrix& wg)
{
    FPTYPE e = 0.0;
    for (size_t i = 0; i < orbital_energy.size(); ++i)
    {
        e += orbital_energy[i] * static_cast<FPTYPE>(wg(ik, i));
    }
    return e;
}

} // anonymous namespace

template <typename FPTYPE>
void write_Vxc_LIP(int nspin,
                   int naos,
                   int drank,
                   const psi::Psi<std::complex<FPTYPE>>& psi_pw,
                   const UnitCell& ucell,
                   Structure_Factor& sf,
                   surchem& solvent,
                   const ModulePW::PW_Basis_K& wfc_basis,
                   const ModulePW::PW_Basis& rho_basis,
                   const ModulePW::PW_Basis& rhod_basis,
                   const ModuleBase::matrix& vloc,
                   const Charge& chg,
                   const K_Vectors& kv,
                   const ModuleBase::matrix& wg,
                   const bool gamma_only,
                   const std::string& global_out_dir,
                   const int out_ndigits,
                   const std::string& ks_solver,
                   bool cal_exx,
                   double hybrid_alpha
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
    double etxc = 0.0;
    double vtxc = 0.0;
    std::unique_ptr<elecstate::Potential> potxc(
        new elecstate::Potential(&rhod_basis, &rho_basis, &ucell, &vloc, &sf, &solvent, &etxc, &vtxc));
    std::vector<std::string> compnents_list = {"xc"};

    potxc->pot_register(compnents_list);
    potxc->update_from_charge(&chg, &ucell);

    // 2. allocate xc operator
    psi::Psi<T> hpsi_localxc(psi_pw.get_nk(), psi_pw.get_nbands(), psi_pw.get_nbasis(), kv.ngk, true);
    hpsi_localxc.zero_out();
    std::unique_ptr<hamilt::Veff<hamilt::OperatorPW<T>>> vxcs_op_pw;

    std::vector<std::vector<FPTYPE>> e_orb_locxc; // orbital energy (local XC)
    std::vector<std::vector<FPTYPE>> e_orb_tot;   // orbital energy (total)
    std::vector<std::vector<FPTYPE>> e_orb_exx;   // orbital energy (EXX)
    Parallel_2D p2d_serial;
    p2d_serial.set_serial(nbands, nbands);

    for (int ik = 0; ik < kv.get_nks(); ++ik)
    {
        // 2.1 local xc
        vxcs_op_pw = std::unique_ptr<hamilt::Veff<hamilt::OperatorPW<T>>>(
            new hamilt::Veff<hamilt::OperatorPW<T>>(kv.isk.data(),
                potxc->get_veff_smooth_data<FPTYPE>(), potxc->get_veff_smooth().nr, potxc->get_veff_smooth().nc, &wfc_basis));
        vxcs_op_pw->init(ik);   // set k-point index
        psi_pw.fix_k(ik);
        hpsi_localxc.fix_k(ik);
#ifdef __DEBUG
        assert(hpsi_localxc.get_current_nbas() == psi_pw.get_current_nbas());
        assert(hpsi_localxc.get_current_nbas() == hpsi_localxc.get_ngk(ik));
#endif
        vxcs_op_pw->act(psi_pw.get_nbands(), psi_pw.get_nbasis(), psi_pw.get_npol(), &psi_pw(ik, 0, 0), &hpsi_localxc(ik, 0, 0), psi_pw.get_ngk(ik));
        vxcs_op_pw.reset();
        std::vector<T> vxc_local_k_mo = psi_Hpsi(&psi_pw(ik, 0, 0), &hpsi_localxc(ik, 0, 0), psi_pw.get_nbasis(), psi_pw.get_nbands());
        Parallel_Reduce::reduce_pool(vxc_local_k_mo.data(), nbands * nbands);
        e_orb_locxc.emplace_back(orbital_energy_lip(ik, nbands, vxc_local_k_mo));

        // 2.2 exx
        std::vector<T> vxc_tot_k_mo(std::move(vxc_local_k_mo));
        std::vector<T> vexx_k_ao(naos * naos);
#if((defined __LCAO)&&(defined __EXX) && !(defined __CUDA)&& !(defined __ROCM))
        if (cal_exx)
        {
            for (int n = 0; n < naos; ++n)
            {
                for (int m = 0; m < naos; ++m)
                {
                    vexx_k_ao[n * naos + m] += static_cast<T>(hybrid_alpha)
                        * exx_lip.get_exx_matrix()[ik][m][n];
                }
            }
            std::vector<T> vexx_k_mo = cVc_lip(vexx_k_ao.data(), &(exx_lip.get_hvec()(ik, 0, 0)), naos, nbands);
            Parallel_Reduce::reduce_pool(vexx_k_mo.data(), nbands * nbands);
            e_orb_exx.emplace_back(orbital_energy_lip(ik, nbands, vexx_k_mo));
            container::BlasConnector::axpy(nbands * nbands, 1.0, vexx_k_mo.data(), 1, vxc_tot_k_mo.data(), 1);
        }
#endif

        // add-up and write
        const int istep = -1;
        const int out_label = 1; // 1 means .txt while 2 means .dat
        const bool out_app_flag = 0;
        std::string vxc_file = ModuleIO::filename_output(
            global_out_dir,
            "vxc", "nao", ik, kv.ik2iktot, nspin, kv.get_nkstot(),
            out_label, out_app_flag, gamma_only, istep);

        ModuleIO::save_mat(istep, vxc_tot_k_mo.data(), nbands,
            false, out_ndigits, true,
            out_app_flag, vxc_file,
            p2d_serial, drank, ks_solver, false);

        e_orb_tot.emplace_back(orbital_energy_lip(ik, nbands, vxc_tot_k_mo));
    }

    // write the orbital energy for xc and exx in LibRPA format
    const int nspin0 = (nspin == 2) ? 2 : 1;
    auto write_orb_energy_lip = [&kv, &nspin0, &nbands, &global_out_dir](const std::vector<std::vector<FPTYPE>>& e_orb,
        const std::string& label,
        const bool app = false) {
            assert(e_orb.size() == kv.get_nks());
            const int nk = kv.get_nks() / nspin0;
            std::ofstream ofs;
            const std::string out_name = (label == "") ? ".txt" : "_" + label + ".txt";
            ofs.open(global_out_dir + "vxc" + out_name,
                app ? std::ios::app : std::ios::out);
            ofs << nk << "\n" << nspin0 << "\n" << nbands << "\n";
            ofs << std::scientific << std::setprecision(16);
            for (int ik = 0; ik < nk; ++ik)
            {
                for (int is = 0; is < nspin0; ++is)
                {
                    for (auto e : e_orb[is * nk + ik])
                    { // Hartree and eV
                        ofs << e / 2. << "\t" << e * static_cast<FPTYPE>(ModuleBase::Ry_to_eV) << "\n";
                    }
                }
            }
        };

    if (GlobalV::MY_RANK == 0)
    {
        write_orb_energy_lip(e_orb_tot, "");
#if((defined __LCAO)&&(defined __EXX) && !(defined __CUDA)&& !(defined __ROCM))
        if (cal_exx)
        {
            write_orb_energy_lip(e_orb_locxc, "local");
            write_orb_energy_lip(e_orb_exx, "exx");
        }
#endif
    }
}

// Explicit template instantiations
template void write_Vxc_LIP<float>(
    int, int, int, const psi::Psi<std::complex<float>>&, const UnitCell&, Structure_Factor&, surchem&,
    const ModulePW::PW_Basis_K&, const ModulePW::PW_Basis&, const ModulePW::PW_Basis&,
    const ModuleBase::matrix&, const Charge&, const K_Vectors&, const ModuleBase::matrix&,
    const bool, const std::string&, const int, const std::string&, bool, double
#ifdef __EXX
    , const Exx_Lip<std::complex<float>>&
#endif
);

template void write_Vxc_LIP<double>(
    int, int, int, const psi::Psi<std::complex<double>>&, const UnitCell&, Structure_Factor&, surchem&,
    const ModulePW::PW_Basis_K&, const ModulePW::PW_Basis&, const ModulePW::PW_Basis&,
    const ModuleBase::matrix&, const Charge&, const K_Vectors&, const ModuleBase::matrix&,
    const bool, const std::string&, const int, const std::string&, bool, double
#ifdef __EXX
    , const Exx_Lip<std::complex<double>>&
#endif
);

} // namespace ModuleIO
