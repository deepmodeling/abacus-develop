#include "vxc_op_mat.h"

#include "source_base/constants.h"
#include "source_base/module_out/filename.h"
#include "source_base/parallel_reduce.h"
#include "source_io/module_hs/hs_dense_io.h"
#include "source_io/module_output/frame_append.h"
#include "source_io/module_output/provenance_header.h"
#include "source_io/module_output/write_notice.h"

namespace ModuleIO
{

template <typename TK, typename TR>
void write_Vxc(const int nspin,
               const int istep,
               const double etot,
               const int nbasis,
               const int drank,
               const Parallel_Orbitals* pv,
               const psi::Psi<TK>& psi,
               const UnitCell& ucell,
               Structure_Factor& sf,
               surchem& solvent,
               const ModulePW::PW_Basis& rho_basis,
               const ModulePW::PW_Basis& rhod_basis,
               const ModuleBase::matrix& vloc,
               const Charge& chg,
               const K_Vectors& kv,
               const std::vector<double>& orb_cutoff,
               const ModuleBase::matrix& wg,
               Grid_Driver& gd,
               const bool dft_plus_u,
               const bool gamma_only,
               const bool out_mat,
               const int out_app_flag,
               const std::string& mat_dir,
               const std::string& global_out_dir,
               const int out_ndigits,
               const std::string& ks_solver,
               bool cal_exx,
               const Exx_Info& exx_info,
               std::ofstream& ofs_running
#ifdef __EXX
               ,
               std::vector<std::map<int, std::map<hamilt::TAC, RI::Tensor<double>>>>* Hexxd,
               std::vector<std::map<int, std::map<hamilt::TAC, RI::Tensor<std::complex<double>>>>>* Hexxc
#endif
)
{
    ModuleBase::TITLE("ModuleIO", "write_Vxc");
    int nbands = wg.nc;

    // 1. real-space xc potential
    double etxc = 0.0;
    double vtxc = 0.0;
    std::unique_ptr<elecstate::Potential> potxc(
        new elecstate::Potential(&rhod_basis, &rho_basis, &ucell, &vloc, &sf, &solvent, &etxc, &vtxc));
    std::vector<std::string> compnents_list = {"xc"};
    potxc->pot_register(compnents_list);
    potxc->update_from_charge(&chg, &ucell);

    // 2. allocate AO-matrix
    // R (the number of hR: 1 for nspin=1, 4; 2 for nspin=2)
    int nspin0 = (nspin == 2) ? 2 : 1;
    std::vector<hamilt::HContainer<TR>> vxcs_R_ao(nspin0, hamilt::HContainer<TR>(ucell, pv));
    for (int is = 0; is < nspin0; ++is) {
        vxcs_R_ao[is].set_zero();
        if (std::is_same<TK, double>::value) { vxcs_R_ao[is].fix_gamma(); }
    }
    // k (size for each k-point)
    hamilt::HS_Matrix_K<TK> vxc_k_ao(pv, 1); // only hk is needed, sk is skipped

    // 3. allocate operators and contribute HR
    std::vector<std::unique_ptr<hamilt::Veff<hamilt::OperatorLCAO<TK, TR>>>> vxcs_op_ao(nspin0);
    for (int is = 0; is < nspin0; ++is)
    {
        vxcs_op_ao[is] = std::unique_ptr<hamilt::Veff<hamilt::OperatorLCAO<TK, TR>>>(
            new hamilt::Veff<hamilt::OperatorLCAO<TK, TR>>(
                &vxc_k_ao, kv.kvec_d, potxc.get(), &vxcs_R_ao[is], &ucell, orb_cutoff, &gd, nspin));
        vxcs_op_ao[is]->set_current_spin(is);
        vxcs_op_ao[is]->contributeHR();
    }
    std::vector<std::vector<double>> e_orb_locxc; // orbital energy (local XC)
    std::vector<std::vector<double>> e_orb_tot;   // orbital energy (total)
#ifdef __EXX
    hamilt::OperatorEXX<hamilt::OperatorLCAO<TK, TR>> vexx_op_ao(&vxc_k_ao,
        &vxcs_R_ao[0], ucell, kv, Hexxd, Hexxc, &exx_info, hamilt::Add_Hexx_Type::k);
    hamilt::HS_Matrix_K<TK> vexxonly_k_ao(pv, 1); // only hk is needed, sk is skipped
    hamilt::OperatorEXX<hamilt::OperatorLCAO<TK, TR>> vexxonly_op_ao(&vexxonly_k_ao,
        &vxcs_R_ao[0], ucell, kv, Hexxd, Hexxc, &exx_info, hamilt::Add_Hexx_Type::k);
    std::vector<std::vector<double>> e_orb_exx; // orbital energy (EXX)
#endif
    hamilt::DFTU_firstzeta<hamilt::OperatorLCAO<TK, TR>> vdftu_op_ao(&vxc_k_ao, kv.kvec_d, nullptr, ucell, nullptr, kv.isk);

    // 4. calculate and write the MO-matrix Exc
    Parallel_Blkcyclic p2d;
    set_para2d_MO(*pv, nbands, p2d);

    for (int ik = 0; ik < kv.get_nks(); ++ik)
    {
        vxc_k_ao.set_zero_hk();
        int is = kv.isk[ik];
        dynamic_cast<hamilt::OperatorLCAO<TK, TR>*>(vxcs_op_ao[is].get())->contributeHk(ik);
        const std::vector<TK>& vlocxc_k_mo = cVc(vxc_k_ao.get_hk(), &psi(ik, 0, 0), nbasis, nbands, *pv, p2d);

#ifdef __EXX
        if (cal_exx)
        {
            e_orb_locxc.emplace_back(orbital_energy(ik, nbands, vlocxc_k_mo, p2d));
            ModuleBase::GlobalFunc::ZEROS(vexxonly_k_ao.get_hk(), pv->para_2d().size());
            vexx_op_ao.contributeHk(ik);
            vexxonly_op_ao.contributeHk(ik);
            std::vector<TK> vexx_k_mo = cVc(vexxonly_k_ao.get_hk(), &psi(ik, 0, 0), nbasis, nbands, *pv, p2d);
            e_orb_exx.emplace_back(orbital_energy(ik, nbands, vexx_k_mo, p2d));
        }
#endif
        if (dft_plus_u)
        {
            vdftu_op_ao.contributeHk(ik);
        }
        const std::vector<TK>& vxc_tot_k_mo = cVc(vxc_k_ao.get_hk(), &psi(ik, 0, 0), nbasis, nbands, *pv, p2d);
        e_orb_tot.emplace_back(orbital_energy(ik, nbands, vxc_tot_k_mo, p2d));

        // write the full MO-basis matrix only when requested
        if (out_mat)
        {
            const int out_label = 1; // 1 means .txt while 2 means .dat
            // Mode 1 appends the frames written after the first ionic step;
            // mode 0 truncates to the latest snapshot and mode 2 writes a
            // separate g{istep+1} file for each step.
            const bool append_vxc = frame_appends(out_app_flag, istep);
            std::string vxc_file = ModuleIO::filename_output(
                    mat_dir,
                    "vxcmo", "nao", ik, kv.ik2iktot, nspin, kv.get_nkstot(),
                    out_label, out_app_flag, gamma_only, istep);

            ModuleIO::save_mat(istep,
                               vxc_tot_k_mo.data(),
                               nbands,
                               false /*binary*/,
                               out_ndigits,
                               true /*triangle*/,
                               append_vxc,
                               vxc_file,
                               p2d,
                               drank,
                               ks_solver);

            if (GlobalV::MY_RANK == 0)
            {
                ModuleIO::write_notice(ofs_running, "Vxc(MO)", is, nspin, vxc_file);
            }
        }
    }

    if (GlobalV::MY_RANK == 0)
    {
        const double etot_eV = etot * ModuleBase::Ry_to_eV;
        const std::string time_stamp = local_timestamp();
        // Mode 1 appends the frames written after the first ionic step; modes
        // 0 and 2 truncate (mode 2 writes a separate g{step} file).
        const bool append = frame_appends(out_app_flag, istep);
        const std::string header_tot = build_eband_header(istep, etot_eV, "V_xc", time_stamp);
        const std::string name_tot = eband_filename("eband_vxc", istep, out_app_flag);
        write_orb_energy(kv, nspin0, nbands, e_orb_tot, name_tot, header_tot, global_out_dir, append);
        const std::string property_tot = "Band-decomposed V_xc energy";
        const std::string path_tot = global_out_dir + name_tot;
        ModuleIO::write_notice(ofs_running, property_tot, path_tot);
#ifdef __EXX
        if (cal_exx)
        {
            const std::string header_loc = build_eband_header(istep, etot_eV, "V_xc_local", time_stamp);
            const std::string header_exx = build_eband_header(istep, etot_eV, "V_EXX", time_stamp);
            const std::string name_loc = eband_filename("eband_vxc_local", istep, out_app_flag);
            const std::string name_exx = eband_filename("eband_vxc_exx", istep, out_app_flag);
            write_orb_energy(kv, nspin0, nbands, e_orb_locxc, name_loc, header_loc, global_out_dir, append);
            write_orb_energy(kv, nspin0, nbands, e_orb_exx, name_exx, header_exx, global_out_dir, append);
            const std::string property_loc = "Band-decomposed V_xc_local energy";
            const std::string property_exx = "Band-decomposed V_EXX energy";
            const std::string path_loc = global_out_dir + name_loc;
            const std::string path_exx = global_out_dir + name_exx;
            ModuleIO::write_notice(ofs_running, property_loc, path_loc);
            ModuleIO::write_notice(ofs_running, property_exx, path_exx);
        }
#endif
    }
}

// Explicit template instantiations
template void write_Vxc<double, double>(
    const int, const int, const double, const int, const int, const Parallel_Orbitals*, const psi::Psi<double>&,
    const UnitCell&, Structure_Factor&, surchem&, const ModulePW::PW_Basis&, const ModulePW::PW_Basis&,
    const ModuleBase::matrix&, const Charge&, const K_Vectors&, const std::vector<double>&,
    const ModuleBase::matrix&, Grid_Driver&, const bool, const bool, const bool, const int, const std::string&,
    const std::string&, const int, const std::string&, bool, const Exx_Info&, std::ofstream&
#ifdef __EXX
    , std::vector<std::map<int, std::map<hamilt::TAC, RI::Tensor<double>>>>*,
    std::vector<std::map<int, std::map<hamilt::TAC, RI::Tensor<std::complex<double>>>>>*
#endif
);

template void write_Vxc<std::complex<double>, double>(
    const int, const int, const double, const int, const int, const Parallel_Orbitals*, const psi::Psi<std::complex<double>>&,
    const UnitCell&, Structure_Factor&, surchem&, const ModulePW::PW_Basis&, const ModulePW::PW_Basis&,
    const ModuleBase::matrix&, const Charge&, const K_Vectors&, const std::vector<double>&,
    const ModuleBase::matrix&, Grid_Driver&, const bool, const bool, const bool, const int, const std::string&,
    const std::string&, const int, const std::string&, bool, const Exx_Info&, std::ofstream&
#ifdef __EXX
    , std::vector<std::map<int, std::map<hamilt::TAC, RI::Tensor<double>>>>*,
    std::vector<std::map<int, std::map<hamilt::TAC, RI::Tensor<std::complex<double>>>>>*
#endif
);

template void write_Vxc<std::complex<double>, std::complex<double>>(
    const int, const int, const double, const int, const int, const Parallel_Orbitals*, const psi::Psi<std::complex<double>>&,
    const UnitCell&, Structure_Factor&, surchem&, const ModulePW::PW_Basis&, const ModulePW::PW_Basis&,
    const ModuleBase::matrix&, const Charge&, const K_Vectors&, const std::vector<double>&,
    const ModuleBase::matrix&, Grid_Driver&, const bool, const bool, const bool, const int, const std::string&,
    const std::string&, const int, const std::string&, bool, const Exx_Info&, std::ofstream&
#ifdef __EXX
    , std::vector<std::map<int, std::map<hamilt::TAC, RI::Tensor<double>>>>*,
    std::vector<std::map<int, std::map<hamilt::TAC, RI::Tensor<std::complex<double>>>>>*
#endif
);

} // namespace ModuleIO
