#include "ctrl_scf_lcao.h" // use ctrl_scf_lcao()

#include "source_base/formatter.h"
#include "source_estate/elecstate_lcao.h" // use elecstate::ElecState
#include "source_hamilt/hamilt.h"         // use Hamilt<T>
#include "source_lcao/hamilt_lcao.h"      // use hamilt::HamiltLCAO<TK, TR>

#include <complex>

// functions
#include "../module_unk/berryphase.h"                          // use berryphase
#include "../module_hs/cal_pLpR.h"                            // use AngularMomentumCalculator()
#include "source_io/module_hs/output_mat_sparse.h"                   // use ModuleIO::output_mat_sparse()
#include "source_io/module_ml/io_npz.h"                       // use ModuleIO::output_mat_npz()
#include "source_io/module_dhs/write_dH.h"                    // use ModuleIO::write_dH_components()
#include "source_io/module_hs/write_H_terms.h"         // use ModuleIO::write_h_*
#include "../module_hs/write_HS_R.h"                          // use ModuleIO::write_hsr()
#include "../module_mulliken/cal_mag.h"                          // use cal_mag()
#include "../module_wannier/to_wannier90_lcao.h"                   // use toWannier90_LCAO
#include "../module_wannier/to_wannier90_lcao_in_pw.h"             // use toWannier90_LCAO_IN_PW
#include "../module_hs/write_HS.h"                            // use ModuleIO::write_hsk()
#include "../module_dm/write_dmk.h"                           // use ModuleIO::write_dmk()
#include "../module_dm/write_dmr.h"                           // use ModuleIO::write_dmr()
#include "../module_dos/write_dos_lcao.h"                      // use ModuleIO::write_dos_lcao()
#include "../module_wf/write_wfc_nao.h"                       // use ModuleIO::write_wfc_nao()
#include "source_lcao/module_deltaspin/spin_constrain.h"   // use spinconstrain::SpinConstrain<TK>
#include "source_lcao/module_operator_lcao/ekinetic.h" // use hamilt::EKinetic
#ifdef __MLALGO
#include "source_lcao/module_deepks/LCAO_deepks.h"
#include "source_lcao/module_deepks/LCAO_deepks_interface.h"
#endif
#ifdef __EXX
#include "source_lcao/module_ri/Exx_LRI_interface.h" // use EXX codes
#include "source_lcao/module_ri/RPA_LRI.h"           // use RPA code
#endif
#include "../module_qo/to_qo.h"                // use toQO
#include "source_lcao/module_rdmft/rdmft.h" // use RDMFT codes
#include "source_lcao/rho_tau_lcao.h"       // mohan add 2025-10-24
#include "source_lcao/module_operator_lcao/overlap.h" // use hamilt::Overlap for NAMD

#ifdef __EXX
template <typename TK>
void setup_exx_dh_params(ModuleIO::WriteDHParams& dh_params, Exx_NAO<TK>& exx_nao)
{}

template <>
void setup_exx_dh_params<double>(ModuleIO::WriteDHParams& dh_params, Exx_NAO<double>& exx_nao)
{
    if (GlobalC::exx_info.info_global.cal_exx)
    {
        if (exx_nao.exd) { dh_params.exd = exx_nao.exd.get(); }
        if (exx_nao.exc) { dh_params.exc = exx_nao.exc.get(); }
    }
}

template <typename TK>
void setup_exx_h_params(ModuleIO::WriteHParams& h_params, Exx_NAO<TK>& exx_nao)
{}

template <>
void setup_exx_h_params<double>(ModuleIO::WriteHParams& h_params, Exx_NAO<double>& exx_nao)
{
    if (GlobalC::exx_info.info_global.cal_exx)
    {
        if (exx_nao.exd) { h_params.exd = exx_nao.exd.get(); }
        if (exx_nao.exc) { h_params.exc = exx_nao.exc.get(); }
        ModuleIO::write_h_exx(h_params);
    }
}
#endif

template <typename TK, typename TR>
void ModuleIO::ctrl_scf_lcao(UnitCell& ucell,
                             const Input_para& inp,
                             K_Vectors& kv,
                             elecstate::ElecState* pelec,
                             elecstate::DensityMatrix<TK, double>* dm, // mohan add 2025-11-04
                             Parallel_Orbitals& pv,
                             Grid_Driver& gd,
                             psi::Psi<TK>* psi,
                             hamilt::HamiltLCAO<TK, TR>* p_hamilt,
                             Plus_U& dftu, // mohan add 2025-11-07
                             TwoCenterBundle& two_center_bundle,
                             LCAO_Orbitals& orb,
                             const ModulePW::PW_Basis_K* pw_wfc,   // for berryphase
                             const ModulePW::PW_Basis* pw_rho,     // for berryphase
                             const ModulePW::PW_Basis_Big* pw_big, // for Wannier90
                             const Structure_Factor& sf,           // for Wannier90
                             const ModulePW::PW_Basis* pw_rhod,    // dense charge grid (for dH veff pots)
                             const ModuleBase::matrix& vloc,       // local pseudopotential (for dH veff pots)
                             surchem& solvent,                     // solvent model (for dH veff pots)
                             rdmft::RDMFT<TK, TR>& rdmft_solver,   // for RDMFT
                             Setup_DeePKS<TK>& deepks,
                             Exx_NAO<TK>& exx_nao,
                             const bool conv_esolver,
                             const bool scf_nmax_flag,
                             const int istep)
{
    ModuleBase::TITLE("ModuleIO", "ctrl_scf_lcao");
    ModuleBase::timer::start("ModuleIO", "ctrl_scf_lcao");

    //*****
    // if istep_in = -1, istep will not appear in file name
    // if iter_in = -1, iter will not appear in file name
    int istep_in = -1;
    int iter_in = -1;
    bool out_flag = false;
    if (PARAM.inp.esolver_type != "tddft" && inp.out_freq_ion > 0) // default value of out_freq_ion is 0
    {
        if (istep % inp.out_freq_ion == 0)
        {
            istep_in = istep;
            out_flag = true;
        }
    }
    else if (PARAM.inp.esolver_type == "tddft" && inp.out_freq_td > 0) // default value of out_freq_td is 0
    {
        if (istep % inp.out_freq_td == 0)
        {
            istep_in = istep;
            out_flag = true;
        }
    }
    else if (conv_esolver || scf_nmax_flag) // mohan add scf_nmax_flag on 20250921
    {
        out_flag = true;
    }

    if (!out_flag)
    {
        ModuleBase::timer::end("ModuleIO", "ctrl_scf_lcao");
        return;
    }

    //*****

    const bool out_app_flag = inp.out_app_flag;
    const bool gamma_only = PARAM.globalv.gamma_only_local;
    const int nspin = inp.nspin;
    const std::string global_out_dir = PARAM.globalv.global_out_dir;

    //------------------------------------------------------------------
    //! 1) print out density of states (DOS)
    //------------------------------------------------------------------
    if (inp.out_dos)
    {
        ModuleIO::write_dos_lcao(psi,
                                 p_hamilt,
                                 pv,
                                 ucell,
                                 kv,
                                 inp.nbands,
                                 pelec->eferm,
                                 pelec->ekb,
                                 pelec->wg,
                                 inp.dos_edelta_ev,
                                 inp.dos_scale,
                                 inp.dos_sigma,
                                 out_app_flag,
                                 istep,
                                 GlobalV::ofs_running);
    }

    //------------------------------------------------------------------
    //! 2) Output density matrix DM(R)
    //------------------------------------------------------------------
    if (inp.out_dmr[0])
    {
        const int precision = inp.out_dmr[1];

        ModuleIO::write_dmr(dm->get_DMR_vector(), &ucell, precision, pv, out_app_flag, 
			ucell.get_iat2iwt(), ucell.nat, istep);
    }

    //------------------------------------------------------------------
    //! 3) Output density matrix DM(k)
    //------------------------------------------------------------------
    if (inp.out_dmk[0])
    {
        std::vector<double> efermis(nspin == 2 ? 2 : 1);
        for (int ispin = 0; ispin < efermis.size(); ispin++)
        {
            efermis[ispin] = pelec->eferm.get_efval(ispin);
        }
        const int precision = inp.out_dmk[1];

        ModuleIO::write_dmk(dm->get_DMK_vector(), kv, precision, efermis, &(ucell), pv, istep);
    }

    //------------------------------------------------------------------
    // 4) Output H(k) and S(k) matrices for each k-point
    //------------------------------------------------------------------
    if (inp.out_mat_hs[0])
    {
        ModuleIO::write_hsk(global_out_dir,
                            nspin,
                            kv.get_nks(),
                            kv.get_nkstot(),
                            kv.ik2iktot,
                            kv.isk,
                            p_hamilt,
                            pv,
                            gamma_only,
                            out_app_flag,
                            istep,
                            GlobalV::ofs_running);
    }

    //------------------------------------------------------------------
    //! 5) Output electronic wavefunctions Psi(k)
    //------------------------------------------------------------------
    if (elecstate::ElecStateLCAO<TK>::out_wfc_lcao)
    {
        ModuleIO::write_wfc_nao(elecstate::ElecStateLCAO<TK>::out_wfc_lcao,
                                out_app_flag,
                                psi[0],
                                pelec->ekb,
                                pelec->wg,
                                kv.kvec_c,
                                kv.ik2iktot,
                                kv.get_nkstot(),
                                pv,
                                nspin,
                                istep);
    }

    //------------------------------------------------------------------
    //! 6) Output DeePKS information
    //------------------------------------------------------------------
#ifdef __MLALGO
    // need control parameter
    hamilt::HamiltLCAO<TK, TR>* p_ham_deepks = p_hamilt;
    LCAO_Deepks_Interface<TK, TR> deepks_interface(&deepks.ld);

    deepks_interface.out_deepks_labels(pelec->f_en.etot,
                                       kv.get_nks(),
                                       ucell.nat,
                                       PARAM.globalv.nlocal,
                                       pelec->ekb,
                                       kv.kvec_d,
                                       ucell,
                                       orb,
                                       gd,
                                       &pv,
                                       *psi,
                                       dm,
                                       p_ham_deepks,
                                       -1,   // -1 when called in after scf
                                       true, // no used when after scf
                                       GlobalV::MY_RANK,
                                       GlobalV::ofs_running);
#endif

    //------------------------------------------------------------------
    //! 7a) Output H(R) and S(R) matrices in CSR format
    //------------------------------------------------------------------
    if (inp.out_mat_hs2[0])
    {
        const int precision = inp.out_mat_hs2[1];
        std::vector<hamilt::HContainer<TR>*> hr_vec = p_hamilt->getHR_vector();
        const hamilt::HContainer<TR>* sr = p_hamilt->getSR();

        ModuleIO::write_hsr(hr_vec, sr, &ucell, precision, pv,
                            out_app_flag, ucell.get_iat2iwt(), ucell.nat, istep);
    }

    //------------------------------------------------------------------
    //! 7a.1) Output H(R), S(R), and DM(R) matrices in NPZ format
    //------------------------------------------------------------------
    if (inp.out_hsr_npz)
    {
        std::string zipname = PARAM.globalv.global_out_dir + "output_SR.npz";
        ModuleIO::output_mat_npz(ucell, zipname, *(p_hamilt->getSR()));
    }

    if (inp.out_hr_npz || inp.out_hsr_npz)
    {
        std::vector<hamilt::HContainer<TR>*> hr_vec = p_hamilt->getHR_vector();
        for (int ispin = 0; ispin < hr_vec.size(); ++ispin)
        {
            std::string zipname
                = PARAM.globalv.global_out_dir + "output_HR" + std::to_string(ispin) + ".npz";
            ModuleIO::output_mat_npz(ucell, zipname, *(hr_vec[ispin]));
        }
    }

    if (inp.out_dm_npz)
    {
        const std::vector<hamilt::HContainer<double>*>& dmr_vec = dm->get_DMR_vector();
        for (int ispin = 0; ispin < dmr_vec.size(); ++ispin)
        {
            std::string zipname
                = PARAM.globalv.global_out_dir + "output_DM" + std::to_string(ispin) + ".npz";
            ModuleIO::output_mat_npz(ucell, zipname, *(dmr_vec[ispin]));
        }
    }

    //------------------------------------------------------------------
    //! 7b) Output dH, dS, T, r matrices (old sparse path, without H/S), only for multi-k
    //------------------------------------------------------------------
    hamilt::Hamilt<TK>* p_ham_tk = static_cast<hamilt::Hamilt<TK>*>(p_hamilt);

    ModuleIO::MatSparseOutputOptions mat_sparse_options;
    mat_sparse_options.out_mat_dh = inp.out_mat_dh[0];
    mat_sparse_options.out_mat_ds = inp.out_mat_ds[0];
    mat_sparse_options.out_mat_t = inp.out_mat_t[0];
    mat_sparse_options.out_mat_r = inp.out_mat_r[0];
    mat_sparse_options.dh_precision = inp.out_mat_dh[1];
    mat_sparse_options.ds_precision = inp.out_mat_ds[1];
    mat_sparse_options.t_precision = inp.out_mat_t[1];
    mat_sparse_options.r_precision = inp.out_mat_r[1];

    if(!PARAM.globalv.gamma_only_local)
    ModuleIO::output_mat_sparse(mat_sparse_options,
                                istep,
                                pelec->pot->get_eff_v(),
                                pv,
                                two_center_bundle,
                                orb,
                                ucell,
                                gd,
                                kv,
                                p_ham_tk,
                                &dftu);

    //------------------------------------------------------------------
    //! 7c) Output atomic dH components (dT/dτ, dV^NL/dτ, dV^L/dτ, dV^H/dτ, dV^XC/dτ), only for nspin =1, 2 now
    //------------------------------------------------------------------
    if( PARAM.inp.nspin < 4 )
    {
        WriteDHParams dh_params;
        dh_params.ucell = &ucell;
        dh_params.gd = &gd;
        dh_params.pv = &pv;
        dh_params.two_center_bundle = &two_center_bundle;
        dh_params.orb = &orb;
        dh_params.kv = &kv;
        dh_params.v_eff = &pelec->pot->get_eff_v();
        dh_params.pot = pelec->pot;
        dh_params.chg = pelec->charge;
        // pelec->pot->get_eff_v() is the SUM V^L + V^H + V^XC; feeding it to cal_dH would
        // give the wrong potential for the separated V^L / V^H / V^XC outputs. Build one
        // dedicated Potential per term with exactly one component registered (see write_vxc.hpp).
        double dh_etxc = 0.0;
        double dh_vtxc = 0.0;
        elecstate::Potential* pot_vl = nullptr;
        elecstate::Potential* pot_vh = nullptr;
        elecstate::Potential* pot_vxc = nullptr;
        // out_mat_dh (total dH = sum of all terms) needs every veff potential regardless of the
        // per-component flags, so allocate all three when it is on; otherwise allocate per flag.
        if (inp.out_mat_dh_vl[0] || inp.out_mat_dh[0])
        {
            pot_vl = new elecstate::Potential(pw_rhod, pw_rho, &ucell, &vloc,
                const_cast<Structure_Factor*>(&sf), &solvent, &dh_etxc, &dh_vtxc);
            pot_vl->pot_register({"local"});
            pot_vl->update_from_charge(pelec->charge, &ucell);
        }
        if (inp.out_mat_dh_vh[0] || inp.out_mat_dh[0])
        {
            pot_vh = new elecstate::Potential(pw_rhod, pw_rho, &ucell, &vloc,
                const_cast<Structure_Factor*>(&sf), &solvent, &dh_etxc, &dh_vtxc);
            pot_vh->pot_register({"hartree"});
            pot_vh->update_from_charge(pelec->charge, &ucell);
        }
        if (inp.out_mat_dh_vxc[0] || inp.out_mat_dh[0])
        {
            pot_vxc = new elecstate::Potential(pw_rhod, pw_rho, &ucell, &vloc,
                const_cast<Structure_Factor*>(&sf), &solvent, &dh_etxc, &dh_vtxc);
            pot_vxc->pot_register({"xc"});
            pot_vxc->update_from_charge(pelec->charge, &ucell);
        }
        dh_params.pot_vl = pot_vl;
        dh_params.pot_vh = pot_vh;
        dh_params.pot_vxc = pot_vxc;
        dh_params.iat2iwt = ucell.get_iat2iwt();
        dh_params.nat = ucell.nat;
        dh_params.nspin = inp.nspin;
        dh_params.istep = istep;
        dh_params.gamma_only = gamma_only;
        dh_params.append = out_app_flag;
        if (PARAM.inp.nspin == 1 || PARAM.inp.nspin == 2)
        {
            // per-spin DM (1-indexed): nspin=1 -> {spin0}, nspin=2 -> {spin-up, spin-down}.
            // The Veff Hellmann-Feynman terms need these (V^H sums spins, V^XC is spin-resolved).
            for (int is = 1; is <= PARAM.inp.nspin; ++is)
            {
                dh_params.dmR.push_back(dm->get_DMR_pointer(is));
            }
        }
#ifdef __EXX
        // dV^EXX/dR output is wired for the gamma (TK==double) exx interfaces. exd/exc are
        // mutually exclusive (real vs complex Hexx); write_dH_exx picks by info_ri.real_number.
        setup_exx_dh_params(dh_params, exx_nao);
#endif
        ModuleIO::write_dH_components(dh_params);
        delete pot_vl;
        delete pot_vh;
        delete pot_vxc;
    }


    //------------------------------------------------------------------
    //! 7d) Output H components (T, Vnl, Vl, Vh, Vxc)
    //------------------------------------------------------------------
    {
        ModuleIO::WriteHParams h_params;
        h_params.ucell = &ucell;
        h_params.gd = &gd;
        h_params.pv = &pv;
        h_params.two_center_bundle = &two_center_bundle;
        h_params.orb = &orb;
        h_params.kv = &kv;
        h_params.pot = pelec->pot;
        h_params.chg = pelec->charge;
        h_params.rho_basis = pw_rho;
        h_params.nrxx = pw_rho->nrxx;
        h_params.nspin = nspin;
        h_params.istep = istep;
        h_params.append = out_app_flag;
        h_params.iat2iwt = ucell.get_iat2iwt();
        h_params.nat = ucell.nat;
        if (inp.out_mat_h_t[0])
        {
            ModuleIO::write_h_t(h_params);
        }
        if (inp.out_mat_h_vnl[0])
        {
            ModuleIO::write_h_vnl(h_params);
        }
        if (inp.out_mat_h_vl[0])
        {
            ModuleIO::write_h_vl(h_params);
        }
        if (inp.out_mat_h_vh[0])
        {
            ModuleIO::write_h_vh(h_params);
        }
        if (inp.out_mat_h_vxc[0])
        {
            ModuleIO::write_h_vxc(h_params);
        }
#ifdef __EXX
        if (inp.out_mat_h_exx[0] && GlobalC::exx_info.info_global.cal_exx)
        {
            // V^EXX(R) output is wired for the gamma (TK==double) exx interfaces.
            setup_exx_h_params(h_params, exx_nao);
        }
#endif
    }
    //------------------------------------------------------------------
    //! 8) Output kinetic matrix
    //------------------------------------------------------------------
    if (inp.out_mat_tk[0])
    {
        hamilt::HS_Matrix_K<TK> hsk(&pv, true);
        hamilt::HContainer<TR> hR(&pv);
        hamilt::Operator<TK>* ekinetic
            = new hamilt::EKinetic<hamilt::OperatorLCAO<TK, TR>>(&hsk,
                                                                    kv.kvec_d,
                                                                    &hR,
                                                                    &ucell,
                                                                    orb.cutoffs(),
                                                                    &gd,
                                                                    two_center_bundle.kinetic_orb.get());

        const int nspin_k = (nspin == 2 ? 2 : 1);
        for (int ik = 0; ik < kv.get_nks() / nspin_k; ++ik)
        {
            ekinetic->init(ik);

            const int out_label = 1; // 1: .txt, 2: .dat

            std::string t_fn = ModuleIO::filename_output(global_out_dir,
                                                         "tk",
                                                         "nao",
                                                         ik,
                                                         kv.ik2iktot,
                                                         inp.nspin,
                                                         kv.get_nkstot(),
                                                         out_label,
                                                         out_app_flag,
                                                         gamma_only,
                                                         istep);

            ModuleIO::save_mat(istep,
                               hsk.get_hk(),
                               PARAM.globalv.nlocal,
                               false, // bit
                               inp.out_mat_tk[1],
                               1, // true for upper triangle matrix
                               inp.out_app_flag,
                               t_fn,
                               pv,
                               GlobalV::DRANK);
        }

        delete ekinetic;
    }

    //------------------------------------------------------------------
    //! 9) Output expectation of angular momentum operator
    //------------------------------------------------------------------
    if (inp.out_mat_l[0])
    {
        ModuleIO::AngularMomentumCalculator mylcalculator(inp.orbital_dir,
                                                          ucell,
                                                          orb.get_rcutmax_Phi(),
                                                          inp.test_deconstructor,
                                                          inp.test_grid,
                                                          inp.test_atom_input,
                                                          PARAM.globalv.search_pbc,
                                                          &GlobalV::ofs_running,
                                                          GlobalV::MY_RANK);
        mylcalculator.calculate(inp.suffix, global_out_dir, ucell, inp.out_mat_l[1], GlobalV::MY_RANK);
    }

    //------------------------------------------------------------------
    //! 10) Output Mulliken charge
    //------------------------------------------------------------------
    if (inp.out_mul)
    {
        ModuleIO::cal_mag(&pv,
                          p_hamilt,
                          kv,
                          dm, // mohan add 2025-11-04
                          two_center_bundle,
                          orb,
                          ucell,
                          gd,
                          istep,
                          true);
    }

    //------------------------------------------------------------------
    //! 11) Output atomic magnetization by using 'spin_constraint'
    //------------------------------------------------------------------
    if (inp.sc_mag_switch)
    {
        spinconstrain::SpinConstrain<TK>& sc = spinconstrain::SpinConstrain<TK>::getScInstance();
        sc.cal_mi_lcao(istep);
        sc.print_Mi(GlobalV::ofs_running);
        sc.print_Mag_Force(GlobalV::ofs_running);
    }

    //------------------------------------------------------------------
    //! 12) Output Berry phase
    //------------------------------------------------------------------
    if (inp.calculation == "nscf" && berryphase::berry_phase_flag && ModuleSymmetry::Symmetry::symm_flag != 1)
    {
        std::cout << FmtCore::format("\n * * * * * *\n << Start %s.\n", "Berry phase calculation");
        berryphase bp(&pv);
        bp.lcao_init(ucell, gd, kv, orb);
        // additional step before calling macroscopic_polarization
        bp.Macroscopic_polarization(ucell, pw_wfc->npwk_max, psi, pw_rho, pw_wfc, kv);
        std::cout << FmtCore::format(" >> Finish %s.\n * * * * * *\n", "Berry phase calculation");
    }

    //------------------------------------------------------------------
    //! 13) Wannier90 interface in LCAO basis
    // added by jingan in 2018.11.7
    //------------------------------------------------------------------
    if (inp.calculation == "nscf" && inp.towannier90)
    {
        std::cout << FmtCore::format("\n * * * * * *\n << Start %s.\n", "Wave function to Wannier90");
        if (inp.wannier_method == 1)
        {
            toWannier90_LCAO_IN_PW wan(inp.out_wannier_mmn,
                                       inp.out_wannier_amn,
                                       inp.out_wannier_unk,
                                       inp.out_wannier_eig,
                                       inp.out_wannier_wvfn_formatted,
                                       inp.nnkpfile,
                                       inp.wannier_spin);
            wan.set_tpiba_omega(ucell.tpiba, ucell.omega);
            wan.calculate(ucell, pelec->ekb, pw_wfc, pw_big, sf, kv, psi, &pv);
        }
        else if (inp.wannier_method == 2)
        {
            toWannier90_LCAO wan(inp.out_wannier_mmn,
                                 inp.out_wannier_amn,
                                 inp.out_wannier_unk,
                                 inp.out_wannier_eig,
                                 inp.out_wannier_wvfn_formatted,
                                 inp.nnkpfile,
                                 inp.wannier_spin,
                                 orb);

            wan.calculate(ucell, gd, pelec->ekb, kv, *psi, &pv);
        }
        std::cout << FmtCore::format(" >> Finish %s.\n * * * * * *\n", "Wave function to Wannier90");
    }

    // 14) calculate the kinetic energy density tau
    // mohan add 2025-10-24
    //    if (inp.out_elf[0] > 0)
    //	{
    //		LCAO_domain::dm2tau(pelec->DM->get_DMR_vector(), inp.nspin, pelec->charge);
    //	}

#ifdef __EXX
    //------------------------------------------------------------------
    //! 15) Output Hexx matrix in LCAO basis
    // (see `out_chg` in docs/advanced/input_files/input-main.md)
    //------------------------------------------------------------------
    bool cal_exx = GlobalC::exx_info.info_global.cal_exx;
    bool real_number = GlobalC::exx_info.info_ri.real_number;

    if (inp.out_chg[0])
    {
        if (cal_exx && inp.calculation != "nscf") // Peize Lin add if 2022.11.14
        {
            const std::string file_name_exx = global_out_dir + "HexxR" + std::to_string(GlobalV::MY_RANK);
            if (real_number)
            {
                ModuleIO::write_Hexxs_csr(file_name_exx, ucell, exx_nao.exd->get_Hexxs());
            }
            else
            {
                ModuleIO::write_Hexxs_csr(file_name_exx, ucell, exx_nao.exc->get_Hexxs());
            }
        }
    }

    //------------------------------------------------------------------
    //! 16) Write RPA information in LCAO basis
    //------------------------------------------------------------------
    if (inp.rpa)
    {
        RPA_LRI<TK, double> rpa_lri_double(GlobalC::exx_info.info_ri);
        rpa_lri_double.postSCF(ucell, MPI_COMM_WORLD, *dm, pelec, kv, orb, pv, *psi);
    }
#endif

    //------------------------------------------------------------------
    //! 17) Perform RDMFT calculations, added by jghan, 2024-10-17
    //------------------------------------------------------------------
    if (inp.rdmft == true)
    {
        ModuleBase::matrix occ_num(pelec->wg);
        for (int ik = 0; ik < occ_num.nr; ++ik)
        {
            for (int inb = 0; inb < occ_num.nc; ++inb)
            {
                occ_num(ik, inb) /= kv.wk[ik];
            }
        }
        rdmft_solver.update_elec(ucell, occ_num, *psi);

        //! initialize the gradients of Etotal with respect to occupation numbers and wfc,
        //! and set all elements to 0.
        //! dedocc = d E/d Occ_Num
        ModuleBase::matrix dedocc(pelec->wg.nr, pelec->wg.nc, true);

        //! dedwfc = d E/d wfc
        psi::Psi<TK> dedwfc(psi->get_nk(), psi->get_nbands(), psi->get_nbasis(), kv.ngk, true);
        dedwfc.zero_out();

        double etot_rdmft = rdmft_solver.run(dedocc, dedwfc);
    }

    //------------------------------------------------------------------
    //! 17) Output quasi orbitals
    //------------------------------------------------------------------
    if (inp.qo_switch)
    {
        toQO tqo(inp.qo_basis, inp.qo_strategy, inp.qo_thr, inp.qo_screening_coeff);
        tqo.initialize(global_out_dir,
                       inp.pseudo_dir,
                       inp.orbital_dir,
                       &ucell,
                       kv.kvec_d,
                       GlobalV::ofs_running,
                       GlobalV::MY_RANK,
                       GlobalV::NPROC);
        tqo.calculate();
    }

    //------------------------------------------------------------------
    //! 18) Calculate and output asynchronous overlap matrix for Hefei-NAMD
    //------------------------------------------------------------------
    if (inp.cal_syns[0] > 0 && (istep > 0 || inp.init_vel))
    {
        ModuleBase::TITLE("ModuleIO", "output_namd_async_overlap");
        ModuleBase::timer::start("ModuleIO", "output_namd_async_overlap");

        // Create a new Overlap instance specifically for SR_async calculation
        // This allows SR_async to be initialized with velocity-shifted dtau
        hamilt::Overlap<hamilt::OperatorLCAO<TK, TR>>* overlap_async =
            new hamilt::Overlap<hamilt::OperatorLCAO<TK, TR>>(
                nullptr,  // hsk_in: not needed for SR_async calculation
                kv.kvec_d,
                nullptr,  // hR_in: not needed for SR_async calculation
                nullptr,  // SR_in: not needed for SR_async calculation
                &ucell,
                orb.cutoffs(),
                &gd,
                two_center_bundle.overlap_orb.get());

        // Use precision from cal_syns[1] (default 8 if not specified)
        const int precision = inp.cal_syns[1];
        const Parallel_Orbitals* paraV = p_hamilt->getSR()->get_paraV();
        hamilt::HContainer<TR>* SR_async = overlap_async->calculate_SR_async(ucell, PARAM.mdp.md_dt, paraV);
        overlap_async->output_SR_async_csr(istep, SR_async, precision);

        // Clean up
        delete SR_async;
        delete overlap_async;

        ModuleBase::timer::end("ModuleIO", "output_namd_async_overlap");
    }

    ModuleBase::timer::end("ModuleIO", "ctrl_scf_lcao");
}

// For gamma only
template void ModuleIO::ctrl_scf_lcao<double, double>(
    UnitCell& ucell,
    const Input_para& inp,
    K_Vectors& kv,
    elecstate::ElecState* pelec,
    elecstate::DensityMatrix<double, double>* dm, // mohan add 2025-11-04
    Parallel_Orbitals& pv,
    Grid_Driver& gd,
    psi::Psi<double>* psi,
    hamilt::HamiltLCAO<double, double>* p_hamilt,
    Plus_U& dftu, // mohan add 2025-11-07
    TwoCenterBundle& two_center_bundle,
    LCAO_Orbitals& orb,
    const ModulePW::PW_Basis_K* pw_wfc,         // for berryphase
    const ModulePW::PW_Basis* pw_rho,           // for berryphase
    const ModulePW::PW_Basis_Big* pw_big,       // for Wannier90
    const Structure_Factor& sf,                 // for Wannier90
    const ModulePW::PW_Basis* pw_rhod,          // dense charge grid (for dH veff pots)
    const ModuleBase::matrix& vloc,             // local pseudopotential (for dH veff pots)
    surchem& solvent,                           // solvent model (for dH veff pots)
    rdmft::RDMFT<double, double>& rdmft_solver, // for RDMFT
    Setup_DeePKS<double>& deepks,
    Exx_NAO<double>& exx_nao,
    const bool conv_esolver,
    const bool scf_nmax_flag,
    const int istep);

// For multiple k-points
template void ModuleIO::ctrl_scf_lcao<std::complex<double>, double>(
    UnitCell& ucell,
    const Input_para& inp,
    K_Vectors& kv,
    elecstate::ElecState* pelec,
    elecstate::DensityMatrix<std::complex<double>, double>* dm, // mohan add 2025-11-04
    Parallel_Orbitals& pv,
    Grid_Driver& gd,
    psi::Psi<std::complex<double>>* psi,
    hamilt::HamiltLCAO<std::complex<double>, double>* p_hamilt,
    Plus_U& dftu, // mohan add 2025-11-07
    TwoCenterBundle& two_center_bundle,
    LCAO_Orbitals& orb,
    const ModulePW::PW_Basis_K* pw_wfc,                       // for berryphase
    const ModulePW::PW_Basis* pw_rho,                         // for berryphase
    const ModulePW::PW_Basis_Big* pw_big,                     // for Wannier90
    const Structure_Factor& sf,                               // for Wannier90
    const ModulePW::PW_Basis* pw_rhod,                        // dense charge grid (for dH veff pots)
    const ModuleBase::matrix& vloc,                           // local pseudopotential (for dH veff pots)
    surchem& solvent,                                         // solvent model (for dH veff pots)
    rdmft::RDMFT<std::complex<double>, double>& rdmft_solver, // for RDMFT
    Setup_DeePKS<std::complex<double>>& deepks,
    Exx_NAO<std::complex<double>>& exx_nao,
    const bool conv_esolver,
    const bool scf_nmax_flag,
    const int istep);

template void ModuleIO::ctrl_scf_lcao<std::complex<double>, std::complex<double>>(
    UnitCell& ucell,
    const Input_para& inp,
    K_Vectors& kv,
    elecstate::ElecState* pelec,
    elecstate::DensityMatrix<std::complex<double>, double>* dm, // mohan add 2025-11-04
    Parallel_Orbitals& pv,
    Grid_Driver& gd,
    psi::Psi<std::complex<double>>* psi,
    hamilt::HamiltLCAO<std::complex<double>, std::complex<double>>* p_hamilt,
    Plus_U& dftu, // mohan add 2025-11-07
    TwoCenterBundle& two_center_bundle,
    LCAO_Orbitals& orb,
    const ModulePW::PW_Basis_K* pw_wfc,                                     // for berryphase
    const ModulePW::PW_Basis* pw_rho,                                       // for berryphase
    const ModulePW::PW_Basis_Big* pw_big,                                   // for Wannier90
    const Structure_Factor& sf,                                             // for Wannier90
    const ModulePW::PW_Basis* pw_rhod,                                      // dense charge grid (for dH veff pots)
    const ModuleBase::matrix& vloc,                                         // local pseudopotential (for dH veff pots)
    surchem& solvent,                                                       // solvent model (for dH veff pots)
    rdmft::RDMFT<std::complex<double>, std::complex<double>>& rdmft_solver, // for RDMFT
    Setup_DeePKS<std::complex<double>>& deepks,
    Exx_NAO<std::complex<double>>& exx_nao,
    const bool conv_esolver,
    const bool scf_nmax_flag,
    const int istep);
