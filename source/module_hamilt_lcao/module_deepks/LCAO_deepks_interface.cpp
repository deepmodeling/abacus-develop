#ifdef __DEEPKS
#include "LCAO_deepks_interface.h"
#include "module_parameter/parameter.h"
#include "LCAO_deepks_io.h" // mohan add 2024-07-22

#include "module_base/global_variable.h"
#include "module_base/tool_title.h"
#include "module_elecstate/cal_dm.h"

LCAO_Deepks_Interface::LCAO_Deepks_Interface(std::shared_ptr<LCAO_Deepks> ld_in) : ld(ld_in)
{
}
// gamma-only
void LCAO_Deepks_Interface::out_deepks_labels(const double& etot,
                                              const int& nks,
                                              const int& nat,
                                              const int& nlocal,
                                              const ModuleBase::matrix& ekb,
                                              const std::vector<ModuleBase::Vector3<double>>& kvec_d,
                                              const UnitCell& ucell,
                                              const LCAO_Orbitals& orb,
                                              Grid_Driver& GridD,
                                              const Parallel_Orbitals* ParaV,
                                              const psi::Psi<double>& psid,
                                              const elecstate::DensityMatrix<double, double>* dm,
                                              const int& deepks_v_delta)
{
    ModuleBase::TITLE("LCAO_Deepks_Interface", "out_deepks_labels");

    const int my_rank = GlobalV::MY_RANK;
    const int nspin = PARAM.inp.nspin;

    // calculating deepks correction to bandgap
    // and save the results
    if (PARAM.inp.deepks_out_labels)
    {
        // mohan updated 2024-07-25
        const std::string file_etot = PARAM.globalv.global_out_dir + "deepks_etot.npy";
        const std::string file_ebase = PARAM.globalv.global_out_dir + "deepks_ebase.npy";

        LCAO_deepks_io::save_npy_e(etot, file_etot, my_rank);

        if (PARAM.inp.deepks_scf)
        {
            /// ebase :no deepks E_delta including
            LCAO_deepks_io::save_npy_e(etot - ld->E_delta,
                           file_ebase, my_rank);
        }
        else // deepks_scf = 0; base calculation
        {
            /// no scf, e_tot=e_base
            LCAO_deepks_io::save_npy_e(etot, file_ebase, my_rank);
        }

        if (PARAM.inp.deepks_bandgap)
        {
            const int nocc = PARAM.inp.nelec / 2;
            ModuleBase::matrix deepks_bands;
            deepks_bands.create(nks, 1);
            for (int iks = 0; iks < nks; ++iks)
            {
                // record band gap for each k point (including spin)
                for (int hl = 0; hl < 1; ++hl)
                {
                    deepks_bands(iks, hl) = ekb(iks, nocc + hl) - ekb(iks, nocc - 1 + hl);
                }
            }

            const std::string file_otot = PARAM.globalv.global_out_dir + "deepks_otot.npy";
            LCAO_deepks_io::save_npy_o(deepks_bands, file_otot, nks, my_rank);

            if (PARAM.inp.deepks_scf)
            {
                ModuleBase::matrix wg_hl;
                wg_hl.create(nspin, PARAM.inp.nbands);

                std::vector<std::vector<ModuleBase::matrix>> dm_bandgap_gamma;

                dm_bandgap_gamma.resize(nspin);
                for (int is = 0; is < nspin; ++is)
                {
                    for (int ib = 0; ib < 1; ++ib)
                    {
                        wg_hl.zero_out();
                        wg_hl(is, ib + nocc - 1) = -1.0;
                        wg_hl(is, ib + nocc) = 1.0;
                        dm_bandgap_gamma[ib].resize(nspin);
                        elecstate::cal_dm(ParaV, wg_hl, psid, dm_bandgap_gamma[ib]);
                    }
                }

                ld->cal_orbital_precalc(dm_bandgap_gamma, nat, ucell, orb, GridD);

                LCAO_deepks_io::save_npy_orbital_precalc(nat, 
                                                         nks, 
                                                         ld->des_per_atom, 
                                                         ld->orbital_precalc_tensor, 
                                                         PARAM.globalv.global_out_dir,
                                                         my_rank);

                ld->cal_o_delta(dm_bandgap_gamma);

                const std::string file_obase = PARAM.globalv.global_out_dir + "deepks_obase.npy";
                LCAO_deepks_io::save_npy_o(deepks_bands - ld->o_delta, file_obase, nks, my_rank);
            }     // end deepks_scf == 1
            else  // deepks_scf == 0
            {
                const std::string file_obase = PARAM.globalv.global_out_dir + "deepks_obase.npy";
                LCAO_deepks_io::save_npy_o(deepks_bands, file_obase, nks, my_rank); // no scf, o_tot=o_base
            }                                                    // end deepks_scf == 0
        } // end bandgap label

        if(deepks_v_delta)//gamma only now
        {
            ModuleBase::matrix h_tot;
            h_tot.create(nlocal,nlocal);

            DeePKS_domain::collect_h_mat(*ParaV, ld->h_mat,h_tot,nlocal);

            const std::string file_htot = PARAM.globalv.global_out_dir + "deepks_htot.npy";
            LCAO_deepks_io::save_npy_h(h_tot, file_htot, nlocal, my_rank);

            if(PARAM.inp.deepks_scf)
            {
                ModuleBase::matrix v_delta;
                v_delta.create(nlocal,nlocal);
                DeePKS_domain::collect_h_mat(*ParaV, ld->H_V_delta,v_delta,nlocal);

                const std::string file_hbase = PARAM.globalv.global_out_dir + "deepks_hbase.npy";
                LCAO_deepks_io::save_npy_h(h_tot-v_delta, file_hbase, nlocal, my_rank);

                const std::string file_vdelta = PARAM.globalv.global_out_dir + "deepks_vdelta.npy";
                LCAO_deepks_io::save_npy_h(v_delta, file_vdelta, nlocal, my_rank);

                if(deepks_v_delta==1)//v_delta_precalc storage method 1
                {
                    ld->cal_v_delta_precalc(nlocal,
                            nat,
                            ucell,
                            orb,
                            GridD);
                
                    const int nks_gamma = 1;
                    LCAO_deepks_io::save_npy_v_delta_precalc(
                      nat, 
                      nks_gamma, 
                      nlocal,
                      ld->des_per_atom,
                      ld->v_delta_precalc_tensor,
                      PARAM.globalv.global_out_dir,
                      my_rank);
                }
                else if(deepks_v_delta==2)//v_delta_precalc storage method 2
                {
                    ld->prepare_psialpha(nlocal,
                                nat,
                                ucell,
                                orb,
                                GridD);

                    const int nks_gamma=1;
                    LCAO_deepks_io::save_npy_psialpha(nat, 
                                nks_gamma, 
                                nlocal,
                                ld->inlmax,
                                ld->lmaxd,
                                ld->psialpha_tensor,
                                PARAM.globalv.global_out_dir,
                                my_rank);

                    ld->prepare_gevdm(
                                nat,
                                orb);

                    LCAO_deepks_io::save_npy_gevdm(nat,
                      ld->inlmax,
                      ld->lmaxd,
                      ld->gevdm_tensor,
                      PARAM.globalv.global_out_dir,
                      my_rank);
                }
            }
            else //deepks_scf == 0
            {
                const std::string file_hbase = PARAM.globalv.global_out_dir + "deepks_hbase.npy";
                LCAO_deepks_io::save_npy_h(h_tot, file_hbase, nlocal, my_rank);
            }
        }//end v_delta label
    
    } // end deepks_out_labels


    // DeePKS PDM and descriptor
    if (PARAM.inp.deepks_out_labels || PARAM.inp.deepks_scf)
    {
        // this part is for integrated test of deepks
        // when deepks_scf is on, the init pdm should be same as the out pdm, so we should not recalculate the pdm
		if(!PARAM.inp.deepks_scf) 
		{
			ld->cal_projected_DM(dm, ucell, orb, GridD);
		}

        ld->check_projected_dm(); // print out the projected dm for NSCF calculaiton

        ld->cal_descriptor(nat);     // final descriptor

        ld->check_descriptor(ucell, PARAM.globalv.global_out_dir);

		if (PARAM.inp.deepks_out_labels)
		{
			LCAO_deepks_io::save_npy_d(
					nat, 
					ld->des_per_atom, 
					ld->inlmax, 
                    ld->inl_l,
					PARAM.inp.deepks_equiv, 
					ld->d_tensor, 
                    PARAM.globalv.global_out_dir,
					my_rank); // libnpy needed
		}
    }
    
    /// print out deepks information to the screen
    if (PARAM.inp.deepks_scf)
    {
        ld->cal_e_delta_band(dm->get_DMK_vector());
        std::cout << "E_delta_band = " << std::setprecision(8) << ld->e_delta_band << " Ry"
                  << " = " << std::setprecision(8) << ld->e_delta_band * ModuleBase::Ry_to_eV << " eV"
                  << std::endl;
        std::cout << "E_delta_NN= " << std::setprecision(8) << ld->E_delta << " Ry"
                  << " = " << std::setprecision(8) << ld->E_delta * ModuleBase::Ry_to_eV << " eV" << std::endl;
    }
}

// multi-k
void LCAO_Deepks_Interface::out_deepks_labels(const double& etot,
                                              const int& nks,
                                              const int& nat,
                                              const int& nlocal,
                                              const ModuleBase::matrix& ekb,
                                              const std::vector<ModuleBase::Vector3<double>>& kvec_d,
                                              const UnitCell& ucell,
                                              const LCAO_Orbitals& orb,
                                              Grid_Driver& GridD,
                                              const Parallel_Orbitals* ParaV,
                                              const psi::Psi<std::complex<double>>& psi,
                                              const elecstate::DensityMatrix<std::complex<double>, double>* dm,
                                              const int& deepks_v_delta)
{
    ModuleBase::TITLE("LCAO_Deepks_Interface", "out_deepks_labels");
    ModuleBase::timer::tick("LCAO_Deepks_Interface", "out_deepks_labels");

    const int my_rank = GlobalV::MY_RANK;
    const int nspin = PARAM.inp.nspin;

    /// calculating deepks correction to bandgap and save the results
    if (PARAM.inp.deepks_out_labels)
    {
        // mohan updated 2024-07-25
        const std::string file_etot = PARAM.globalv.global_out_dir + "deepks_etot.npy";
        const std::string file_ebase = PARAM.globalv.global_out_dir + "deepks_ebase.npy";

        LCAO_deepks_io::save_npy_e(etot, file_etot, my_rank);

        if (PARAM.inp.deepks_scf)
        {
            /// ebase :no deepks E_delta including
            LCAO_deepks_io::save_npy_e(etot - ld->E_delta,
                           file_ebase, my_rank);
        }
        else // deepks_scf = 0; base calculation
        {
            LCAO_deepks_io::save_npy_e(etot, file_ebase, my_rank);
        }

        if (PARAM.inp.deepks_bandgap)
        {
            int nocc = PARAM.inp.nelec / 2;
            ModuleBase::matrix deepks_bands;
            deepks_bands.create(nks, 1);
            for (int iks = 0; iks < nks; iks++)
            {
                for (int hl = 0; hl < 1; hl++)
                {
                    deepks_bands(iks, hl) = ekb(iks, nocc + hl) - ekb(iks, nocc - 1 + hl);
                }
            }

            const std::string file_otot = PARAM.globalv.global_out_dir + "deepks_otot.npy";
            LCAO_deepks_io::save_npy_o(deepks_bands, file_otot, nks, my_rank);

            if (PARAM.inp.deepks_scf)
            {
                int nocc = PARAM.inp.nelec / 2; // redundant!
                ModuleBase::matrix wg_hl;
                wg_hl.create(nks, PARAM.inp.nbands);
                std::vector<std::vector<ModuleBase::ComplexMatrix>> dm_bandgap_k;
                dm_bandgap_k.resize(1);

                for (int ib = 0; ib < 1; ib++)
                {
                    wg_hl.zero_out();
                    for (int ik = 0; ik < nks; ik++)
                    {
                        wg_hl(ik, ib + nocc - 1) = -1.0;
                        wg_hl(ik, ib + nocc) = 1.0;
                    }
                    dm_bandgap_k[ib].resize(nks);
                    elecstate::cal_dm(ParaV, wg_hl, psi, dm_bandgap_k[ib]);
                }

                // ld->cal_o_delta_k(dm_bandgap_k, ParaV, nks);
                ld->cal_orbital_precalc_k(dm_bandgap_k, nat, nks, kvec_d, ucell, orb, GridD);

				LCAO_deepks_io::save_npy_orbital_precalc(
						nat, 
						nks, 
						ld->des_per_atom, 
						ld->orbital_precalc_tensor, 
                        PARAM.globalv.global_out_dir,
						GlobalV::MY_RANK);

                ld->cal_o_delta_k(dm_bandgap_k, nks);

                const std::string file_obase = PARAM.globalv.global_out_dir + "deepks_obase.npy";
                LCAO_deepks_io::save_npy_o(deepks_bands - ld->o_delta, file_obase, nks, my_rank);
            }     // end deepks_scf == 1
            else  // deepks_scf == 0
            {
                const std::string file_obase = PARAM.globalv.global_out_dir + "deepks_obase.npy";
                LCAO_deepks_io::save_npy_o(deepks_bands, file_obase, nks, my_rank); // no scf, o_tot=o_base
            }     // end deepks_scf == 0
        } // end bandgap label
        if(deepks_v_delta)
        {
            std::vector<ModuleBase::ComplexMatrix> h_tot(nks);
            for (int ik = 0; ik < nks; ik++)
            {
                h_tot[ik].create(nlocal, nlocal);
            }

            DeePKS_domain::collect_h_mat(*ParaV, ld->h_mat_k,h_tot,nlocal,nks);

            const std::string file_htot = PARAM.globalv.global_out_dir + "deepks_htot.npy";
            LCAO_deepks_io::save_npy_h(h_tot, file_htot, nlocal, nks, my_rank);

            if(PARAM.inp.deepks_scf)
            {
                std::vector<ModuleBase::ComplexMatrix> v_delta(nks);
                std::vector<ModuleBase::ComplexMatrix> hbase(nks);
                for (int ik = 0; ik < nks; ik++)
                {
                    v_delta[ik].create(nlocal, nlocal);
                    hbase[ik].create(nlocal, nlocal);
                }
                DeePKS_domain::collect_h_mat(*ParaV, ld->H_V_delta_k,v_delta,nlocal,nks); 

                const std::string file_hbase = PARAM.globalv.global_out_dir + "deepks_hbase.npy";
                for (int ik = 0; ik < nks; ik++)
                {
                    hbase[ik] = h_tot[ik] - v_delta[ik];
                }
                LCAO_deepks_io::save_npy_h(hbase, file_hbase, nlocal, nks, my_rank);

                const std::string file_vdelta = PARAM.globalv.global_out_dir + "deepks_vdelta.npy";
                LCAO_deepks_io::save_npy_h(v_delta, file_vdelta, nlocal, nks, my_rank);

                if(deepks_v_delta==1)//v_delta_precalc storage method 1
                {
                    ld->cal_v_delta_precalc_k(nlocal,
                            nat,
                            nks,
                            kvec_d,
                            ucell,
                            orb,
                            GridD);
                
                    LCAO_deepks_io::save_npy_v_delta_precalc(
                      nat, 
                      nks, 
                      nlocal,
                      ld->des_per_atom,
                      ld->v_delta_precalc_tensor,
                      PARAM.globalv.global_out_dir,
                      my_rank);

                }
                else if(deepks_v_delta==2)//v_delta_precalc storage method 2
                {
                    ld->prepare_psialpha_k(nlocal,
                                nat,
                                nks,
                                kvec_d,
                                ucell,
                                orb,
                                GridD);
                    
                    LCAO_deepks_io::save_npy_psialpha(nat, 
                                nks, 
                                nlocal,
                                ld->inlmax,
                                ld->lmaxd,
                                ld->psialpha_tensor,
                                PARAM.globalv.global_out_dir,
                                my_rank);

                    ld->prepare_gevdm(
                                nat,
                                orb);

                    LCAO_deepks_io::save_npy_gevdm(nat,
                      ld->inlmax,
                      ld->lmaxd,
                      ld->gevdm_tensor,
                      PARAM.globalv.global_out_dir,
                      my_rank);
                }
            }
            else //deepks_scf == 0
            {
                const std::string file_hbase = PARAM.globalv.global_out_dir + "deepks_hbase.npy";
                LCAO_deepks_io::save_npy_h(h_tot, file_hbase, nlocal, nks, my_rank);
            }
        }
    }  // end deepks_out_labels


    // DeePKS PDM and descriptor
    if (PARAM.inp.deepks_out_labels || PARAM.inp.deepks_scf)
    {
        // this part is for integrated test of deepks
        // so it is printed no matter even if deepks_out_labels is not used
        // when deepks_scf is on, the init pdm should be same as the out pdm, so we should not recalculate the pdm
		if(!PARAM.inp.deepks_scf) 
		{
			ld->cal_projected_DM_k(dm, ucell, orb, GridD);
		}

        ld->check_projected_dm(); // print out the projected dm for NSCF calculaiton

        ld->cal_descriptor(nat);     // final descriptor

        ld->check_descriptor(ucell, PARAM.globalv.global_out_dir);

        if (PARAM.inp.deepks_out_labels)
        {
            LCAO_deepks_io::save_npy_d(nat, 
                                       ld->des_per_atom, 
									   ld->inlmax, 
									   ld->inl_l,
									   PARAM.inp.deepks_equiv, 
                                       ld->d_tensor, 
                                       PARAM.globalv.global_out_dir,
                                       GlobalV::MY_RANK); // libnpy needed
        }
    }
    //
    if (PARAM.inp.deepks_scf)
    {
        ld->cal_e_delta_band_k(dm->get_DMK_vector(), nks);

        std::cout << "E_delta_band = " << std::setprecision(8) << ld->e_delta_band << " Ry"
                  << " = " << std::setprecision(8) << ld->e_delta_band * ModuleBase::Ry_to_eV << " eV"
                  << std::endl;

        std::cout << "E_delta_NN= " << std::setprecision(8) << ld->E_delta << " Ry"
                  << " = " << std::setprecision(8) << ld->E_delta * ModuleBase::Ry_to_eV << " eV" << std::endl;
    }

    ModuleBase::timer::tick("LCAO_Deepks_Interface", "out_deepks_labels");
}

#endif
