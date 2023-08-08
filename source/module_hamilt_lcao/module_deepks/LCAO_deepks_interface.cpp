#ifdef __DEEPKS
#include "LCAO_deepks_interface.h"

#include "module_base/global_variable.h"
#include "module_base/tool_title.h"
#include "module_elecstate/cal_dm.h"

LCAO_Deepks_Interface::LCAO_Deepks_Interface(std::shared_ptr<LCAO_Deepks> ld_in) : ld(ld_in)
{
}

void LCAO_Deepks_Interface::out_deepks_labels(double etot,
                                              int nks,
                                              int nat,
                                              const ModuleBase::matrix& ekb,
                                              const std::vector<ModuleBase::Vector3<double>>& kvec_d,
                                              const UnitCell& ucell,
                                              const LCAO_Orbitals& orb,
                                              Grid_Driver& GridD,
                                              const Parallel_Orbitals* ParaV,
                                              const psi::Psi<std::complex<double>>& psi,
                                              const psi::Psi<double>& psid,
                                              const std::vector<ModuleBase::matrix>& dm_gamma,
                                              const std::vector<ModuleBase::ComplexMatrix>& dm_k)
{
    ModuleBase::TITLE("LCAO_Deepks_Interface", "out_deepks_labels");
    // calculating deepks correction to bandgap
    // and save the results
    if (GlobalV::deepks_out_labels) // caoyu add 2021-06-04
    {
        ld->save_npy_e(etot, "e_tot.npy");
        if (GlobalV::deepks_scf)
        {
            ld->save_npy_e(etot - ld->E_delta,
                           "e_base.npy"); // ebase :no deepks E_delta including
        }
        else // deepks_scf = 0; base calculation
        {
            ld->save_npy_e(etot, "e_base.npy"); // no scf, e_tot=e_base
        }

        if (GlobalV::deepks_bandgap)
        {
            int nocc = GlobalV::nelec / 2;
            ModuleBase::matrix deepks_bands;
            deepks_bands.create(nks, 1);
            for (int iks = 0; iks < nks; iks++)
            {
                for (int hl = 0; hl < 1; hl++)
                {
                    deepks_bands(iks, hl) = ekb(iks, nocc + hl) - ekb(iks, nocc - 1 + hl);
                }
            }
            ld->save_npy_o(deepks_bands, "o_tot.npy", nks);
            if (GlobalV::deepks_scf)
            {
                int nocc = GlobalV::nelec / 2;
                ModuleBase::matrix wg_hl;
                if (GlobalV::GAMMA_ONLY_LOCAL)
                {
                    wg_hl.create(GlobalV::NSPIN, GlobalV::NBANDS);
                    std::vector<std::vector<ModuleBase::matrix>> dm_bandgap_gamma;
                    dm_bandgap_gamma.resize(GlobalV::NSPIN);
                    for (int is = 0; is < GlobalV::NSPIN; is++)
                    {
                        for (int ib = 0; ib < 1; ib++)
                        {
                            wg_hl.zero_out();
                            wg_hl(is, ib + nocc - 1) = -1.0;
                            wg_hl(is, ib + nocc) = 1.0;
                            dm_bandgap_gamma[ib].resize(GlobalV::NSPIN);
                            elecstate::cal_dm(ParaV, wg_hl, psid, dm_bandgap_gamma[ib]);
                        }
                    }

                    ld->cal_orbital_precalc(dm_bandgap_gamma, nat, ucell, orb, GridD);

                    ld->save_npy_orbital_precalc(nat, nks);
                    ld->cal_o_delta(dm_bandgap_gamma);
                    ld->save_npy_o(deepks_bands - ld->o_delta, "o_base.npy", nks);
                }    // end deepks_scf gamma-only;
                else // multi-k bandgap label
                {
                    wg_hl.create(nks, GlobalV::NBANDS);
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
                    ld->save_npy_orbital_precalc(nat, nks);
                    ld->cal_o_delta_k(dm_bandgap_k, nks);
                    ld->save_npy_o(deepks_bands - ld->o_delta, "o_base.npy", nks);
                } // end deepks_scf multi-k
            }     // end deepks_scf == 1
            else  // deepks_scf == 0
            {
                ld->save_npy_o(deepks_bands, "o_base.npy", nks); // no scf, o_tot=o_base
            }                                                    // end deepks_scf == 0
        }                                                        // end bandgap label
        if(GlobalV::deepks_v_delta)//gamma only now
        {
            if(GlobalV::GAMMA_ONLY_LOCAL)
            {
                ModuleBase::matrix h_tot;
                h_tot.create(GlobalV::NLOCAL,GlobalV::NLOCAL);

                ld->collect_h_mat(ld->h_mat,h_tot,GlobalV::NLOCAL);
                ld->save_npy_h(h_tot, "h_tot.npy",GlobalV::NLOCAL);

                if(GlobalV::deepks_scf)
                {
                    ModuleBase::matrix v_delta;
                    v_delta.create(GlobalV::NLOCAL,GlobalV::NLOCAL);
                    ld->collect_h_mat(ld->H_V_delta,v_delta,GlobalV::NLOCAL);
                    ld->save_npy_h(h_tot-v_delta, "h_base.npy",GlobalV::NLOCAL);
                    ld->save_npy_h(v_delta, "v_delta.npy",GlobalV::NLOCAL);

                    if(GlobalV::deepks_v_delta==1)//v_delta_precalc storage method 1
                    {
                        ld->cal_v_delta_precalc(GlobalV::NLOCAL,
                                nat,
                                ucell,
                                orb,
                                GridD);
                    
                        ld->save_npy_v_delta_precalc(nat, 1, GlobalV::NLOCAL);
                    }
                    else if(GlobalV::deepks_v_delta==2)//v_delta_precalc storage method 2:overlap_part
                    {
                        ld->prepare_psialpha(GlobalV::NLOCAL,
                                    nat,
                                    ucell,
                                    orb,
                                    GridD);
                        
                        ld->save_npy_psialpha(nat, 1, GlobalV::NLOCAL);

                        ld->prepare_gevdm(
                                    nat,
                                    orb);

                        ld->save_npy_gevdm(nat);
                    }
                }
                else //deepks_scf == 0
                {
                    ld->save_npy_h(h_tot, "h_base.npy",GlobalV::NLOCAL);
                }
            }
            else
            {
                ModuleBase::WARNING_QUIT("ESolver_KS_LCAO", "V_delta label has not been developed for multi-k now!");
            }
        }//end v_delta label
    }                                                            // end deepks_out_labels

    // DeePKS PDM and descriptor
    if (GlobalV::deepks_out_labels || GlobalV::deepks_scf)
    {
        // this part is for integrated test of deepks
        // so it is printed no matter even if deepks_out_labels is not used
        if (GlobalV::GAMMA_ONLY_LOCAL)
        {
            ld->cal_projected_DM(dm_gamma, ucell, orb, GridD);
        }
        else
        {
            ld->cal_projected_DM_k(dm_k, ucell, orb, GridD, nks, kvec_d);
        }
        ld->check_projected_dm(); // print out the projected dm for NSCF calculaiton
        ld->cal_descriptor();     // final descriptor
        ld->check_descriptor(ucell);

        if (GlobalV::deepks_out_labels)
            ld->save_npy_d(nat); // libnpy needed
    }
    //
    if (GlobalV::deepks_scf)
    {
        if (GlobalV::GAMMA_ONLY_LOCAL)
        {
            ld->cal_e_delta_band(dm_gamma);
        }
        else
        {
            ld->cal_e_delta_band_k(dm_k, nks);
        }
        std::cout << "E_delta_band = " << std::setprecision(8) << ld->e_delta_band << " Ry"
                  << " = " << std::setprecision(8) << ld->e_delta_band * ModuleBase::Ry_to_eV << " eV"
                  << std::endl;
        std::cout << "E_delta_NN= " << std::setprecision(8) << ld->E_delta << " Ry"
                  << " = " << std::setprecision(8) << ld->E_delta * ModuleBase::Ry_to_eV << " eV" << std::endl;
    }
}

#endif
