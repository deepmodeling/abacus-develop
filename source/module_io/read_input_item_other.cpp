#include "read_input.h"
#include "read_input_tool.h"

namespace ModuleIO
{
void ReadInput::item_others()
{
    // 10. Electric field and dipole correction
    {
        Input_Item item("efield_flag");
        item.annotation = "add electric field";
        read_sync_bool(efield_flag);
        this->add_item(item);
    }
    {
        Input_Item item("dip_cor_flag");
        item.annotation = "dipole correction";
        read_sync_bool(dip_cor_flag);
        this->add_item(item);
    }
    {
        Input_Item item("efield_dir");
        item.annotation = "the direction of the electric field or dipole correction";
        read_sync_int(efield_dir);
        this->add_item(item);
    }
    {
        Input_Item item("efield_pos_max");
        item.annotation = "position of the maximum of the saw-like potential along crystal axis efield_dir";
        read_sync_double(efield_pos_max);
        this->add_item(item);
    }
    {
        Input_Item item("efield_pos_dec");
        item.annotation = "zone in the unit cell where the saw-like potential decreases";
        read_sync_double(efield_pos_dec);
        this->add_item(item);
    }
    {
        Input_Item item("efield_amp");
        item.annotation = "amplitude of the electric field";
        read_sync_double(efield_amp);
        this->add_item(item);
    }

    // 11. Gate field
    {
        Input_Item item("gate_flag");
        item.annotation = "compensating charge or not";
        read_sync_bool(gate_flag);
        this->add_item(item);
    }
    {
        Input_Item item("zgate");
        item.annotation = "position of charged plate";
        read_sync_double(zgate);
        this->add_item(item);
    }
    {
        Input_Item item("relax");
        item.annotation = "allow relaxation along the specific direction";
        read_sync_bool(relax);
        this->add_item(item);
    }
    {
        Input_Item item("block");
        item.annotation = "add a block potential or not";
        read_sync_bool(block);
        this->add_item(item);
    }
    {
        Input_Item item("block_down");
        item.annotation = "low bound of the block";
        read_sync_double(block_down);
        this->add_item(item);
    }
    {
        Input_Item item("block_up");
        item.annotation = "high bound of the block";
        read_sync_double(block_up);
        this->add_item(item);
    }
    {
        Input_Item item("block_height");
        item.annotation = "height of the block";
        read_sync_double(block_height);
        this->add_item(item);
    }

    // 12. Test
    {
        Input_Item item("out_alllog");
        item.annotation = "output information for each processor, when parallel";
        read_sync_bool(out_alllog);
        this->add_item(item);
    }
    {
        Input_Item item("nurse");
        item.annotation = "for coders";
        read_sync_int(nurse);
        this->add_item(item);
    }
    {
        Input_Item item("colour");
        item.annotation = "for coders, make their live colourful";
        read_sync_bool(colour);
        this->add_item(item);
    }
    {
        Input_Item item("t_in_h");
        item.annotation = "calculate the kinetic energy or not";
        read_sync_bool(t_in_h);
        this->add_item(item);
    }
    {
        Input_Item item("vl_in_h");
        item.annotation = "calculate the local potential or not";
        read_sync_bool(vl_in_h);
        this->add_item(item);
    }
    {
        Input_Item item("vnl_in_h");
        item.annotation = "calculate the nonlocal potential or not";
        read_sync_bool(vnl_in_h);
        this->add_item(item);
    }
    {
        Input_Item item("vh_in_h");
        item.annotation = "calculate the hartree potential or not";
        read_sync_bool(vh_in_h);
        this->add_item(item);
    }
    {
        Input_Item item("vion_in_h");
        item.annotation = "calculate the local ionic potential or not";
        read_sync_bool(vion_in_h);
        this->add_item(item);
    }
    {
        Input_Item item("test_force");
        item.annotation = "test the force";
        read_sync_bool(test_force);
        this->add_item(item);
    }
    {
        Input_Item item("test_stress");
        item.annotation = "test the stress";
        read_sync_bool(test_stress);
        this->add_item(item);
    }

    // 13. vdW Correction
    {
        Input_Item item("vdw_method");
        item.annotation = "the method of calculating vdw (none ; d2 ; d3_0 ; d3_bj";
        read_sync_string(vdw_method);
        this->add_item(item);
    }
    {
        Input_Item item("vdw_s6");
        item.annotation = "scale parameter of d2/d3_0/d3_bj";
        read_sync_string(vdw_s6);
        this->add_item(item);
    }
    {
        Input_Item item("vdw_s8");
        item.annotation = "scale parameter of d3_0/d3_bj";
        read_sync_string(vdw_s8);
        this->add_item(item);
    }
    {
        Input_Item item("vdw_a1");
        item.annotation = "damping parameter of d3_0/d3_bj";
        read_sync_string(vdw_a1);
        this->add_item(item);
    }
    {
        Input_Item item("vdw_a2");
        item.annotation = "damping parameter of d3_bj";
        read_sync_string(vdw_a2);
        this->add_item(item);
    }
    {
        Input_Item item("vdw_d");
        item.annotation = "damping parameter of d2";
        read_sync_double(vdw_d);
        this->add_item(item);
    }
    {
        Input_Item item("vdw_abc");
        item.annotation = "third-order term?";
        read_sync_bool(vdw_abc);
        this->add_item(item);
    }
    {
        Input_Item item("vdw_C6_file");
        item.annotation = "filename of C6";
        read_sync_string(vdw_C6_file);
        this->add_item(item);
    }
    {
        Input_Item item("vdw_C6_unit");
        item.annotation = "unit of C6, Jnm6/mol or eVA6";
        read_sync_string(vdw_C6_unit);
        this->add_item(item);
    }
    {
        Input_Item item("vdw_R0_file");
        item.annotation = "filename of R0";
        read_sync_string(vdw_R0_file);
        this->add_item(item);
    }
    {
        Input_Item item("vdw_R0_unit");
        item.annotation = "unit of R0, A or Bohr";
        read_sync_string(vdw_R0_unit);
        this->add_item(item);
    }
    {
        Input_Item item("vdw_cutoff_type");
        item.annotation = "expression model of periodic structure, radius or period";
        read_sync_string(vdw_cutoff_type);
        this->add_item(item);
    }
    {
        Input_Item item("vdw_cutoff_radius");
        item.annotation = "radius cutoff for periodic structure";
        read_sync_string(vdw_cutoff_radius);
        this->add_item(item);
    }
    {
        Input_Item item("vdw_radius_unit");
        item.annotation = "unit of radius cutoff for periodic structure";
        read_sync_string(vdw_radius_unit);
        this->add_item(item);
    }
    {
        Input_Item item("vdw_cn_thr");
        item.annotation = "radius cutoff for cn";
        read_sync_double(vdw_cn_thr);
        this->add_item(item);
    }
    {
        Input_Item item("vdw_cn_thr_unit");
        item.annotation = "unit of cn_thr, Bohr or Angstrom";
        read_sync_string(vdw_cn_thr_unit);
        this->add_item(item);
    }
    {
        Input_Item item("vdw_cutoff_period");
        item.annotation = "periods of periodic structure";
        item.readvalue = [](const Input_Item& item, Parameter& para) {
            int count = item.str_values.size();
            if (count == 3)
            {
                para.input.vdw_cutoff_period[0] = convertstr<int>(item.str_values[0]);
                para.input.vdw_cutoff_period[1] = convertstr<int>(item.str_values[1]);
                para.input.vdw_cutoff_period[2] = convertstr<int>(item.str_values[2]);
            }
            else
            {
                throw std::runtime_error("vdw_cutoff_period should have 3 values");
            }
        };
        item.getfinalvalue = [](Input_Item& item, const Parameter& para) {
            item.final_value << para.input.vdw_cutoff_period[0] << " " << para.input.vdw_cutoff_period[1] << " "
                             << para.input.vdw_cutoff_period[2];
        };
#ifdef __MPI
        bcastfuncs.push_back(
            [](Parameter& para) { Parallel_Common::bcast_int((int*)&para.input.vdw_cutoff_period, 3); });
#endif
        this->add_item(item);
    }

    // 14. exx
    {
        Input_Item item("exx_hybrid_alpha");
        item.annotation = "fraction of Fock exchange in hybrid functionals";
        read_sync_string(exx_hybrid_alpha);
        this->add_item(item);
    }
    {
        Input_Item item("exx_hse_omega");
        item.annotation = "range-separation parameter in HSE functional";
        read_sync_double(exx_hse_omega);
        this->add_item(item);
    }
    {
        Input_Item item("exx_separate_loop");
        item.annotation
            = "if 1, a two-step method is employed, else it will start with a GGA-Loop, and then Hybrid-Loop";
        read_sync_bool(exx_separate_loop);
        this->add_item(item);
    }
    {
        Input_Item item("exx_hybrid_step");
        item.annotation = "the maximal electronic iteration number in the evaluation of Fock exchange";
        read_sync_int(exx_hybrid_step);
        this->add_item(item);
    }
    {
        Input_Item item("exx_mixing_beta");
        item.annotation = "mixing_beta for outer-loop when exx_separate_loop=1";
        read_sync_double(exx_mixing_beta);
        this->add_item(item);
    }
    {
        Input_Item item("exx_lambda");
        item.annotation = "used to compensate for divergence points at G=0 in the evaluation of Fock exchange using "
                          "lcao_in_pw method";
        read_sync_double(exx_lambda);
        this->add_item(item);
    }
    {
        Input_Item item("exx_real_number");
        item.annotation = "exx calculated in real or complex";
        read_sync_string(exx_real_number);
        this->add_item(item);
    }
    {
        Input_Item item("exx_pca_threshold");
        item.annotation = "threshold to screen on-site ABFs in exx";
        read_sync_double(exx_pca_threshold);
        this->add_item(item);
    }
    {
        Input_Item item("exx_c_threshold");
        item.annotation = "threshold to screen C matrix in exx";
        read_sync_double(exx_c_threshold);
        this->add_item(item);
    }
    {
        Input_Item item("exx_v_threshold");
        item.annotation = "threshold to screen C matrix in exx";
        read_sync_double(exx_v_threshold);
        this->add_item(item);
    }
    {
        Input_Item item("exx_dm_threshold");
        item.annotation = "threshold to screen density matrix in exx";
        read_sync_double(exx_dm_threshold);
        this->add_item(item);
    }
    {
        Input_Item item("exx_cauchy_threshold");
        item.annotation = "threshold to screen exx using Cauchy-Schwartz inequality";
        read_sync_double(exx_cauchy_threshold);
        this->add_item(item);
    }
    {
        Input_Item item("exx_c_grad_threshold");
        item.annotation = "threshold to screen nabla C matrix in exx";
        read_sync_double(exx_c_grad_threshold);
        this->add_item(item);
    }
    {
        Input_Item item("exx_v_grad_threshold");
        item.annotation = "threshold to screen nabla V matrix in exx";
        read_sync_double(exx_v_grad_threshold);
        this->add_item(item);
    }
    {
        Input_Item item("exx_cauchy_force_threshold");
        item.annotation = "threshold to screen exx force using Cauchy-Schwartz inequality";
        read_sync_double(exx_cauchy_force_threshold);
        this->add_item(item);
    }
    {
        Input_Item item("exx_cauchy_stress_threshold");
        item.annotation = "threshold to screen exx stress using Cauchy-Schwartz inequality";
        read_sync_double(exx_cauchy_stress_threshold);
        this->add_item(item);
    }
    {
        Input_Item item("exx_ccp_rmesh_times");
        item.annotation = "how many times larger the radial mesh required for calculating Columb potential is to that "
                          "of atomic orbitals";
        read_sync_string(exx_ccp_rmesh_times);
        this->add_item(item);
    }
    {
        Input_Item item("exx_distribute_type");
        item.annotation = "exx_distribute_type";
        read_sync_string(exx_distribute_type);
        this->add_item(item);
    }
    {
        Input_Item item("exx_opt_orb_lmax");
        item.annotation = "the maximum l of the spherical Bessel functions for opt ABFs";
        read_sync_int(exx_opt_orb_lmax);
        this->add_item(item);
    }
    {
        Input_Item item("exx_opt_orb_ecut");
        item.annotation = "the cut-off of plane wave expansion for opt ABFs";
        read_sync_double(exx_opt_orb_ecut);
        this->add_item(item);
    }
    {
        Input_Item item("exx_opt_orb_tolerence");
        item.annotation = "the threshold when solving for the zeros of spherical Bessel functions for opt ABFs";
        read_sync_double(exx_opt_orb_tolerence);
        this->add_item(item);
    }
    {
        Input_Item item("rpa_ccp_rmesh_times");
        item.annotation = "how many times larger the radial mesh required for calculating Columb potential is to that "
                          "of atomic orbitals";
        read_sync_double(rpa_ccp_rmesh_times);
        this->add_item(item);
    }

    // 16. tddft
    {
        Input_Item item("td_force_dt");
        item.annotation = "time of force change";
        read_sync_double(td_force_dt);
        this->add_item(item);
    }
    {
        Input_Item item("td_vext");
        item.annotation = "add extern potential or not";
        read_sync_bool(td_vext);
        this->add_item(item);
    }
    {
        Input_Item item("td_vext_dire");
        item.annotation = "extern potential direction";
        item.readvalue = [](const Input_Item& item, Parameter& para) {
            int count = item.str_values.size();
            para.input.td_nvext_dire = count;
            for (auto& str: item.str_values)
            {
                para.input.td_vext_dire.push_back(convertstr<int>(str));
            }
        };
        add_int_bcast(td_nvext_dire); // Since "td_nvext_dire" has been assigned a value, it needs to be broadcasted
        // We must firt bcast td_nvext_direr, then bcast td_vext_dire
        sync_intvec(td_vext_dire, para.input.td_nvext_dire);
        this->add_item(item);
    }
    {
        Input_Item item("out_dipole");
        item.annotation = "output dipole or not";
        read_sync_bool(out_dipole);
        this->add_item(item);
    }
    {
        Input_Item item("out_efield");
        item.annotation = "output dipole or not";
        read_sync_bool(out_efield);
        this->add_item(item);
    }
    {
        Input_Item item("out_current");
        item.annotation = "output current or not";
        read_sync_bool(out_current);
        this->add_item(item);
    }
    {
        Input_Item item("out_vecpot");
        item.annotation = "output TDDFT vector potential or not";
        read_sync_bool(out_vecpot);
        this->add_item(item);
    }
    {
        Input_Item item("init_vecpot_file");
        item.annotation = "init vector potential through file or not";
        read_sync_bool(init_vecpot_file);
        this->add_item(item);
    }
    {
        Input_Item item("td_print_eij");
        item.annotation = "print eij or not";
        read_sync_double(td_print_eij);
        this->add_item(item);
    }
    {
        Input_Item item("td_edm");
        item.annotation = "the method to calculate the energy density matrix";
        read_sync_int(td_edm);
        this->add_item(item);
    }
    {
        Input_Item item("td_propagator");
        item.annotation = "method of propagator";
        read_sync_int(propagator);
        this->add_item(item);
    }
    {
        Input_Item item("td_stype");
        item.annotation = "type of electric field in space domain";
        read_sync_int(td_stype);
        this->add_item(item);
    }
    {
        Input_Item item("td_ttype");
        item.annotation = "type of electric field in time domain";
        read_sync_string(td_ttype);
        this->add_item(item);
    }
    {
        Input_Item item("td_tstart");
        item.annotation = " number of steps where electric field starts";
        read_sync_int(td_tstart);
        this->add_item(item);
    }
    {
        Input_Item item("td_tend");
        item.annotation = "number of steps where electric field ends";
        read_sync_int(td_tend);
        this->add_item(item);
    }
    {
        Input_Item item("td_lcut1");
        item.annotation = "cut1 of interval in length gauge";
        read_sync_double(td_lcut1);
        this->add_item(item);
    }
    {
        Input_Item item("td_lcut2");
        item.annotation = "cut2 of interval in length gauge";
        read_sync_double(td_lcut2);
        this->add_item(item);
    }
    {
        Input_Item item("td_gauss_freq");
        item.annotation = "frequency (freq) of Gauss type electric field";
        read_sync_double(td_gauss_freq);
        this->add_item(item);
    }
    {
        Input_Item item("td_gauss_phase");
        item.annotation = "phase of Gauss type electric field";
        read_sync_double(td_gauss_phase);
        this->add_item(item);
    }
    {
        Input_Item item("td_gauss_sigma");
        item.annotation = "sigma of Gauss type electric field";
        read_sync_double(td_gauss_sigma);
        this->add_item(item);
    }
    {
        Input_Item item("td_gauss_t0");
        item.annotation = "step number of time center (t0) of Gauss type electric field";
        read_sync_double(td_gauss_t0);
        this->add_item(item);
    }
    {
        Input_Item item("td_gauss_amp");
        item.annotation = "amplitude of Gauss type electric field";
        read_sync_double(td_gauss_amp);
        this->add_item(item);
    }
    {
        Input_Item item("td_trape_freq");
        item.annotation = "frequency of Trapezoid type electric field";
        read_sync_double(td_trape_freq);
        this->add_item(item);
    }
    {
        Input_Item item("td_trape_phase");
        item.annotation = "phase of Trapezoid type electric field";
        read_sync_double(td_trape_phase);
        this->add_item(item);
    }
    {
        Input_Item item("td_trape_t1");
        item.annotation = "t1 of Trapezoid type electric field";
        read_sync_double(td_trape_t1);
        this->add_item(item);
    }
    {
        Input_Item item("td_trape_t2");
        item.annotation = "t2 of Trapezoid type electric field";
        read_sync_double(td_trape_t2);
        this->add_item(item);
    }
    {
        Input_Item item("td_trape_t3");
        item.annotation = "t3 of Trapezoid type electric field";
        read_sync_double(td_trape_t3);
        this->add_item(item);
    }
    {
        Input_Item item("td_trape_amp");
        item.annotation = "amplitude of Trapezoid type electric field";
        read_sync_double(td_trape_amp);
        this->add_item(item);
    }
    {
        Input_Item item("td_trigo_freq1");
        item.annotation = "frequency 1 of Trigonometric type electric field";
        read_sync_double(td_trigo_freq1);
        this->add_item(item);
    }
    {
        Input_Item item("td_trigo_freq2");
        item.annotation = "frequency 2 of Trigonometric type electric field";
        read_sync_double(td_trigo_freq2);
        this->add_item(item);
    }
    {
        Input_Item item("td_trigo_phase1");
        item.annotation = "phase 1 of Trigonometric type electric field";
        read_sync_double(td_trigo_phase1);
        this->add_item(item);
    }
    {
        Input_Item item("td_trigo_phase2");
        item.annotation = "phase 2 of Trigonometric type electric field";
        read_sync_double(td_trigo_phase2);
        this->add_item(item);
    }
    {
        Input_Item item("td_trigo_amp");
        item.annotation = "amplitude of Trigonometric type electric field";
        read_sync_double(td_trigo_amp);
        this->add_item(item);
    }
    {
        Input_Item item("td_heavi_t0");
        item.annotation = "t0 of Heaviside type electric field";
        read_sync_double(td_heavi_t0);
        this->add_item(item);
    }
    {
        Input_Item item("td_heavi_amp");
        item.annotation = "amplitude of Heaviside type electric field";
        read_sync_double(td_heavi_amp);
        this->add_item(item);
    }
    {
        Input_Item item("ocp");
        item.annotation = "change occupation or not";
        read_sync_bool(ocp);
        this->add_item(item);
    }
    {
        Input_Item item("ocp_set");
        item.annotation = "set occupation";
        read_sync_string(ocp_set);
        this->add_item(item);
    }

    // 17. berry_wannier
    {
        Input_Item item("berry_phase");
        item.annotation = "calculate berry phase or not";
        read_sync_bool(berry_phase);
        this->add_item(item);
    }
    {
        Input_Item item("gdir");
        item.annotation = "calculate the polarization in the direction of the lattice vector";
        read_sync_int(gdir);
        this->add_item(item);
    }
    {
        Input_Item item("towannier90");
        item.annotation = "use wannier90 code interface or not";
        read_sync_bool(towannier90);
        this->add_item(item);
    }
    {
        Input_Item item("nnkpfile");
        item.annotation = "the wannier90 code nnkp file name";
        read_sync_string(nnkpfile);
        this->add_item(item);
    }
    {
        Input_Item item("wannier_spin");
        item.annotation = "calculate spin in wannier90 code interface";
        read_sync_string(wannier_spin);
        this->add_item(item);
    }
    {
        Input_Item item("wannier_method");
        item.annotation = "different implementation methods under Lcao basis set";
        read_sync_int(wannier_method);
        this->add_item(item);
    }
    {
        Input_Item item("out_wannier_mmn");
        item.annotation = "output .mmn file or not";
        read_sync_bool(out_wannier_mmn);
        this->add_item(item);
    }
    {
        Input_Item item("out_wannier_amn");
        item.annotation = "output .amn file or not";
        read_sync_bool(out_wannier_amn);
        this->add_item(item);
    }
    {
        Input_Item item("out_wannier_unk");
        item.annotation = "output UNK. file or not";
        read_sync_bool(out_wannier_unk);
        this->add_item(item);
    }
    {
        Input_Item item("out_wannier_eig");
        item.annotation = "output .eig file or not";
        read_sync_bool(out_wannier_eig);
        this->add_item(item);
    }
    {
        Input_Item item("out_wannier_wvfn_formatted");
        item.annotation = "output UNK. file in text format or in binary format";
        read_sync_bool(out_wannier_wvfn_formatted);
        this->add_item(item);
    }

    // 18. imlicit_solvation
    {
        Input_Item item("imp_sol");
        item.annotation = "calculate implicit solvation correction or not";
        read_sync_bool(imp_sol);
        this->add_item(item);
    }
    {
        Input_Item item("eb_k");
        item.annotation = "the relative permittivity of the bulk solvent";
        read_sync_double(eb_k);
        this->add_item(item);
    }
    {
        Input_Item item("tau");
        item.annotation = "the effective surface tension parameter";
        read_sync_double(tau);
        this->add_item(item);
    }
    {
        Input_Item item("sigma_k");
        item.annotation = "the width of the diffuse cavity";
        read_sync_double(sigma_k);
        this->add_item(item);
    }
    {
        Input_Item item("nc_k");
        item.annotation = "the cut-off charge density";
        read_sync_double(nc_k);
        this->add_item(item);
    }

    // 19. OFDFT
    {
        Input_Item item("of_kinetic");
        item.annotation = "kinetic energy functional, such as tf, vw, wt";
        read_sync_string(of_kinetic);
        this->add_item(item);
    }
    {
        Input_Item item("of_method");
        item.annotation = "optimization method used in OFDFT, including cg1, cg2, tn (default)";
        read_sync_string(of_method);
        this->add_item(item);
    }
    {
        Input_Item item("of_conv");
        item.annotation = "the convergence criterion, potential, energy (default), or both";
        read_sync_string(of_conv);
        this->add_item(item);
    }
    {
        Input_Item item("of_tole");
        item.annotation = "tolerance of the energy change (in Ry) for determining the convergence, default=2e-6 Ry";
        read_sync_double(of_tole);
        this->add_item(item);
    }
    {
        Input_Item item("of_tolp");
        item.annotation = "tolerance of potential for determining the convergence, default=1e-5 in a.u.";
        read_sync_double(of_tolp);
        this->add_item(item);
    }
    {
        Input_Item item("of_tf_weight");
        item.annotation = "weight of TF KEDF";
        read_sync_double(of_tf_weight);
        this->add_item(item);
    }
    {
        Input_Item item("of_vw_weight");
        item.annotation = "weight of vW KEDF";
        read_sync_double(of_vw_weight);
        this->add_item(item);
    }
    {
        Input_Item item("of_wt_alpha");
        item.annotation = "parameter alpha of WT KEDF";
        read_sync_double(of_wt_alpha);
        this->add_item(item);
    }
    {
        Input_Item item("of_wt_beta");
        item.annotation = "parameter beta of WT KEDF";
        read_sync_double(of_wt_beta);
        this->add_item(item);
    }
    {
        Input_Item item("of_wt_rho0");
        item.annotation = "the average density of system, used in WT KEDF, in Bohr^-3";
        read_sync_double(of_wt_rho0);
        this->add_item(item);
    }
    {
        Input_Item item("of_hold_rho0");
        item.annotation = "If set to 1, the rho0 will be fixed even if the volume of system has changed, it will be "
                          "set to 1 automaticly if of_wt_rho0 is not zero";
        read_sync_bool(of_hold_rho0);
        this->add_item(item);
    }
    {
        Input_Item item("of_lkt_a");
        item.annotation = "parameter a of LKT KEDF";
        read_sync_double(of_lkt_a);
        this->add_item(item);
    }
    {
        Input_Item item("of_full_pw");
        item.annotation
            = "If set to 1, ecut will be ignored when collect planewaves, so that all planewaves will be used";
        read_sync_bool(of_full_pw);
        this->add_item(item);
    }
    {
        Input_Item item("of_full_pw_dim");
        item.annotation = "If of_full_pw = true, dimention of FFT is testricted to be (0) either odd or even; (1) odd "
                          "only; (2) even only";
        read_sync_int(of_full_pw_dim);
        this->add_item(item);
    }
    {
        Input_Item item("of_read_kernel");
        item.annotation = "If set to 1, the kernel of WT KEDF will be filled from file of_kernel_file, not from "
                          "formula. Only usable for WT KEDF";
        read_sync_bool(of_read_kernel);
        this->add_item(item);
    }
    {
        Input_Item item("of_kernel_file");
        item.annotation = "The name of WT kernel file.";
        read_sync_string(of_kernel_file);
        this->add_item(item);
    }

    // 20. dft+u
    {
        Input_Item item("dft_plus_u");
        item.annotation = "new/old DFT+U correction method; 0: standard DFT calcullation(default)";
        read_sync_int(dft_plus_u);
        this->add_item(item);
    }
    {
        Input_Item item("yukawa_lambda");
        item.annotation = "default:0.0";
        read_sync_double(yukawa_lambda);
        this->add_item(item);
    }
    {
        Input_Item item("yukawa_potential");
        item.annotation = "default: false";
        read_sync_bool(yukawa_potential);
        this->add_item(item);
    }
    {
        Input_Item item("uramping");
        item.annotation = "increasing U values during SCF";
        read_sync_double(uramping);
        this->add_item(item);
    }
    {
        Input_Item item("omc");
        item.annotation = "the mode of occupation matrix control";
        read_sync_int(omc);
        this->add_item(item);
    }
    {
        Input_Item item("onsite_radius");
        item.annotation = "radius of the sphere for onsite projection (Bohr)";
        read_sync_double(onsite_radius);
        this->add_item(item);
    }
    {
        Input_Item item("hubbard_u");
        item.annotation = "Hubbard Coulomb interaction parameter U(ev)";
        item.readvalue = [](const Input_Item& item, Parameter& para) {
            int count = item.str_values.size();
            for (int i = 0; i < count; i++)
            {
                para.input.hubbard_u.push_back(convertstr<double>(item.str_values[i]));
            }
        };
        item.checkvalue = [](const Input_Item& item, const Parameter& para) {
            if (para.input.hubbard_u.size() != para.input.ntype)
            {
                throw std::runtime_error("hubbard_u should have the same number of elements as ntype");
            }
        };
        // We must firt bcast ntype (in item_general), then bcast hubbard_u
        sync_doublevec(hubbard_u, para.input.ntype);
        this->add_item(item);
    }
    {
        Input_Item item("orbital_corr");
        item.annotation = "which correlated orbitals need corrected ; d:2 ,f:3, do not need correction:-1";
        item.readvalue = [](const Input_Item& item, Parameter& para) {
            int count = item.str_values.size();
            for (int i = 0; i < count; i++)
            {
                para.input.orbital_corr.push_back(convertstr<int>(item.str_values[i]));
            }
        };
        item.checkvalue = [](const Input_Item& item, const Parameter& para) {
            if (para.input.orbital_corr.size() != para.input.ntype)
            {
                throw std::runtime_error("orbital_corr should have the same number of elements as ntype");
            }
        };
        // We must firt bcast ntype (in item_general), then bcast orbital_corr
        sync_intvec(orbital_corr, para.input.ntype);
        this->add_item(item);
    }

    // 21. spherical bessel
    {
        Input_Item item("bessel_nao_ecut");
        item.annotation = "energy cutoff for spherical bessel functions(Ry)";
        read_sync_string(bessel_nao_ecut);
        this->add_item(item);
    }
    {
        Input_Item item("bessel_nao_tolerence");
        item.annotation = "tolerence for spherical bessel root";
        read_sync_double(bessel_nao_tolerence);
        this->add_item(item);
    }
    {
        Input_Item item("bessel_nao_rcut");
        item.annotation = "radial cutoff for spherical bessel functions(a.u.)";
        item.readvalue = [](const Input_Item& item, Parameter& para) {
            int count = item.str_values.size();
            for (int i = 0; i < count; i++)
            {
                para.input.bessel_nao_rcuts.push_back(convertstr<double>(item.str_values[i]));
            }
            para.input.bessel_nao_rcut = para.input.bessel_nao_rcuts[0]; // also compatible with old input file
            para.input.nrcut = count;
        };
        // Since nrcut and bessel_nao_rcut are also valued, we need to broadcast them
        add_int_bcast(nrcut);
        add_double_bcast(bessel_nao_rcut);
        // We must firt bcast nrcut, then bcast bessel_nao_rcut
        sync_doublevec(bessel_nao_rcuts, para.input.nrcut);
        this->add_item(item);
    }
    {
        Input_Item item("bessel_nao_smooth");
        item.annotation = "spherical bessel smooth or not";
        read_sync_bool(bessel_nao_smooth);
        this->add_item(item);
    }
    {
        Input_Item item("bessel_nao_sigma");
        item.annotation = "spherical bessel smearing_sigma";
        read_sync_double(bessel_nao_sigma);
        this->add_item(item);
    }
    {
        Input_Item item("bessel_descriptor_lmax");
        item.annotation = "lmax used in generating spherical bessel functions";
        read_sync_int(bessel_descriptor_lmax);
        this->add_item(item);
    }
    {
        Input_Item item("bessel_descriptor_ecut");
        item.annotation = "energy cutoff for spherical bessel functions(Ry)";
        read_sync_string(bessel_descriptor_ecut);
        this->add_item(item);
    }
    {
        Input_Item item("bessel_descriptor_tolerence");
        item.annotation = "tolerence for spherical bessel root";
        read_sync_double(bessel_descriptor_tolerence);
        this->add_item(item);
    }
    {
        Input_Item item("bessel_descriptor_rcut");
        item.annotation = "radial cutoff for spherical bessel functions(a.u.)";
        read_sync_double(bessel_descriptor_rcut);
        this->add_item(item);
    }
    {
        Input_Item item("bessel_descriptor_smooth");
        item.annotation = "spherical bessel smooth or not";
        read_sync_bool(bessel_descriptor_smooth);
        this->add_item(item);
    }
    {
        Input_Item item("bessel_descriptor_sigma");
        item.annotation = "sphereical bessel smearing_sigma";
        read_sync_double(bessel_descriptor_sigma);
        this->add_item(item);
    }

    // 22. non-collinear spin-constrained
    {
        Input_Item item("sc_mag_switch");
        item.annotation = "switch to control spin-constrained DFT";
        read_sync_bool(sc_mag_switch);
        this->add_item(item);
    }
    {
        Input_Item item("decay_grad_switch");
        item.annotation = "switch to control gradient break condition";
        read_sync_bool(decay_grad_switch);
        this->add_item(item);
    }
    {
        Input_Item item("sc_thr");
        item.annotation = "Convergence criterion of spin-constrained iteration (RMS) in uB";
        read_sync_double(sc_thr);
        this->add_item(item);
    }
    {
        Input_Item item("nsc");
        item.annotation = "Maximal number of spin-constrained iteration";
        read_sync_int(nsc);
        this->add_item(item);
    }
    {
        Input_Item item("nsc_min");
        item.annotation = "Minimum number of spin-constrained iteration";
        read_sync_int(nsc_min);
        this->add_item(item);
    }
    {
        Input_Item item("sc_scf_nmin");
        item.annotation = "Minimum number of outer scf loop before initializing lambda loop";
        read_sync_int(sc_scf_nmin);
        this->add_item(item);
    }
    {
        Input_Item item("alpha_trial");
        item.annotation = "Initial trial step size for lambda in eV/uB^2";
        read_sync_double(alpha_trial);
        this->add_item(item);
    }
    {
        Input_Item item("sccut");
        item.annotation = "Maximal step size for lambda in eV/uB";
        read_sync_double(sccut);
        this->add_item(item);
    }
    {
        Input_Item item("sc_file");
        item.annotation = "file name for parameters used in non-collinear spin-constrained DFT (json format)";
        read_sync_string(sc_file);
        this->add_item(item);
    }

    // 23. Quasiatomic Orbital analysis
    {
        Input_Item item("qo_switch");
        item.annotation = "switch to control quasiatomic orbital analysis";
        read_sync_bool(qo_switch);
        this->add_item(item);
    }
    {
        Input_Item item("qo_basis");
        item.annotation
            = "type of QO basis function: hydrogen: hydrogen-like basis, pswfc: read basis from pseudopotential";
        read_sync_string(qo_basis);
        this->add_item(item);
    }
    {
        Input_Item item("qo_thr");
        item.annotation = "accuracy for evaluating cutoff radius of QO basis function";
        read_sync_double(qo_thr);
        this->add_item(item);
    }
    {
        Input_Item item("qo_strategy");
        item.annotation = "strategy to generate generate radial orbitals";
        item.readvalue = [](const Input_Item& item, Parameter& para) {
            int count = item.str_values.size();
            for (int i = 0; i < count; i++)
            {
                para.input.qo_strategy.push_back(item.str_values[i]);
            }
        };
        item.resetvalue = [](const Input_Item& item, Parameter& para) {
            if (para.input.qo_strategy.size() != para.input.ntype)
            {
                if (para.input.qo_strategy.size() == 1)
                {
                    para.input.qo_strategy.resize(para.input.ntype, para.input.qo_strategy[0]);
                }
                else
                {
                    std::string default_strategy;
                    if (para.input.qo_basis == "hydrogen")
                        default_strategy = "energy-valence";
                    else if ((para.input.qo_basis == "pswfc") || (para.input.qo_basis == "szv"))
                        default_strategy = "all";
                    else
                    {
                        throw std::runtime_error("When setting default values for qo_strategy, unexpected/unknown "
                                                 "qo_basis is found. Please check it.");
                    }
                    para.input.qo_strategy.resize(para.input.ntype, default_strategy);
                }
            }
        };
        // We must firt bcast ntype (in item_general), then bcast qo_strategy
        sync_stringvec(qo_strategy, para.input.ntype);
        this->add_item(item);
    };

    {
        Input_Item item("qo_screening_coeff");
        item.annotation = "rescale the shape of radial orbitals";
        item.readvalue = [](const Input_Item& item, Parameter& para) {
            int count = item.str_values.size();
            for (int i = 0; i < count; i++)
            {
                para.input.qo_screening_coeff.push_back(convertstr<double>(item.str_values[i]));
            }
        };
        item.resetvalue = [](const Input_Item& item, Parameter& para) {
            if (para.input.qo_screening_coeff.size() != para.input.ntype)
            {
                if (para.input.qo_basis == "pswfc")
                {
                    double default_screening_coeff
                        = (para.input.qo_screening_coeff.size() == 1) ? para.input.qo_screening_coeff[0] : 0.1;
                    para.input.qo_screening_coeff.resize(para.input.ntype, default_screening_coeff);
                }
                else
                {
                    throw std::runtime_error("qo_screening_coeff should have the same number of elements as ntype");
                }
            }
        };
        item.checkvalue = [](const Input_Item& item, const Parameter& para) {
            for (auto screen_coeff: para.input.qo_screening_coeff)
            {
                if (screen_coeff < 0)
                {
                    throw std::runtime_error("screening coefficient must >= 0 to tune the pswfc decay");
                }
                if (std::fabs(screen_coeff) < 1e-6)
                {
                    throw std::runtime_error(
                        "every low screening coefficient might yield very high computational cost");
                }
            }
        };
        // We must firt bcast ntype (in item_general), then bcast qo_screening_coeff
        sync_doublevec(qo_screening_coeff, para.input.ntype);
        this->add_item(item);
    }

    // 24. PEXSI
    {
        Input_Item item("pexsi_npole");
        item.annotation = "Number of poles in expansion";
        read_sync_int(pexsi_npole);
        this->add_item(item);
    }
    {
        Input_Item item("pexsi_inertia");
        item.annotation = "Whether inertia counting is used at the very beginning of PEXSI process";
        read_sync_bool(pexsi_inertia);
        this->add_item(item);
    }
    {
        Input_Item item("pexsi_nmax");
        item.annotation = "Maximum number of PEXSI iterations after each inertia counting procedure";
        read_sync_int(pexsi_nmax);
        this->add_item(item);
    }
    {
        Input_Item item("pexsi_comm");
        item.annotation = "Whether to construct PSelInv communication pattern";
        read_sync_bool(pexsi_comm);
        this->add_item(item);
    }
    {
        Input_Item item("pexsi_storage");
        item.annotation = "Storage space used by the Selected Inversion algorithm for symmetric matrices";
        read_sync_bool(pexsi_storage);
        this->add_item(item);
    }
    {
        Input_Item item("pexsi_ordering");
        item.annotation = "Ordering strategy for factorization and selected inversion";
        read_sync_int(pexsi_ordering);
        this->add_item(item);
    }
    {
        Input_Item item("pexsi_row_ordering");
        item.annotation
            = "Row permutation strategy for factorization and selected inversion, 0: NoRowPerm, 1: LargeDiag";
        read_sync_int(pexsi_row_ordering);
        this->add_item(item);
    }
    {
        Input_Item item("pexsi_nproc");
        item.annotation = "Number of processors for parmetis";
        read_sync_int(pexsi_nproc);
        this->add_item(item);
    }
    {
        Input_Item item("pexsi_symm");
        item.annotation = "Matrix symmetry";
        read_sync_bool(pexsi_symm);
        this->add_item(item);
    }
    {
        Input_Item item("pexsi_trans");
        item.annotation = "Whether to transpose";
        read_sync_bool(pexsi_trans);
        this->add_item(item);
    }
    {
        Input_Item item("pexsi_method");
        item.annotation = "pole expansion method, 1: Cauchy Contour Integral, 2: Moussa optimized method";
        read_sync_int(pexsi_method);
        this->add_item(item);
    }
    {
        Input_Item item("pexsi_nproc_pole");
        item.annotation = "Number of processes used by each pole";
        read_sync_int(pexsi_nproc_pole);
        this->add_item(item);
    }
    {
        Input_Item item("pexsi_temp");
        item.annotation = "Temperature, in the same unit as H";
        read_sync_double(pexsi_temp);
        this->add_item(item);
    }
    {
        Input_Item item("pexsi_gap");
        item.annotation = "Spectral gap";
        read_sync_double(pexsi_gap);
        this->add_item(item);
    }
    {
        Input_Item item("pexsi_delta_e");
        item.annotation = "An upper bound for the spectral radius of S^{-1} H";
        read_sync_double(pexsi_delta_e);
        this->add_item(item);
    }
    {
        Input_Item item("pexsi_mu_lower");
        item.annotation = "Initial guess of lower bound for mu";
        read_sync_double(pexsi_mu_lower);
        this->add_item(item);
    }
    {
        Input_Item item("pexsi_mu_upper");
        item.annotation = "Initial guess of upper bound for mu";
        read_sync_double(pexsi_mu_upper);
        this->add_item(item);
    }
    {
        Input_Item item("pexsi_mu");
        item.annotation = "Initial guess for mu (for the solver)";
        read_sync_double(pexsi_mu);
        this->add_item(item);
    }
    {
        Input_Item item("pexsi_mu_thr");
        item.annotation = "Stopping criterion in terms of the chemical potential for the inertia counting procedure";
        read_sync_double(pexsi_mu_thr);
        this->add_item(item);
    }
    {
        Input_Item item("pexsi_mu_expand");
        item.annotation = "If the chemical potential is not in the initial interval, the interval is expanded by "
                          "muInertiaExpansion";
        read_sync_double(pexsi_mu_expand);
        this->add_item(item);
    }
    {
        Input_Item item("pexsi_mu_guard");
        item.annotation
            = "Safe guard criterion in terms of the chemical potential to reinvoke the inertia counting procedure";
        read_sync_double(pexsi_mu_guard);
        this->add_item(item);
    }
    {
        Input_Item item("pexsi_elec_thr");
        item.annotation = "Stopping criterion of the PEXSI iteration in terms of the number of electrons compared to "
                          "numElectronExact";
        read_sync_double(pexsi_elec_thr);
        this->add_item(item);
    }
    {
        Input_Item item("pexsi_zero_thr");
        item.annotation = "if the absolute value of matrix element is less than ZERO_Limit, it will be considered as 0";
        read_sync_double(pexsi_zero_thr);
        this->add_item(item);
    }
}
} // namespace ModuleIO