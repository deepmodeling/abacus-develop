
#include "para_json.h"
#include "module_base/global_variable.h"


namespace Para_Json
{

    // @param doc: the output json file
    rapidjson::Document doc;
    rapidjson::Value abacus(rapidjson::kObjectType);
  
    // @param general_info ：
    rapidjson::Value general_info(rapidjson::kObjectType);
    rapidjson::Value version;
     
    rapidjson::Value commit;
    rapidjson::Value begin_time;
    rapidjson::Value begin_time_epoch;
    rapidjson::Value device_g;
    // @param general_info -- parallel：
    rapidjson::Value parallel(rapidjson::kObjectType);
    rapidjson::Value drank;
    rapidjson::Value dsize;
    rapidjson::Value dcolor ;
    // @param general_info -- path
    rapidjson::Value path(rapidjson::kObjectType);
    rapidjson::Value global_out_dir;
    rapidjson::Value global_in_card;
    rapidjson::Value pseudo_dir_path ;
    rapidjson::Value orbital_dir_path;

    
    // @param reading_information：
    rapidjson::Value readin_info(rapidjson::kObjectType);
    // @param reading_information -- input_file：
    rapidjson::Value input_file(rapidjson::kObjectType);


    // @param reading_information -- input_file -- system_variables：
    rapidjson::Value input_suffix;
    rapidjson::Value ntype;
    rapidjson::Value calculation;
    rapidjson::Value esolver_type;
    rapidjson::Value symmetry;
    rapidjson::Value symmetry_precfield;
    rapidjson::Value symmetry_autoclose;
    rapidjson::Value kpar;
    rapidjson::Value bndpar;
    rapidjson::Value latname;
    rapidjson::Value init_wfc;
    rapidjson::Value init_chg;
    rapidjson::Value init_vel;
    rapidjson::Value nelec;
    rapidjson::Value nupdown;
    rapidjson::Value dft_functional;
    rapidjson::Value xc_temperature;
    rapidjson::Value pseudo_rcut;
    rapidjson::Value pseudo_mesh;
    rapidjson::Value mem_saver;
    rapidjson::Value diago_proc;
    rapidjson::Value nbspline;
    rapidjson::Value kspacing;
    rapidjson::Value min_dist_coef;
    rapidjson::Value device;
    // @param reading_information -- input_file -- files_related
    rapidjson::Value stru_file;
    rapidjson::Value kpoint_file;
    rapidjson::Value pseudo_dir;
    rapidjson::Value orbital_dir;
    rapidjson::Value read_file_dir;
    rapidjson::Value wannier_card;

    // @param reading_information -- input_file -- planewave_related
    rapidjson::Value ecutwfc;
    rapidjson::Value nx;
    rapidjson::Value ny;
    rapidjson::Value nz;
    rapidjson::Value pw_seed;
    rapidjson::Value pw_diag_thr;
    rapidjson::Value pw_diag_nmax;
    rapidjson::Value pw_diag_ndim;
    // @param reading_information -- input_file -- numerical_atomic_orbitals_related
    rapidjson::Value nb2d;
    rapidjson::Value lmaxmax;
    rapidjson::Value lcao_ecut;
    rapidjson::Value lcao_dk;
    rapidjson::Value lcao_dr;
    rapidjson::Value lcao_rmax;
    rapidjson::Value search_radius;
    rapidjson::Value search_pbc;
    rapidjson::Value bx;
    rapidjson::Value by;
    rapidjson::Value bz;
    // @param reading_information -- input_file -- electronic_structure
    rapidjson::Value basis_type;
    rapidjson::Value ks_solver;
    rapidjson::Value nbands;
    rapidjson::Value nbands_istate;
    rapidjson::Value nspin;
    rapidjson::Value smearing_method;
    rapidjson::Value smearing_sigma;
    rapidjson::Value smearing_sigma_temp;
    rapidjson::Value mixing_type;
    rapidjson::Value mixing_beta;
    rapidjson::Value mixing_ndim;
    rapidjson::Value mixing_gg0;
    rapidjson::Value mixing_tau;
    rapidjson::Value mixing_dftu;
    rapidjson::Value gamma_only;
    rapidjson::Value printe;
    rapidjson::Value scf_nmax;
    rapidjson::Value scf_thr;
    rapidjson::Value scf_thr_type;
    rapidjson::Value chg_extrap;
    rapidjson::Value lspinorb;
    rapidjson::Value noncolin;
    rapidjson::Value soc_lambda;
    // @param reading_information -- input_file -- electronic_structure_SDFT
    rapidjson::Value method_sto;
    rapidjson::Value nbands_sto;
    rapidjson::Value nche_sto;
    rapidjson::Value emin_sto;
    rapidjson::Value emax_sto;
    rapidjson::Value seed_sto;
    rapidjson::Value initsto_freq;
    rapidjson::Value npart_sto;
    // @param reading_information -- input_file -- geometry_relaxation
    rapidjson::Value relax_method;
    rapidjson::Value relax_new;
    rapidjson::Value relax_scale_force;
    rapidjson::Value relax_nmax;
    rapidjson::Value relax_cg_thr;
    rapidjson::Value cal_force;
    rapidjson::Value force_thr;
    rapidjson::Value force_thr_ev;
    rapidjson::Value force_thr_ev2;
    rapidjson::Value relax_bfgs_w1;
    rapidjson::Value relax_bfgs_w2;
    rapidjson::Value relax_bfgs_rmax;
    rapidjson::Value relax_bfgs_rmin;
    rapidjson::Value relax_bfgs_init;
    rapidjson::Value cal_stress;
    rapidjson::Value stress_thr;
    rapidjson::Value press1;
    rapidjson::Value press2;
    rapidjson::Value press3;
    rapidjson::Value fixed_axes;
    rapidjson::Value fixed_ibrav;
    rapidjson::Value fixed_atoms;
    rapidjson::Value cell_factor;

    // @param reading_information -- input_file -- output_information_related
    rapidjson::Value out_mul;
    rapidjson::Value out_freq_elec;
    rapidjson::Value out_freq_ion;
    rapidjson::Value out_chg;
    rapidjson::Value out_pot;
    rapidjson::Value out_dm;
    rapidjson::Value out_dm1;
    rapidjson::Value out_wfc_pw;
    rapidjson::Value out_wfc_r;
    rapidjson::Value out_wfc_lcao;
    rapidjson::Value out_dos;
    rapidjson::Value out_band;
    rapidjson::Value out_proj_band;
    rapidjson::Value out_stru;
    rapidjson::Value out_bandgap;
    rapidjson::Value out_level;
    rapidjson::Value out_alllog;
    rapidjson::Value out_mat_hs;
    rapidjson::Value out_mat_r;
    rapidjson::Value out_mat_hs2;
    rapidjson::Value out_mat_t;
    rapidjson::Value out_mat_dh;
    rapidjson::Value out_app_flag;
    rapidjson::Value out_interval;
    rapidjson::Value out_element_info;
    rapidjson::Value restart_save;
    rapidjson::Value restart_load;
    rapidjson::Value rpa;

    // @param reading_information -- input_file -- density_of_states
    rapidjson::Value dos_edelta_ev;
    rapidjson::Value dos_sigma;
    rapidjson::Value dos_scale;
    rapidjson::Value dos_emin_ev;
    rapidjson::Value dos_emax_ev;
    rapidjson::Value dos_nche;
    // @param reading_information -- input_file -- naos
    rapidjson::Value bessel_nao_ecut;
    rapidjson::Value bessel_nao_tolerence;
    rapidjson::Value bessel_nao_rcut;
    rapidjson::Value bessel_nao_smooth;
    rapidjson::Value bessel_nao_sigma;
    // @param reading_information -- input_file -- deepks
    rapidjson::Value input_file_out_labels;
    rapidjson::Value input_file_scf;
    rapidjson::Value input_file_model;
    rapidjson::Value bessel_descriptor_lmax;
    rapidjson::Value bessel_descriptor_ecut;
    rapidjson::Value bessel_descriptor_tolerence;
    rapidjson::Value bessel_descriptor_rcut;
    rapidjson::Value bessel_descriptor_smooth;
    rapidjson::Value bessel_descriptor_sigma;
    rapidjson::Value input_file_bandgap;
    rapidjson::Value input_file_out_unittest;
    // @param reading_information -- input_file -- ofdft
    rapidjson::Value of_kinetic;
    rapidjson::Value of_method;
    rapidjson::Value of_conv;
    rapidjson::Value of_tole;
    rapidjson::Value of_tolp;
    rapidjson::Value of_tf_weight;
    rapidjson::Value of_vw_weight;
    rapidjson::Value of_wt_alpha;
    rapidjson::Value of_wt_beta;
    rapidjson::Value of_wt_rho0;
    rapidjson::Value of_hold_rho0;
    rapidjson::Value of_lkt_a;
    rapidjson::Value of_read_kernel;
    rapidjson::Value of_kernel_file;
    rapidjson::Value of_full_pw;
    rapidjson::Value of_full_pw_dim;

    // @param reading_information -- input_file -- electric_field_and_dipole_correction
    rapidjson::Value efield_flag;
    rapidjson::Value dip_cor_flag;
    rapidjson::Value efield_dir;
    rapidjson::Value efield_pos_max;
    rapidjson::Value efield_pos_dec;
    rapidjson::Value efield_amp;
    // @param reading_information -- input_file -- gate_field 
    rapidjson::Value gate_flag;
    rapidjson::Value zgate;
    rapidjson::Value block;
    rapidjson::Value block_down;
    rapidjson::Value block_up;
    rapidjson::Value block_height;
    // @param reading_information -- input_file -- exact_exchange
    rapidjson::Value exx_hybrid_alpha;
    rapidjson::Value exx_hse_omega;
    rapidjson::Value exx_separate_loop;
    rapidjson::Value exx_hybrid_step;
    rapidjson::Value exx_mixing_beta;
    rapidjson::Value exx_lambda;
    rapidjson::Value exx_pca_threshold;
    rapidjson::Value exx_c_threshold;
    rapidjson::Value exx_v_threshold;
    rapidjson::Value exx_dm_threshold;
    rapidjson::Value exx_c_grad_threshold;
    rapidjson::Value exx_v_grad_threshold;
    rapidjson::Value exx_schwarz_threshold;
    rapidjson::Value exx_cauchy_threshold;
    rapidjson::Value exx_cauchy_force_threshold;
    rapidjson::Value exx_cauchy_stress_threshold;
    rapidjson::Value exx_ccp_threshold;
    rapidjson::Value exx_ccp_rmesh_times;
    rapidjson::Value exx_distribute_type;
    rapidjson::Value exx_opt_orb_lmax;
    rapidjson::Value exx_opt_orb_ecut;
    rapidjson::Value exx_opt_orb_tolerence;
    rapidjson::Value exx_real_number;

    // @param reading_information -- input_file -- molecular_dynamics
    rapidjson::Value md_type;
    rapidjson::Value md_nstep;
    rapidjson::Value md_dt;
    rapidjson::Value md_thermostat;
    rapidjson::Value md_tlast;
    rapidjson::Value md_tfirst;
    rapidjson::Value md_restart;
    rapidjson::Value md_restartfreq;
    rapidjson::Value md_dumpfreq;
    rapidjson::Value dump_force;
    rapidjson::Value dump_vel;
    rapidjson::Value dump_virial;
    rapidjson::Value md_seed;
    rapidjson::Value md_tfreq;
    rapidjson::Value md_tchain;
    rapidjson::Value md_pmode;
    rapidjson::Value md_prec_level;
    rapidjson::Value ref_cell_factor;
    rapidjson::Value md_pcouple;
    rapidjson::Value md_pfirst;
    rapidjson::Value md_plast;
    rapidjson::Value md_pfreq;
    rapidjson::Value md_pchain;
    rapidjson::Value lj_rcut;
    rapidjson::Value lj_epsilon;
    rapidjson::Value lj_sigma;
    rapidjson::Value pot_file;
    rapidjson::Value msst_direction;
    rapidjson::Value msst_vel;
    rapidjson::Value msst_vis;
    rapidjson::Value msst_tscale;
    rapidjson::Value msst_qmass;
    rapidjson::Value md_damp;
    rapidjson::Value md_tolerance;
    rapidjson::Value md_nraise;
    rapidjson::Value cal_syns;
    rapidjson::Value dmax;

    // @param reading_information -- input_file -- dft_plus_u
    rapidjson::Value orbital_corr;
    rapidjson::Value hubbard_u;
    rapidjson::Value yukawa_potential;
    rapidjson::Value yukawa_lambda;
    rapidjson::Value omc;

    // @param reading_information -- input_file -- vdw_correction
    rapidjson::Value vdw_method;
    rapidjson::Value vdw_s6;
    rapidjson::Value vdw_s8;
    rapidjson::Value vdw_a1;
    rapidjson::Value vdw_a2;
    rapidjson::Value vdw_d;
    rapidjson::Value vdw_abc;
    rapidjson::Value vdw_C6_file;
    rapidjson::Value vdw_C6_unit;
    rapidjson::Value vdw_R0_file;
    rapidjson::Value vdw_R0_unit;
    rapidjson::Value vdw_cutoff_type;
    rapidjson::Value vdw_cutoff_radius;
    rapidjson::Value vdw_radius_unit;
    rapidjson::Value vdw_cutoff_period;
    rapidjson::Value vdw_cn_thr;
    rapidjson::Value vdw_cn_thr_unit;

    // @param reading_information -- input_file -- berry_phase_and_wannier90_interface
    rapidjson::Value berry_phase;
    rapidjson::Value gdir;
    rapidjson::Value towannier90;
    rapidjson::Value nnkpfile;
    rapidjson::Value wannier_spin;

    // @param reading_information -- input_file -- tddft
    rapidjson::Value td_edm;
    rapidjson::Value td_print_eij;
    rapidjson::Value td_propagator;
    rapidjson::Value td_vext;
    rapidjson::Value td_vext_dire;
    rapidjson::Value td_stype;
    rapidjson::Value td_ttype;
    rapidjson::Value td_tstart;
    rapidjson::Value td_tend;
    rapidjson::Value td_lcut1;
    rapidjson::Value td_lcut2;
    rapidjson::Value td_gauss_freq;
    rapidjson::Value td_gauss_phase;
    rapidjson::Value td_gauss_sigma;
    rapidjson::Value td_gauss_t0;
    rapidjson::Value td_gauss_amp;
    rapidjson::Value td_trape_freq;
    rapidjson::Value td_trape_phase;
    rapidjson::Value td_trape_t1;
    rapidjson::Value td_trape_t2;
    rapidjson::Value td_trape_t3;
    rapidjson::Value td_trape_amp;
    rapidjson::Value td_trigo_freq1;
    rapidjson::Value td_trigo_freq2;
    rapidjson::Value td_trigo_phase1;
    rapidjson::Value td_trigo_phase2;
    rapidjson::Value td_trigo_amp;
    rapidjson::Value td_heavi_t0;
    rapidjson::Value td_heavi_amp;
    rapidjson::Value td_out_dipole;
    rapidjson::Value td_out_efield;
    rapidjson::Value ocp;
    rapidjson::Value ocp_set;

    // @param reading_information -- input_file -- debuging_related
    rapidjson::Value t_in_h;
    rapidjson::Value vl_in_h;
    rapidjson::Value vnl_in_h;
    rapidjson::Value vh_in_h;
    rapidjson::Value vion_in_h;
    rapidjson::Value test_force;
    rapidjson::Value test_stress;
    rapidjson::Value colour;
    rapidjson::Value test_skip_ewald;

    // @param reading_information -- input_file -- electronic_conductivities
    rapidjson::Value cal_cond;
    rapidjson::Value cond_nche;
    rapidjson::Value cond_dw;
    rapidjson::Value cond_wcut;
    rapidjson::Value cond_dt;
    rapidjson::Value cond_dtbatch;
    rapidjson::Value cond_fwhm;
    rapidjson::Value cond_nonlocal;
    // @param reading_information -- input_file -- implicit_solvation_model
    rapidjson::Value imp_sol;
    rapidjson::Value eb_k;
    rapidjson::Value tau;
    rapidjson::Value sigma_k;
    rapidjson::Value nc_k;

    // @param reading_information -- stru_infos：
    rapidjson::Value stru_infos(rapidjson::kObjectType);
    // rapidjson::Value ATOMIC_SPECIES(rapidjson::kArrayType);
    // rapidjson::Value NUMERICAL_ORBITAL;
    // rapidjson::Value LATTICE_CONSTANT(rapidjson::kArrayType);
    // rapidjson::Value ATOMIC_POSITIONS(rapidjson::kArrayType);

    // @param reading_information -- KPT_infos
    rapidjson::Value KPT_infos(rapidjson::kObjectType);
    // rapidjson::Value total_number;
    // rapidjson::Value mode;
    // rapidjson::Value vectors(rapidjson::kArrayType);

    // @param reading_information -- orb_infos
    rapidjson::Value orb_infos(rapidjson::kObjectType);

    // @param reading_information -- pp
    rapidjson::Value pp(rapidjson::kObjectType);

    // @param init
    rapidjson::Value init(rapidjson::kObjectType);
    // @param init -- general
    // rapidjson::Value calculation;
    // rapidjson::Value esolver_type;
    // rapidjson::Value basis_type;
    // rapidjson::Value gamma_only;
    // rapidjson::Value ks_solver;
    // rapidjson::Value ntype;
    // rapidjson::Value nspin;
    // rapidjson::Value ecutwfc;
    // rapidjson::Value scf_thr;
    // rapidjson::Value scf_nmax;

    // @param init -- symmetry
    // rapidjson::Value symmetry(rapidjson::kObjectType);
    // rapidjson::Value BRAVAIS_TYPE;
    // rapidjson::Value BRAVAIS_LATTICE_NAME;
    // rapidjson::Value IBRAV;
    // rapidjson::Value LATTICE_CONSTANT_A;
    // rapidjson::Value right_hand_lattice;

    // @param init -- Kpoints
    rapidjson::Value kpoints(rapidjson::kObjectType);
    rapidjson::Value nkstot;
    rapidjson::Value nkstot_ibz;
    rapidjson::Value coordinates(rapidjson::kArrayType);
    rapidjson::Value weight(rapidjson::kArrayType);

    // @param init -- grid
    rapidjson::Value grid(rapidjson::kObjectType);
    rapidjson::Value energy_cutoff_for_wavefunc;
    rapidjson::Value fft_grid_for_wave_functions(rapidjson::kArrayType);
    rapidjson::Value number_of_plane_waves;
    rapidjson::Value number_of_sticks;

    // @param init -- Smearing
    // rapidjson::Value smearing_method;
    // rapidjson::Value smearing_sigma;

    // @param init -- mixing
    rapidjson::Value mixing;


    // @param output
    rapidjson::Value output(rapidjson::kArrayType);



    // @param final_stru
    rapidjson::Value final_stru(rapidjson::kObjectType);
    rapidjson::Value cell;
    rapidjson::Value coordinate;



    /**
     *  The functions below initialize the json output parameter 
     *  tree to connect the nodes of the module
    */

    /**
     * @brief   add Top stage：parameter in Abacus:
     */
    void Init_json_abacus()
    {


        // add First stage：parameter in abcus:

        abacus.AddMember("general_info", general_info, doc.GetAllocator());

        abacus.AddMember("readin_info", readin_info, doc.GetAllocator());
        
        abacus.AddMember("init", init, doc.GetAllocator());

        abacus.AddMember("output", output, doc.GetAllocator());

        abacus.AddMember("final_stru", final_stru, doc.GetAllocator());

        doc.SetObject();
        // abacus.SetObject();
        doc.AddMember("ABACUS", abacus, doc.GetAllocator());
        /**
         * .
         * .
         * .
         * .
         * .
         * .
         * .
         * */
    }
    /**
     * @brief   add Second stage：parameter in Abacus - general_info:
     */
    void Init_json_abacus_generalInfo(){
        general_info.AddMember("version", version, doc.GetAllocator());

        general_info.AddMember("commit", commit, doc.GetAllocator());      

        general_info.AddMember("begin_time", begin_time, doc.GetAllocator());      

        general_info.AddMember("begin_date", begin_date, doc.GetAllocator());     

        general_info.AddMember("device", device_g, doc.GetAllocator());                


        
        parallel.AddMember("drank", drank, doc.GetAllocator());

        parallel.AddMember("dsize", dsize, doc.GetAllocator());
                        
        parallel.AddMember("dcolor", dcolor, doc.GetAllocator());
    

        // add Third stage：parameter in parallel:
        general_info.AddMember("parallel", parallel, doc.GetAllocator());
                
    }
    /**
     * @brief   delete null node 
     */
    void RemoveNullValues(rapidjson::Value& parent) {
        if (parent.IsObject()) {
            for (rapidjson::Value::MemberIterator itr = parent.MemberBegin(); itr != parent.MemberEnd(); ) {
                if (itr->value.IsNull()) {
                    itr = parent.EraseMember(itr);
                } else {
                    // delet son null node
                    RemoveNullValues(itr->value);
                    ++itr;
                }
            }
        } else if (parent.IsArray()) {
            for (int i = 0; i < parent.Size(); ) {
                if (parent[i].IsNull()) {
                    parent.Erase(parent.Begin() + i);
                } else {
                    // delet son null node
                    RemoveNullValues(parent[i]);
                    ++i;
                }
            }
        }
    }

    /**
     * @brief   add Second stage：parameter in Abacus - readin_info:
     */
    void Init_json_abacus_readinInfo(){
        //add Third stage：parameter in system_variables:
        input_file.AddMember("suffix", input_suffix, doc.GetAllocator());
        input_file.AddMember("ntype", ntype, doc.GetAllocator());
        input_file.AddMember("calculation", calculation, doc.GetAllocator());
        input_file.AddMember("esolver_type", esolver_type, doc.GetAllocator());
        input_file.AddMember("symmetry", symmetry, doc.GetAllocator());
        input_file.AddMember("symmetry_precfield", symmetry_precfield, doc.GetAllocator());
        input_file.AddMember("symmetry_autoclose", symmetry_autoclose, doc.GetAllocator());
        input_file.AddMember("kpar", kpar, doc.GetAllocator());
        input_file.AddMember("bndpar", bndpar, doc.GetAllocator());
        input_file.AddMember("latname", latname, doc.GetAllocator());
        input_file.AddMember("init_wfc", init_wfc, doc.GetAllocator());
        input_file.AddMember("init_chg", init_chg, doc.GetAllocator());
        input_file.AddMember("init_vel", init_vel, doc.GetAllocator());
        input_file.AddMember("nelec", nelec, doc.GetAllocator());
        input_file.AddMember("nupdown", nupdown, doc.GetAllocator());
        input_file.AddMember("dft_functional", dft_functional, doc.GetAllocator());
        input_file.AddMember("xc_temperature", xc_temperature, doc.GetAllocator());
        input_file.AddMember("pseudo_rcut", pseudo_rcut, doc.GetAllocator());
        input_file.AddMember("pseudo_mesh", pseudo_mesh, doc.GetAllocator());
        input_file.AddMember("mem_saver", mem_saver, doc.GetAllocator());
        input_file.AddMember("diago_proc", diago_proc, doc.GetAllocator());
        input_file.AddMember("nbspline", nbspline, doc.GetAllocator());
        input_file.AddMember("kspacing", kspacing, doc.GetAllocator());
        input_file.AddMember("min_dist_coef", min_dist_coef, doc.GetAllocator());
        input_file.AddMember("device", device, doc.GetAllocator());

        //add Third stage：parameter in files_related:
        input_file.AddMember("stru_file", stru_file, doc.GetAllocator());
        input_file.AddMember("kpoint_file", kpoint_file, doc.GetAllocator());
        input_file.AddMember("pseudo_dir", pseudo_dir, doc.GetAllocator());
        input_file.AddMember("orbital_dir", orbital_dir, doc.GetAllocator());
        input_file.AddMember("read_file_dir", read_file_dir, doc.GetAllocator());
        input_file.AddMember("wannier_card", wannier_card, doc.GetAllocator());
    
        //add Third stage：parameter in planewave_related:
        input_file.AddMember("ecutwfc", ecutwfc, doc.GetAllocator());
        input_file.AddMember("nx", nx, doc.GetAllocator());
        input_file.AddMember("ny", ny, doc.GetAllocator());
        input_file.AddMember("nz", nz, doc.GetAllocator());
        input_file.AddMember("pw_seed", pw_seed, doc.GetAllocator());
        input_file.AddMember("pw_diag_thr", pw_diag_thr, doc.GetAllocator());
        input_file.AddMember("pw_diag_nmax", pw_diag_nmax, doc.GetAllocator());
        input_file.AddMember("pw_diag_ndim", pw_diag_ndim, doc.GetAllocator());    
    
    
        //add Third stage：parameter in numerical_atomic_orbitals_related:
        input_file.AddMember("nb2d", nb2d, doc.GetAllocator());
        input_file.AddMember("lmaxmax", lmaxmax, doc.GetAllocator());
        input_file.AddMember("lcao_ecut", lcao_ecut, doc.GetAllocator());
        input_file.AddMember("lcao_dk", lcao_dk, doc.GetAllocator());
        input_file.AddMember("lcao_dr", lcao_dr, doc.GetAllocator());
        input_file.AddMember("lcao_rmax", lcao_rmax, doc.GetAllocator());
        input_file.AddMember("search_radius", search_radius, doc.GetAllocator());
        input_file.AddMember("search_pbc", search_pbc, doc.GetAllocator());
        input_file.AddMember("bx", bx, doc.GetAllocator());
        input_file.AddMember("by", by, doc.GetAllocator());
        input_file.AddMember("bz", bz, doc.GetAllocator());        
    
        //add Third stage：parameter in electronic_structure:
        input_file.AddMember("basis_type", basis_type, doc.GetAllocator());
        input_file.AddMember("ks_solver", ks_solver, doc.GetAllocator());
        input_file.AddMember("nbands", nbands, doc.GetAllocator());
        input_file.AddMember("nbands_istate", nbands_istate, doc.GetAllocator());
        input_file.AddMember("nspin", nspin, doc.GetAllocator());
        input_file.AddMember("smearing_method", smearing_method, doc.GetAllocator());
        input_file.AddMember("smearing_sigma", smearing_sigma, doc.GetAllocator());
        input_file.AddMember("smearing_sigma_temp", smearing_sigma_temp, doc.GetAllocator());
        input_file.AddMember("mixing_type", mixing_type, doc.GetAllocator());
        input_file.AddMember("mixing_beta", mixing_beta, doc.GetAllocator());
        input_file.AddMember("mixing_ndim", mixing_ndim, doc.GetAllocator());
        input_file.AddMember("mixing_gg0", mixing_gg0, doc.GetAllocator());
        input_file.AddMember("mixing_tau", mixing_tau, doc.GetAllocator());
        input_file.AddMember("mixing_dftu", mixing_dftu, doc.GetAllocator());
        input_file.AddMember("gamma_only", gamma_only, doc.GetAllocator());
        input_file.AddMember("printe", printe, doc.GetAllocator());
        input_file.AddMember("scf_nmax", scf_nmax, doc.GetAllocator());
        input_file.AddMember("scf_thr", scf_thr, doc.GetAllocator());
        input_file.AddMember("scf_thr_type", scf_thr_type, doc.GetAllocator());
        input_file.AddMember("chg_extrap", chg_extrap, doc.GetAllocator());
        input_file.AddMember("lspinorb", lspinorb, doc.GetAllocator());
        input_file.AddMember("noncolin", noncolin, doc.GetAllocator());
        input_file.AddMember("soc_lambda", soc_lambda, doc.GetAllocator());    


        //add Third stage：parameter in electronic_structure_SDFT:
        input_file.AddMember("method_sto", method_sto, doc.GetAllocator());
        input_file.AddMember("nbands_sto", nbands_sto, doc.GetAllocator());
        input_file.AddMember("nche_sto", nche_sto, doc.GetAllocator());
        input_file.AddMember("emin_sto", emin_sto, doc.GetAllocator());
        input_file.AddMember("emax_sto", emax_sto, doc.GetAllocator());
        input_file.AddMember("seed_sto", seed_sto, doc.GetAllocator());
        input_file.AddMember("initsto_freq", initsto_freq, doc.GetAllocator());
        input_file.AddMember("npart_sto", npart_sto, doc.GetAllocator());
        
        
        //add Third stage：parameter in geometry_relaxation:
        input_file.AddMember("relax_method", relax_method, doc.GetAllocator());
        input_file.AddMember("relax_new", relax_new, doc.GetAllocator());
        input_file.AddMember("relax_scale_force", relax_scale_force, doc.GetAllocator());
        input_file.AddMember("relax_nmax", relax_nmax, doc.GetAllocator());
        input_file.AddMember("relax_cg_thr", relax_cg_thr, doc.GetAllocator());
        input_file.AddMember("cal_force", cal_force, doc.GetAllocator());
        input_file.AddMember("force_thr", force_thr, doc.GetAllocator());
        input_file.AddMember("force_thr_ev", force_thr_ev, doc.GetAllocator());
        input_file.AddMember("force_thr_ev2", force_thr_ev2, doc.GetAllocator());
        input_file.AddMember("relax_bfgs_w1", relax_bfgs_w1, doc.GetAllocator());
        input_file.AddMember("relax_bfgs_w2", relax_bfgs_w2, doc.GetAllocator());
        input_file.AddMember("relax_bfgs_rmax", relax_bfgs_rmax, doc.GetAllocator());
        input_file.AddMember("relax_bfgs_rmin", relax_bfgs_rmin, doc.GetAllocator());
        input_file.AddMember("relax_bfgs_init", relax_bfgs_init, doc.GetAllocator());
        input_file.AddMember("cal_stress", cal_stress, doc.GetAllocator());
        input_file.AddMember("stress_thr", stress_thr, doc.GetAllocator());
        input_file.AddMember("press1", press1, doc.GetAllocator());
        input_file.AddMember("press2", press2, doc.GetAllocator());
        input_file.AddMember("press3", press3, doc.GetAllocator());
        input_file.AddMember("fixed_axes", fixed_axes, doc.GetAllocator());
        input_file.AddMember("fixed_ibrav", fixed_ibrav, doc.GetAllocator());
        input_file.AddMember("fixed_atoms", fixed_atoms, doc.GetAllocator());
        input_file.AddMember("cell_factor", cell_factor, doc.GetAllocator());
        
        
        //add Third stage：parameter in output_information_related:
        input_file.AddMember("out_mul", out_mul, doc.GetAllocator());
        input_file.AddMember("out_freq_elec", out_freq_elec, doc.GetAllocator());
        input_file.AddMember("out_freq_ion", out_freq_ion, doc.GetAllocator());        
        input_file.AddMember("out_chg", out_chg, doc.GetAllocator());
        input_file.AddMember("out_pot", out_pot, doc.GetAllocator());
        input_file.AddMember("out_dm", out_dm, doc.GetAllocator());
        input_file.AddMember("out_dm1", out_dm1, doc.GetAllocator());
        input_file.AddMember("out_wfc_pw", out_wfc_pw, doc.GetAllocator());
        input_file.AddMember("out_wfc_r", out_wfc_r, doc.GetAllocator());
        input_file.AddMember("out_wfc_lcao", out_wfc_lcao, doc.GetAllocator());
        input_file.AddMember("out_dos", out_dos, doc.GetAllocator());
        input_file.AddMember("out_band", out_band, doc.GetAllocator());
        input_file.AddMember("out_proj_band", out_proj_band, doc.GetAllocator());
        input_file.AddMember("out_stru", out_stru, doc.GetAllocator());
        input_file.AddMember("out_bandgap", out_bandgap, doc.GetAllocator());
        input_file.AddMember("out_level", out_level, doc.GetAllocator());
        input_file.AddMember("out_alllog", out_alllog, doc.GetAllocator());
        input_file.AddMember("out_mat_hs", out_mat_hs, doc.GetAllocator());
        input_file.AddMember("out_mat_r", out_mat_r, doc.GetAllocator());
        input_file.AddMember("out_mat_hs2", out_mat_hs2, doc.GetAllocator());
        input_file.AddMember("out_mat_t", out_mat_t, doc.GetAllocator());
        input_file.AddMember("out_mat_dh", out_mat_dh, doc.GetAllocator());
        input_file.AddMember("out_app_flag", out_app_flag, doc.GetAllocator());
        input_file.AddMember("out_interval", out_interval, doc.GetAllocator());
        input_file.AddMember("out_element_info", out_element_info, doc.GetAllocator());
        input_file.AddMember("restart_save", restart_save, doc.GetAllocator());
        input_file.AddMember("restart_load", restart_load, doc.GetAllocator());
        input_file.AddMember("rpa", rpa, doc.GetAllocator());

        //add Third stage：parameter in density_of_states:
        input_file.AddMember("dos_edelta_ev", dos_edelta_ev, doc.GetAllocator());
        input_file.AddMember("dos_sigma", dos_sigma, doc.GetAllocator());
        input_file.AddMember("dos_scale", dos_scale, doc.GetAllocator());
        input_file.AddMember("dos_emin_ev", dos_emin_ev, doc.GetAllocator());
        input_file.AddMember("dos_emax_ev", dos_emax_ev, doc.GetAllocator());
        input_file.AddMember("dos_nche", dos_nche, doc.GetAllocator());
        
        //add Third stage：parameter in naos:
        input_file.AddMember("bessel_nao_ecut", bessel_nao_ecut, doc.GetAllocator());
        input_file.AddMember("bessel_nao_tolerence", bessel_nao_tolerence, doc.GetAllocator());
        input_file.AddMember("bessel_nao_rcut", bessel_nao_rcut, doc.GetAllocator());
        input_file.AddMember("bessel_nao_smooth", bessel_nao_smooth, doc.GetAllocator());
        input_file.AddMember("bessel_nao_sigma", bessel_nao_sigma, doc.GetAllocator());
        
        //add Third stage：parameter in deepks:
        input_file.AddMember("input_file_out_labels", input_file_out_labels, doc.GetAllocator());
        input_file.AddMember("input_file_scf", input_file_scf, doc.GetAllocator());
        input_file.AddMember("input_file_model", input_file_model, doc.GetAllocator());
        input_file.AddMember("bessel_descriptor_lmax", bessel_descriptor_lmax, doc.GetAllocator());
        input_file.AddMember("bessel_descriptor_ecut", bessel_descriptor_ecut, doc.GetAllocator());
        input_file.AddMember("bessel_descriptor_tolerence", bessel_descriptor_tolerence, doc.GetAllocator());
        input_file.AddMember("bessel_descriptor_rcut", bessel_descriptor_rcut, doc.GetAllocator());
        input_file.AddMember("bessel_descriptor_smooth", bessel_descriptor_smooth, doc.GetAllocator());
        input_file.AddMember("bessel_descriptor_sigma", bessel_descriptor_sigma, doc.GetAllocator());
        input_file.AddMember("input_file_bandgap", input_file_bandgap, doc.GetAllocator());
        input_file.AddMember("input_file_out_unittest", input_file_out_unittest, doc.GetAllocator());
        
        //add Third stage：parameter in ofdft:
        input_file.AddMember("of_kinetic", of_kinetic, doc.GetAllocator());
        input_file.AddMember("of_method", of_method, doc.GetAllocator());
        input_file.AddMember("of_conv", of_conv, doc.GetAllocator());
        input_file.AddMember("of_tole", of_tole, doc.GetAllocator());
        input_file.AddMember("of_tolp", of_tolp, doc.GetAllocator());
        input_file.AddMember("of_tf_weight", of_tf_weight, doc.GetAllocator());
        input_file.AddMember("of_vw_weight", of_vw_weight, doc.GetAllocator());
        input_file.AddMember("of_wt_alpha", of_wt_alpha, doc.GetAllocator());
        input_file.AddMember("of_wt_beta", of_wt_beta, doc.GetAllocator());
        input_file.AddMember("of_wt_rho0", of_wt_rho0, doc.GetAllocator());
        input_file.AddMember("of_hold_rho0", of_hold_rho0, doc.GetAllocator());
        input_file.AddMember("of_lkt_a", of_lkt_a, doc.GetAllocator());
        input_file.AddMember("of_read_kernel", of_read_kernel, doc.GetAllocator());
        input_file.AddMember("of_kernel_file", of_kernel_file, doc.GetAllocator());
        input_file.AddMember("of_full_pw", of_full_pw, doc.GetAllocator());
        input_file.AddMember("of_full_pw_dim", of_full_pw_dim, doc.GetAllocator());
        
        
        //add Third stage：parameter in electric_field_and_dipole_correction:
        input_file.AddMember("efield_flag", efield_flag, doc.GetAllocator());
        input_file.AddMember("dip_cor_flag", dip_cor_flag, doc.GetAllocator());
        input_file.AddMember("efield_dir", efield_dir, doc.GetAllocator());
        input_file.AddMember("efield_pos_max", efield_pos_max, doc.GetAllocator());
        input_file.AddMember("efield_pos_dec", efield_pos_dec, doc.GetAllocator());
        input_file.AddMember("efield_amp", efield_amp, doc.GetAllocator());
        
        //add Third stage：parameter in gate_field:
        input_file.AddMember("gate_flag", gate_flag, doc.GetAllocator());
        input_file.AddMember("zgate", zgate, doc.GetAllocator());
        input_file.AddMember("block", block, doc.GetAllocator());
        input_file.AddMember("block_down", block_down, doc.GetAllocator());
        input_file.AddMember("block_up", block_up, doc.GetAllocator());
        input_file.AddMember("block_height", block_height, doc.GetAllocator());
    
        //add Third stage：parameter in exact_exchange:
        input_file.AddMember("exx_hybrid_alpha", exx_hybrid_alpha, doc.GetAllocator());
        input_file.AddMember("exx_hse_omega", exx_hse_omega, doc.GetAllocator());
        input_file.AddMember("exx_separate_loop", exx_separate_loop, doc.GetAllocator());
        input_file.AddMember("exx_hybrid_step", exx_hybrid_step, doc.GetAllocator());
        input_file.AddMember("exx_mixing_beta", exx_mixing_beta, doc.GetAllocator());
        input_file.AddMember("exx_lambda", exx_lambda, doc.GetAllocator());
        input_file.AddMember("exx_pca_threshold", exx_pca_threshold, doc.GetAllocator());
        input_file.AddMember("exx_c_threshold", exx_c_threshold, doc.GetAllocator());
        input_file.AddMember("exx_v_threshold", exx_v_threshold, doc.GetAllocator());
        input_file.AddMember("exx_dm_threshold", exx_dm_threshold, doc.GetAllocator());
        input_file.AddMember("exx_c_grad_threshold", exx_c_grad_threshold, doc.GetAllocator());
        input_file.AddMember("exx_v_grad_threshold", exx_v_grad_threshold, doc.GetAllocator());
        input_file.AddMember("exx_schwarz_threshold", exx_schwarz_threshold, doc.GetAllocator());
        input_file.AddMember("exx_cauchy_threshold", exx_cauchy_threshold, doc.GetAllocator());
        input_file.AddMember("exx_cauchy_force_threshold", exx_cauchy_force_threshold, doc.GetAllocator());
        input_file.AddMember("exx_cauchy_stress_threshold", exx_cauchy_stress_threshold, doc.GetAllocator());
        input_file.AddMember("exx_ccp_threshold", exx_ccp_threshold, doc.GetAllocator());
        input_file.AddMember("exx_ccp_rmesh_times", exx_ccp_rmesh_times, doc.GetAllocator());
        input_file.AddMember("exx_distribute_type", exx_distribute_type, doc.GetAllocator());
        input_file.AddMember("exx_opt_orb_lmax", exx_opt_orb_lmax, doc.GetAllocator());
        input_file.AddMember("exx_opt_orb_ecut", exx_opt_orb_ecut, doc.GetAllocator());
        input_file.AddMember("exx_opt_orb_tolerence", exx_opt_orb_tolerence, doc.GetAllocator());
        input_file.AddMember("exx_real_number", exx_real_number, doc.GetAllocator());
        
        
        //add Third stage：parameter in molecular_dynamics:
        input_file.AddMember("md_type", md_type, doc.GetAllocator());
        input_file.AddMember("md_nstep", md_nstep, doc.GetAllocator());
        input_file.AddMember("md_dt", md_dt, doc.GetAllocator());
        input_file.AddMember("md_thermostat", md_thermostat, doc.GetAllocator());
        input_file.AddMember("md_tlast", md_tlast, doc.GetAllocator());
        input_file.AddMember("md_tfirst", md_tfirst, doc.GetAllocator());
        input_file.AddMember("md_restart", md_restart, doc.GetAllocator());
        input_file.AddMember("md_restartfreq", md_restartfreq, doc.GetAllocator());
        input_file.AddMember("md_dumpfreq", md_dumpfreq, doc.GetAllocator());
        input_file.AddMember("dump_force", dump_force, doc.GetAllocator());
        input_file.AddMember("dump_vel", dump_vel, doc.GetAllocator());
        input_file.AddMember("dump_virial", dump_virial, doc.GetAllocator());
        input_file.AddMember("md_seed", md_seed, doc.GetAllocator());
        input_file.AddMember("md_tfreq", md_tfreq, doc.GetAllocator());
        input_file.AddMember("md_tchain", md_tchain, doc.GetAllocator());
        input_file.AddMember("md_pmode", md_pmode, doc.GetAllocator());
        input_file.AddMember("md_prec_level", md_prec_level, doc.GetAllocator());
        input_file.AddMember("ref_cell_factor", ref_cell_factor, doc.GetAllocator());
        input_file.AddMember("md_pcouple", md_pcouple, doc.GetAllocator());
        input_file.AddMember("md_pfirst", md_pfirst, doc.GetAllocator());
        input_file.AddMember("md_plast", md_plast, doc.GetAllocator());
        input_file.AddMember("md_pfreq", md_pfreq, doc.GetAllocator());
        input_file.AddMember("md_pchain", md_pchain, doc.GetAllocator());
        input_file.AddMember("lj_rcut", lj_rcut, doc.GetAllocator());
        input_file.AddMember("lj_epsilon", lj_epsilon, doc.GetAllocator());
        input_file.AddMember("lj_sigma", lj_sigma, doc.GetAllocator());
        input_file.AddMember("pot_file", pot_file, doc.GetAllocator());
        input_file.AddMember("msst_direction", msst_direction, doc.GetAllocator());
        input_file.AddMember("msst_vel", msst_vel, doc.GetAllocator());
        input_file.AddMember("msst_vis", msst_vis, doc.GetAllocator());
        input_file.AddMember("msst_tscale", msst_tscale, doc.GetAllocator());
        input_file.AddMember("msst_qmass", msst_qmass, doc.GetAllocator());
        input_file.AddMember("md_damp", md_damp, doc.GetAllocator());
        input_file.AddMember("md_tolerance", md_tolerance, doc.GetAllocator());
        input_file.AddMember("md_nraise", md_nraise, doc.GetAllocator());
        input_file.AddMember("cal_syns", cal_syns, doc.GetAllocator());
        input_file.AddMember("dmax", dmax, doc.GetAllocator());

        //add Third stage：parameter in dft_plus_u:
        input_file.AddMember("orbital_corr", orbital_corr, doc.GetAllocator());
        input_file.AddMember("hubbard_u", hubbard_u, doc.GetAllocator());
        input_file.AddMember("yukawa_potential", yukawa_potential, doc.GetAllocator());
        input_file.AddMember("yukawa_lambda", yukawa_lambda, doc.GetAllocator());
        input_file.AddMember("omc", omc, doc.GetAllocator());

        //add Third stage：parameter in vdw_correction:
        input_file.AddMember("vdw_method", vdw_method, doc.GetAllocator());
        input_file.AddMember("vdw_s6", vdw_s6, doc.GetAllocator());
        input_file.AddMember("vdw_s8", vdw_s8, doc.GetAllocator());
        input_file.AddMember("vdw_a1", vdw_a1, doc.GetAllocator());
        input_file.AddMember("vdw_a2", vdw_a2, doc.GetAllocator());
        input_file.AddMember("vdw_d", vdw_d, doc.GetAllocator());
        input_file.AddMember("vdw_abc", vdw_abc, doc.GetAllocator());
        input_file.AddMember("vdw_C6_file", vdw_C6_file, doc.GetAllocator());
        input_file.AddMember("vdw_C6_unit", vdw_C6_unit, doc.GetAllocator());
        input_file.AddMember("vdw_R0_file", vdw_R0_file, doc.GetAllocator());
        input_file.AddMember("vdw_R0_unit", vdw_R0_unit, doc.GetAllocator());
        input_file.AddMember("vdw_cutoff_type", vdw_cutoff_type, doc.GetAllocator());
        input_file.AddMember("vdw_cutoff_radius", vdw_cutoff_radius, doc.GetAllocator());
        input_file.AddMember("vdw_radius_unit", vdw_radius_unit, doc.GetAllocator());
        input_file.AddMember("vdw_cutoff_period", vdw_cutoff_period, doc.GetAllocator());
        input_file.AddMember("vdw_cn_thr", vdw_cn_thr, doc.GetAllocator());
        input_file.AddMember("vdw_cn_thr_unit", vdw_cn_thr_unit, doc.GetAllocator());

        //add Third stage：parameter in berry_phase_and_wannier90_interface:
        input_file.AddMember("berry_phase", berry_phase, doc.GetAllocator());
        input_file.AddMember("gdir", gdir, doc.GetAllocator());
        input_file.AddMember("towannier90", towannier90, doc.GetAllocator());
        input_file.AddMember("nnkpfile", nnkpfile, doc.GetAllocator());
        input_file.AddMember("wannier_spin", wannier_spin, doc.GetAllocator());    
    
        //add Third stage：parameter in tddft:
        input_file.AddMember("td_edm", td_edm, doc.GetAllocator());
        input_file.AddMember("td_print_eij", td_print_eij, doc.GetAllocator());
        input_file.AddMember("td_propagator", td_propagator, doc.GetAllocator());
        input_file.AddMember("td_vext", td_vext, doc.GetAllocator());
        input_file.AddMember("td_vext_dire", td_vext_dire, doc.GetAllocator());
        input_file.AddMember("td_stype", td_stype, doc.GetAllocator());
        input_file.AddMember("td_ttype", td_ttype, doc.GetAllocator());
        input_file.AddMember("td_tstart", td_tstart, doc.GetAllocator());
        input_file.AddMember("td_tend", td_tend, doc.GetAllocator());
        input_file.AddMember("td_lcut1", td_lcut1, doc.GetAllocator());
        input_file.AddMember("td_lcut2", td_lcut2, doc.GetAllocator());
        input_file.AddMember("td_gauss_freq", td_gauss_freq, doc.GetAllocator());
        input_file.AddMember("td_gauss_phase", td_gauss_phase, doc.GetAllocator());
        input_file.AddMember("td_gauss_sigma", td_gauss_sigma, doc.GetAllocator());
        input_file.AddMember("td_gauss_t0", td_gauss_t0, doc.GetAllocator());
        input_file.AddMember("td_gauss_amp", td_gauss_amp, doc.GetAllocator());
        input_file.AddMember("td_trape_freq", td_trape_freq, doc.GetAllocator());
        input_file.AddMember("td_trape_phase", td_trape_phase, doc.GetAllocator());
        input_file.AddMember("td_trape_t1", td_trape_t1, doc.GetAllocator());
        input_file.AddMember("td_trape_t2", td_trape_t2, doc.GetAllocator());
        input_file.AddMember("td_trape_t3", td_trape_t3, doc.GetAllocator());
        input_file.AddMember("td_trape_amp", td_trape_amp, doc.GetAllocator());
        input_file.AddMember("td_trigo_freq1", td_trigo_freq1, doc.GetAllocator());
        input_file.AddMember("td_trigo_freq2", td_trigo_freq2, doc.GetAllocator());
        input_file.AddMember("td_trigo_phase1", td_trigo_phase1, doc.GetAllocator());
        input_file.AddMember("td_trigo_phase2", td_trigo_phase2, doc.GetAllocator());
        input_file.AddMember("td_trigo_amp", td_trigo_amp, doc.GetAllocator());
        input_file.AddMember("td_heavi_t0", td_heavi_t0, doc.GetAllocator());
        input_file.AddMember("td_heavi_amp", td_heavi_amp, doc.GetAllocator());
        input_file.AddMember("td_out_dipole", td_out_dipole, doc.GetAllocator());
        input_file.AddMember("td_out_efield", td_out_efield, doc.GetAllocator());
        input_file.AddMember("ocp", ocp, doc.GetAllocator());
        input_file.AddMember("ocp_set", ocp_set, doc.GetAllocator());

        //add Third stage：parameter in debuging_related:
        input_file.AddMember("t_in_h", t_in_h, doc.GetAllocator());
        input_file.AddMember("vl_in_h", vl_in_h, doc.GetAllocator());
        input_file.AddMember("vnl_in_h", vnl_in_h, doc.GetAllocator());
        input_file.AddMember("vh_in_h", vh_in_h, doc.GetAllocator());
        input_file.AddMember("vion_in_h", vion_in_h, doc.GetAllocator());
        input_file.AddMember("test_force", test_force, doc.GetAllocator());
        input_file.AddMember("test_stress", test_stress, doc.GetAllocator());
        input_file.AddMember("colour", colour, doc.GetAllocator());
        input_file.AddMember("test_skip_ewald", test_skip_ewald, doc.GetAllocator());

        //add Third stage：parameter in electronic_conductivities:
        input_file.AddMember("cal_cond", cal_cond, doc.GetAllocator());
        input_file.AddMember("cond_nche", cond_nche, doc.GetAllocator());
        input_file.AddMember("cond_dw", cond_dw, doc.GetAllocator());
        input_file.AddMember("cond_wcut", cond_wcut, doc.GetAllocator());
        input_file.AddMember("cond_dt", cond_dt, doc.GetAllocator());
        input_file.AddMember("cond_dtbatch", cond_dtbatch, doc.GetAllocator());
        input_file.AddMember("cond_fwhm", cond_fwhm, doc.GetAllocator());
        input_file.AddMember("cond_nonlocal", cond_nonlocal, doc.GetAllocator());

        //add Third stage：parameter in implicit_solvation_model:
        input_file.AddMember("imp_sol", imp_sol, doc.GetAllocator());
        input_file.AddMember("eb_k", eb_k, doc.GetAllocator());
        input_file.AddMember("tau", tau, doc.GetAllocator());
        input_file.AddMember("sigma_k", sigma_k, doc.GetAllocator());
        input_file.AddMember("nc_k", nc_k, doc.GetAllocator());


        RemoveNullValues(input_file);


        // after add child_node's node in readin_info, add child node
        // add parameters in readin_info:
        readin_info.AddMember("input_file", input_file, doc.GetAllocator());

    }
    std::string getLevelStr(int level)
    {
        std::string levelStr = "";
        for (int i = 0; i < level; i++)
        {
            levelStr += "\t"; 
        }
        return levelStr;

    }

    // format a json string
    std::string formatJson(std::string json)
    {
        std::string result = "";
        int level = 0;
        for (std::string::size_type index = 0; index < json.size(); index++)
        {
            char c = json[index];

            if (level > 0 && '\n' == json[json.size() - 1])
            {
                result += getLevelStr(level);
            }

            switch (c)
            {
            case '{':
            case '[':
                result = result + c + "\n";
                level++;
                result += getLevelStr(level);
                break;
            case ',':
                result = result + c + "\n";
                result += getLevelStr(level);
                break;
            case '}':
            case ']':
                result += "\n";
                level--;
                result += getLevelStr(level);
                result += c;
                break;
            default:
                result += c;
                break;
            }

        }
        return result;
    }


    // @brief   finish json tree and output to json file
    void Finish_json_tree(){

        // init json tree
        Init_json_abacus_readinInfo();
        Init_json_abacus_generalInfo();
        Init_json_abacus();

        
        // Converts a json object to a string
        rapidjson::StringBuffer buffer;
        rapidjson::Writer<rapidjson::StringBuffer> writer(buffer);
        doc.Accept(writer);

        // Output the json string to a file
        std::string json_path;
        json_path.append("abacus.json");

        // format json
        std::string strJson= buffer.GetString();
        std::string beauterJson = formatJson(strJson);

        std::ofstream ofs(json_path);
        ofs << beauterJson << std::endl;
        ofs.close();
    }




}