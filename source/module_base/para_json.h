
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <vector>

#include <rapidjson/document.h>
#include <rapidjson/writer.h>
#include <rapidjson/stringbuffer.h>


/**
 * @brief   This namespace is used to store the json object of the 
 *          abacus parameter and its handlers. Used to read the parameters 
 *          at run time and finally organize them into json format files
 * 
*/
namespace Para_Json
{
   
    extern int test;
    // @param doc: the output json file
    extern rapidjson::Document doc;
    extern rapidjson::Value abacus;

    // @param general_info ：
    extern rapidjson::Value general_info;
    extern rapidjson::Value version;
    extern rapidjson::Value commit;
    extern rapidjson::Value begin_time;
    extern rapidjson::Value begin_time_epoch;
    extern rapidjson::Value device_g;
    // @param general_info -- parallel：
    extern rapidjson::Value parallel;
    extern rapidjson::Value drank;
    extern rapidjson::Value dsize;
    extern rapidjson::Value dcolor ;
    // @param general_info -- path
    extern rapidjson::Value path;
    extern rapidjson::Value global_out_dir;
    extern rapidjson::Value global_in_card;
    extern rapidjson::Value pseudo_dir_path ;
    extern rapidjson::Value orbital_dir_path;

    
    // @param reading_information：
    extern rapidjson::Value readin_info;
    // @param reading_information -- input_para：
    
    // @param reading_information -- input_para -- system_variables：
    extern rapidjson::Value system_variables;

    extern rapidjson::Value input_file;
    extern rapidjson::Value input_suffix;
    extern rapidjson::Value ntype;
    extern rapidjson::Value calculation;
    extern rapidjson::Value esolver_type;
    extern rapidjson::Value symmetry;
    extern rapidjson::Value symmetry_precfield;
    extern rapidjson::Value symmetry_autoclose;
    extern rapidjson::Value kpar;
    extern rapidjson::Value bndpar;
    extern rapidjson::Value latname;
    extern rapidjson::Value init_wfc;
    extern rapidjson::Value init_chg;
    extern rapidjson::Value init_vel;
    extern rapidjson::Value nelec;
    extern rapidjson::Value nupdown;
    extern rapidjson::Value dft_functional;
    extern rapidjson::Value xc_temperature;
    extern rapidjson::Value pseudo_rcut;
    extern rapidjson::Value pseudo_mesh;
    extern rapidjson::Value mem_saver;
    extern rapidjson::Value diago_proc;
    extern rapidjson::Value nbspline;
    extern rapidjson::Value kspacing;
    extern rapidjson::Value min_dist_coef;
    extern rapidjson::Value device;
    // @param reading_information -- input_para -- files_related

    extern rapidjson::Value stru_file;
    extern rapidjson::Value kpoint_file;
    extern rapidjson::Value pseudo_dir;
    extern rapidjson::Value orbital_dir;
    extern rapidjson::Value read_file_dir;
    extern rapidjson::Value wannier_card;
    // @param reading_information -- input_para -- planewave_related

    extern rapidjson::Value ecutwfc;
    extern rapidjson::Value nx;
    extern rapidjson::Value ny;
    extern rapidjson::Value nz;
    extern rapidjson::Value pw_seed;
    extern rapidjson::Value pw_diag_thr;
    extern rapidjson::Value pw_diag_nmax;
    extern rapidjson::Value pw_diag_ndim;
    // @param reading_information -- input_para -- numerical_atomic_orbitals_related
    
    extern rapidjson::Value nb2d;
    extern rapidjson::Value lmaxmax;
    extern rapidjson::Value lcao_ecut;
    extern rapidjson::Value lcao_dk;
    extern rapidjson::Value lcao_dr;
    extern rapidjson::Value lcao_rmax;
    extern rapidjson::Value search_radius;
    extern rapidjson::Value search_pbc;
    extern rapidjson::Value bx;
    extern rapidjson::Value by;
    extern rapidjson::Value bz;
    // @param reading_information -- input_para -- electronic_structure
    
    extern rapidjson::Value basis_type;
    extern rapidjson::Value ks_solver;
    extern rapidjson::Value nbands;
    extern rapidjson::Value nbands_istate;
    extern rapidjson::Value nspin;
    extern rapidjson::Value smearing_method;
    extern rapidjson::Value smearing_sigma;
    extern rapidjson::Value smearing_sigma_temp;
    extern rapidjson::Value mixing_type;
    extern rapidjson::Value mixing_beta;
    extern rapidjson::Value mixing_ndim;
    extern rapidjson::Value mixing_gg0;
    extern rapidjson::Value mixing_tau;
    extern rapidjson::Value mixing_dftu;
    extern rapidjson::Value gamma_only;
    extern rapidjson::Value printe;
    extern rapidjson::Value scf_nmax;
    extern rapidjson::Value scf_thr;
    extern rapidjson::Value scf_thr_type;
    extern rapidjson::Value chg_extrap;
    extern rapidjson::Value lspinorb;
    extern rapidjson::Value noncolin;
    extern rapidjson::Value soc_lambda;
    // @param reading_information -- input_para -- electronic_structure_SDFT

    extern rapidjson::Value method_sto;
    extern rapidjson::Value nbands_sto;
    extern rapidjson::Value nche_sto;
    extern rapidjson::Value emin_sto;
    extern rapidjson::Value emax_sto;
    extern rapidjson::Value seed_sto;
    extern rapidjson::Value initsto_freq;
    extern rapidjson::Value npart_sto;
    // @param reading_information -- input_para -- geometry_relaxation

    extern rapidjson::Value relax_method;
    extern rapidjson::Value relax_new;
    extern rapidjson::Value relax_scale_force;
    extern rapidjson::Value relax_nmax;
    extern rapidjson::Value relax_cg_thr;
    extern rapidjson::Value cal_force;
    extern rapidjson::Value force_thr;
    extern rapidjson::Value force_thr_ev;
    extern rapidjson::Value force_thr_ev2;
    extern rapidjson::Value relax_bfgs_w1;
    extern rapidjson::Value relax_bfgs_w2;
    extern rapidjson::Value relax_bfgs_rmax;
    extern rapidjson::Value relax_bfgs_rmin;
    extern rapidjson::Value relax_bfgs_init;
    extern rapidjson::Value cal_stress;
    extern rapidjson::Value stress_thr;
    extern rapidjson::Value press1;
    extern rapidjson::Value press2;
    extern rapidjson::Value press3;
    extern rapidjson::Value fixed_axes;
    extern rapidjson::Value fixed_ibrav;
    extern rapidjson::Value fixed_atoms;
    extern rapidjson::Value cell_factor;

    // @param reading_information -- input_para -- output_information_related

    extern rapidjson::Value out_mul;
    extern rapidjson::Value out_freq_elec;
    extern rapidjson::Value out_freq_ion;
    extern rapidjson::Value out_chg;
    extern rapidjson::Value out_pot;
    extern rapidjson::Value out_dm;
    extern rapidjson::Value out_dm1;
    extern rapidjson::Value out_wfc_pw;
    extern rapidjson::Value out_wfc_r;
    extern rapidjson::Value out_wfc_lcao;
    extern rapidjson::Value out_dos;
    extern rapidjson::Value out_band;
    extern rapidjson::Value out_proj_band;
    extern rapidjson::Value out_stru;
    extern rapidjson::Value out_bandgap;
    extern rapidjson::Value out_level;
    extern rapidjson::Value out_alllog;
    extern rapidjson::Value out_mat_hs;
    extern rapidjson::Value out_mat_r;
    extern rapidjson::Value out_mat_hs2;
    extern rapidjson::Value out_mat_t;
    extern rapidjson::Value out_mat_dh;
    extern rapidjson::Value out_app_flag;
    extern rapidjson::Value out_interval;
    extern rapidjson::Value out_element_info;
    extern rapidjson::Value restart_save;
    extern rapidjson::Value restart_load;
    extern rapidjson::Value rpa;

    // @param reading_information -- input_para -- density_of_states

    extern rapidjson::Value dos_edelta_ev;
    extern rapidjson::Value dos_sigma;
    extern rapidjson::Value dos_scale;
    extern rapidjson::Value dos_emin_ev;
    extern rapidjson::Value dos_emax_ev;
    extern rapidjson::Value dos_nche;
    // @param reading_information -- input_para -- naos
    extern rapidjson::Value bessel_nao_ecut;
    extern rapidjson::Value bessel_nao_tolerence;
    extern rapidjson::Value bessel_nao_rcut;
    extern rapidjson::Value bessel_nao_smooth;
    extern rapidjson::Value bessel_nao_sigma;
    // @param reading_information -- input_para -- deepks

    extern rapidjson::Value deepks_out_labels;
    extern rapidjson::Value deepks_scf;
    extern rapidjson::Value deepks_model;
    extern rapidjson::Value bessel_descriptor_lmax;
    extern rapidjson::Value bessel_descriptor_ecut;
    extern rapidjson::Value bessel_descriptor_tolerence;
    extern rapidjson::Value bessel_descriptor_rcut;
    extern rapidjson::Value bessel_descriptor_smooth;
    extern rapidjson::Value bessel_descriptor_sigma;
    extern rapidjson::Value deepks_bandgap;
    extern rapidjson::Value deepks_out_unittest;
    // @param reading_information -- input_para -- ofdft
    extern rapidjson::Value of_kinetic;
    extern rapidjson::Value of_method;
    extern rapidjson::Value of_conv;
    extern rapidjson::Value of_tole;
    extern rapidjson::Value of_tolp;
    extern rapidjson::Value of_tf_weight;
    extern rapidjson::Value of_vw_weight;
    extern rapidjson::Value of_wt_alpha;
    extern rapidjson::Value of_wt_beta;
    extern rapidjson::Value of_wt_rho0;
    extern rapidjson::Value of_hold_rho0;
    extern rapidjson::Value of_lkt_a;
    extern rapidjson::Value of_read_kernel;
    extern rapidjson::Value of_kernel_file;
    extern rapidjson::Value of_full_pw;
    extern rapidjson::Value of_full_pw_dim;

    // @param reading_information -- input_para -- electric_field_and_dipole_correction
    
    extern rapidjson::Value efield_flag;
    extern rapidjson::Value dip_cor_flag;
    extern rapidjson::Value efield_dir;
    extern rapidjson::Value efield_pos_max;
    extern rapidjson::Value efield_pos_dec;
    extern rapidjson::Value efield_amp;
    // @param reading_information -- input_para -- gate_field 
    
    extern rapidjson::Value gate_flag;
    extern rapidjson::Value zgate;
    extern rapidjson::Value block;
    extern rapidjson::Value block_down;
    extern rapidjson::Value block_up;
    extern rapidjson::Value block_height;
    // @param reading_information -- input_para -- exact_exchange
    extern rapidjson::Value exx_hybrid_alpha;
    extern rapidjson::Value exx_hse_omega;
    extern rapidjson::Value exx_separate_loop;
    extern rapidjson::Value exx_hybrid_step;
    extern rapidjson::Value exx_mixing_beta;
    extern rapidjson::Value exx_lambda;
    extern rapidjson::Value exx_pca_threshold;
    extern rapidjson::Value exx_c_threshold;
    extern rapidjson::Value exx_v_threshold;
    extern rapidjson::Value exx_dm_threshold;
    extern rapidjson::Value exx_c_grad_threshold;
    extern rapidjson::Value exx_v_grad_threshold;
    extern rapidjson::Value exx_schwarz_threshold;
    extern rapidjson::Value exx_cauchy_threshold;
    extern rapidjson::Value exx_cauchy_force_threshold;
    extern rapidjson::Value exx_cauchy_stress_threshold;
    extern rapidjson::Value exx_ccp_threshold;
    extern rapidjson::Value exx_ccp_rmesh_times;
    extern rapidjson::Value exx_distribute_type;
    extern rapidjson::Value exx_opt_orb_lmax;
    extern rapidjson::Value exx_opt_orb_ecut;
    extern rapidjson::Value exx_opt_orb_tolerence;
    extern rapidjson::Value exx_real_number;

    // @param reading_information -- input_para -- molecular_dynamics
    extern rapidjson::Value md_type;
    extern rapidjson::Value md_nstep;
    extern rapidjson::Value md_dt;
    extern rapidjson::Value md_thermostat;
    extern rapidjson::Value md_tlast;
    extern rapidjson::Value md_tfirst;
    extern rapidjson::Value md_restart;
    extern rapidjson::Value md_restartfreq;
    extern rapidjson::Value md_dumpfreq;
    extern rapidjson::Value dump_force;
    extern rapidjson::Value dump_vel;
    extern rapidjson::Value dump_virial;
    extern rapidjson::Value md_seed;
    extern rapidjson::Value md_tfreq;
    extern rapidjson::Value md_tchain;
    extern rapidjson::Value md_pmode;
    extern rapidjson::Value md_prec_level;
    extern rapidjson::Value ref_cell_factor;
    extern rapidjson::Value md_pcouple;
    extern rapidjson::Value md_pfirst;
    extern rapidjson::Value md_plast;
    extern rapidjson::Value md_pfreq;
    extern rapidjson::Value md_pchain;
    extern rapidjson::Value lj_rcut;
    extern rapidjson::Value lj_epsilon;
    extern rapidjson::Value lj_sigma;
    extern rapidjson::Value pot_file;
    extern rapidjson::Value msst_direction;
    extern rapidjson::Value msst_vel;
    extern rapidjson::Value msst_vis;
    extern rapidjson::Value msst_tscale;
    extern rapidjson::Value msst_qmass;
    extern rapidjson::Value md_damp;
    extern rapidjson::Value md_tolerance;
    extern rapidjson::Value md_nraise;
    extern rapidjson::Value cal_syns;
    extern rapidjson::Value dmax;

    // @param reading_information -- input_para -- dft_plus_u
    extern rapidjson::Value orbital_corr;
    extern rapidjson::Value hubbard_u;
    extern rapidjson::Value yukawa_potential;
    extern rapidjson::Value yukawa_lambda;
    extern rapidjson::Value omc;

    // @param reading_information -- input_para -- vdw_correction
    extern rapidjson::Value vdw_method;
    extern rapidjson::Value vdw_s6;
    extern rapidjson::Value vdw_s8;
    extern rapidjson::Value vdw_a1;
    extern rapidjson::Value vdw_a2;
    extern rapidjson::Value vdw_d;
    extern rapidjson::Value vdw_abc;
    extern rapidjson::Value vdw_C6_file;
    extern rapidjson::Value vdw_C6_unit;
    extern rapidjson::Value vdw_R0_file;
    extern rapidjson::Value vdw_R0_unit;
    extern rapidjson::Value vdw_cutoff_type;
    extern rapidjson::Value vdw_cutoff_radius;
    extern rapidjson::Value vdw_radius_unit;
    extern rapidjson::Value vdw_cutoff_period;
    extern rapidjson::Value vdw_cn_thr;
    extern rapidjson::Value vdw_cn_thr_unit;

    // @param reading_information -- input_para -- berry_phase_and_wannier90_interface
    extern rapidjson::Value berry_phase;
    extern rapidjson::Value gdir;
    extern rapidjson::Value towannier90;
    extern rapidjson::Value nnkpfile;
    extern rapidjson::Value wannier_spin;

    // @param reading_information -- input_para -- tddft
    extern rapidjson::Value td_edm;
    extern rapidjson::Value td_print_eij;
    extern rapidjson::Value td_propagator;
    extern rapidjson::Value td_vext;
    extern rapidjson::Value td_vext_dire;
    extern rapidjson::Value td_stype;
    extern rapidjson::Value td_ttype;
    extern rapidjson::Value td_tstart;
    extern rapidjson::Value td_tend;
    extern rapidjson::Value td_lcut1;
    extern rapidjson::Value td_lcut2;
    extern rapidjson::Value td_gauss_freq;
    extern rapidjson::Value td_gauss_phase;
    extern rapidjson::Value td_gauss_sigma;
    extern rapidjson::Value td_gauss_t0;
    extern rapidjson::Value td_gauss_amp;
    extern rapidjson::Value td_trape_freq;
    extern rapidjson::Value td_trape_phase;
    extern rapidjson::Value td_trape_t1;
    extern rapidjson::Value td_trape_t2;
    extern rapidjson::Value td_trape_t3;
    extern rapidjson::Value td_trape_amp;
    extern rapidjson::Value td_trigo_freq1;
    extern rapidjson::Value td_trigo_freq2;
    extern rapidjson::Value td_trigo_phase1;
    extern rapidjson::Value td_trigo_phase2;
    extern rapidjson::Value td_trigo_amp;
    extern rapidjson::Value td_heavi_t0;
    extern rapidjson::Value td_heavi_amp;
    extern rapidjson::Value td_out_dipole;
    extern rapidjson::Value td_out_efield;
    extern rapidjson::Value ocp;
    extern rapidjson::Value ocp_set;

    // @param reading_information -- input_para -- debuging_related
    extern rapidjson::Value t_in_h;
    extern rapidjson::Value vl_in_h;
    extern rapidjson::Value vnl_in_h;
    extern rapidjson::Value vh_in_h;
    extern rapidjson::Value vion_in_h;
    extern rapidjson::Value test_force;
    extern rapidjson::Value test_stress;
    extern rapidjson::Value colour;
    extern rapidjson::Value test_skip_ewald;

    // @param reading_information -- input_para -- electronic_conductivities
    extern rapidjson::Value cal_cond;
    extern rapidjson::Value cond_nche;
    extern rapidjson::Value cond_dw;
    extern rapidjson::Value cond_wcut;
    extern rapidjson::Value cond_dt;
    extern rapidjson::Value cond_dtbatch;
    extern rapidjson::Value cond_fwhm;
    extern rapidjson::Value cond_nonlocal;

    // @param reading_information -- input_para -- implicit_solvation_model
    extern rapidjson::Value imp_sol;
    extern rapidjson::Value eb_k;
    extern rapidjson::Value tau;
    extern rapidjson::Value sigma_k;
    extern rapidjson::Value nc_k;

    // @param reading_information -- stru_infos：
    extern rapidjson::Value stru_infos;
    // extern rapidjson::Value ATOMIC_SPECIES;
    // extern rapidjson::Value NUMERICAL_ORBITAL;
    // extern rapidjson::Value LATTICE_CONSTANT;
    // extern rapidjson::Value ATOMIC_POSITIONS;

    // @param reading_information -- KPT_infos
    extern rapidjson::Value KPT_infos;
    // extern rapidjson::Value total_number;
    // extern rapidjson::Value mode;
    // extern rapidjson::Value vectors;

    // @param reading_information -- orb_infos
    extern rapidjson::Value orb_infos;

    // @param reading_information -- pp
    extern rapidjson::Value pp;

    // @param init
    extern rapidjson::Value init;
    // @param init -- general
    // extern rapidjson::Value calculation;
    // extern rapidjson::Value esolver_type;
    // extern rapidjson::Value basis_type;
    // extern rapidjson::Value gamma_only;
    // extern rapidjson::Value ks_solver;
    // extern rapidjson::Value ntype;
    // extern rapidjson::Value nspin;
    // extern rapidjson::Value ecutwfc;
    // extern rapidjson::Value scf_thr;
    // extern rapidjson::Value scf_nmax;

    // @param init -- symmetry
    // extern rapidjson::Value symmetry;
    // extern rapidjson::Value BRAVAIS_TYPE;
    // extern rapidjson::Value BRAVAIS_LATTICE_NAME;
    // extern rapidjson::Value IBRAV;
    // extern rapidjson::Value LATTICE_CONSTANT_A;
    // extern rapidjson::Value right_hand_lattice;

    // @param init -- Kpoints
    extern rapidjson::Value kpoints;
    extern rapidjson::Value nkstot;
    extern rapidjson::Value nkstot_ibz;
    extern rapidjson::Value coordinates;
    extern rapidjson::Value weight;

    // @param init -- grid
    extern rapidjson::Value grid;
    extern rapidjson::Value energy_cutoff_for_wavefunc;
    extern rapidjson::Value fft_grid_for_wave_functions;
    extern rapidjson::Value number_of_plane_waves;
    extern rapidjson::Value number_of_sticks;

    // @param init -- Smearing
    // extern rapidjson::Value smearing_method;
    // extern rapidjson::Value smearing_sigma;

    // @param init -- mixing
    extern rapidjson::Value mixing;


    // @param output
    extern rapidjson::Value output;



    // @param final_stru
    extern rapidjson::Value final_stru;
    extern rapidjson::Value cell;
    extern rapidjson::Value coordinate;




    /**
     *  The functions below initialize the json output parameter 
     *  tree to connect the nodes of the module
    */

    /**
     * @brief   add Top stage：parameter in Abacus:
     */
    void Init_json_abacus();


    /**
     * @brief   add Second stage：parameter in Abacus - general_info:
     */
    void Init_json_abacus_generalInfo();


    /**
     * @brief   add Second stage：parameter in Abacus - readin_info:
     */
    void Init_json_abacus_readinInfo();


    /**
     * @brief   finish json tree build
     */
    void Finish_json_tree();



    /**
     * @brief   This function is used to populate the template type parameter 
     *          values into rapidjson's Value object
     */
    template <typename T> 
    void set_json_value(rapidjson::Value &json_v,T *para){
        if(std::is_same<T,int>::value)
        {
            json_v.SetInt(*reinterpret_cast<int*>(para)); 
        }
        else if(std::is_same<T,double>::value)
        {
            json_v.SetDouble(*reinterpret_cast<double*>(para));
        }
        else if(std::is_same<T,bool>::value)
        {
            json_v.SetBool(*reinterpret_cast<bool*>(para));
        }
        else if(std::is_same<T,std::string>::value)
        {
            // json_v.SetString(rapidjson::StringRef((*reinterpret_cast<std::string*>(para)).c_str()));

            json_v.SetString((*reinterpret_cast<std::string*>(para)).c_str(), std::strlen((*reinterpret_cast<std::string*>(para)).c_str()), doc.GetAllocator());
            //printf("exx_real_number = %s\n",(*reinterpret_cast<std::string*>(para)).c_str());
        }
    }
}

