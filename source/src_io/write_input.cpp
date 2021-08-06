#include "../input.h"
#include "../src_pw/tools.h"

void Input::Print(const string &fn)const
{
    if (GlobalV::MY_RANK!=0) return;

    TITLE("Input","Print");

    ofstream ofs(fn.c_str());

	//----------------------------------
	// output the information in INPUT.
	//----------------------------------
    ofs << "INPUT_PARAMETERS" << endl;
	ofs << setiosflags(ios::left);
	
	ofs << "#Parameters (1.General)" << endl;
	OUTP(ofs,"suffix",suffix,"the name of main output directory");
	OUTP(ofs,"latname",latname,"the name of lattice name");
	OUTP(ofs,"atom_file",GlobalV::global_atom_card,"the filename of file containing atom positions");
	OUTP(ofs,"kpoint_file",GlobalV::global_kpoint_card,"the name of file containing k points");
	OUTP(ofs,"pseudo_dir",GlobalV::global_pseudo_dir,"the directory containing pseudo files");
	OUTP(ofs,"pseudo_type",GlobalV::global_pseudo_type,"the type pseudo files");
	OUTP(ofs,"pseudo_rcut",pseudo_rcut,"cut-off radius for radial integration");
	OUTP(ofs,"renormwithmesh",renormwithmesh,"0: use our own mesh to do radial renormalization; 1: use mesh as in QE");
	OUTP(ofs,"lmaxmax",lmaxmax,"maximum of l channels used");
	OUTP(ofs,"dft_functional",dft_functional,"exchange correlation functional");
	OUTP(ofs,"calculation",calculation,"test; scf; relax; nscf; ienvelope; istate;");
	OUTP(ofs,"ntype",ntype,"atom species number");
	OUTP(ofs,"nspin",nspin,"1: single spin; 2: up and down spin; 4: noncollinear spin");
	OUTP(ofs,"nbands",nbands,"number of bands");
	OUTP(ofs,"nbands_sto",nbands_sto,"number of stochastic bands");
	OUTP(ofs,"nbands_istate",nbands_istate,"number of bands around Fermi level for istate calulation");
	OUTP(ofs,"nche_sto",nche_sto,"number of orders for Chebyshev expansion in stochastic DFT");
	OUTP(ofs,"symmetry",symmetry,"turn symmetry on or off");	
	OUTP(ofs,"set_vel",set_vel,"read velocity from STRU or not");
	OUTP(ofs,"nelec",nelec,"input number of electrons");
	OUTP(ofs,"tot_magnetization",tot_magnetization,"total magnetization of the system");

	ofs << "\n#Parameters (2.PW)" << endl;
	OUTP(ofs,"ecutwfc",ecutwfc,"#energy cutoff for wave functions");
	if(ks_solver=="cg")
	{
		OUTP(ofs,"diago_cg_maxiter",diago_cg_maxiter,"max iteration number for cg");
		OUTP(ofs,"diago_cg_prec",diago_cg_prec,"diago_cg_prec");
	}
	else if(ks_solver=="dav")
	{
		OUTP(ofs,"diago_david_ndim",diago_david_ndim,"max dimension for davidson");
	}
	OUTP(ofs,"ethr",ethr,"threshold for eigenvalues is cg electron iterations");
	OUTP(ofs,"dr2",dr2,"charge density error");
	OUTP(ofs,"start_wfc",start_wfc,"start wave functions are from 'atomic' or 'file'");
	OUTP(ofs,"start_charge",start_pot,"start charge is from 'atomic' or file");
	OUTP(ofs,"charge_extrap",charge_extrap,"atomic; first-order; second-order; dm:coefficients of SIA");
	OUTP(ofs,"out_charge",out_charge,">0 output charge density for selected electron steps");
	OUTP(ofs,"out_potential",out_potential,"output realspace potential");
	OUTP(ofs,"out_wf",out_wf,"output wave functions");
	OUTP(ofs,"out_dos",out_dos,"output energy and dos");
	OUTP(ofs,"out_band",out_band,"output energy and band structure");
	OUTP(ofs,"restart_save",restart_save,"print to disk every step for restart");
	OUTP(ofs,"restart_load",restart_load,"restart from disk");
	OUTP(ofs,"read_file_dir",read_file_dir,"directory of files for reading");
	OUTP(ofs,"nx",nx,"number of points along x axis for FFT grid");
	OUTP(ofs,"ny",ny,"number of points along y axis for FFT grid");
	OUTP(ofs,"nz",nz,"number of points along z axis for FFT grid");	
	OUTP(ofs,"cell_factor",cell_factor,"used in the construction of the pseudopotential tables");	
	
	ofs << "\n#Parameters (3.Relaxation)" << endl;
	OUTP(ofs,"ks_solver",GlobalV::KS_SOLVER,"cg; dav; lapack; genelpa; hpseps; scalapack_gvx");
	OUTP(ofs,"niter",niter,"#number of electron iterations");
	OUTP(ofs,"force_set",force_set,"output the force_set or not"); 
	OUTP(ofs,"nstep",nstep,"number of ion iteration steps");
	OUTP(ofs,"out_stru",out_stru,"output the structure files after each ion step");
	OUTP(ofs,"force_thr",force_thr,"force threshold, unit: Ry/Bohr");
	OUTP(ofs,"force_thr_ev",force_thr*13.6058/0.529177,"force threshold, unit: eV/Angstrom");
	OUTP(ofs,"force_thr_ev2",force_thr_ev2,"force invalid threshold, unit: eV/Angstrom");
	OUTP(ofs,"cg_threshold",cg_threshold,"threshold for switching from cg to bfgs, unit: eV/Angstrom");
	OUTP(ofs,"stress_thr",stress_thr,"stress threshold");
	OUTP(ofs,"press1",press1,"target pressure, unit: KBar");
	OUTP(ofs,"press2",press2,"target pressure, unit: KBar");
	OUTP(ofs,"press3",press3,"target pressure, unit: KBar");
	OUTP(ofs,"bfgs_w1",bfgs_w1,"wolfe condition 1 for bfgs");
	OUTP(ofs,"bfgs_w2",bfgs_w2,"wolfe condition 2 for bfgs");
	OUTP(ofs,"trust_radius_max", trust_radius_max,"maximal trust radius, unit: Bohr");
	OUTP(ofs,"trust_radius_min", trust_radius_min,"minimal trust radius, unit: Bohr");
	OUTP(ofs,"trust_radius_ini", trust_radius_ini,"initial trust radius, unit: Bohr");
	OUTP(ofs,"stress",stress,"calculate the stress or not");
	OUTP(ofs,"fixed_axes",fixed_axes,"which axes are fixed");
	OUTP(ofs,"move_method",ion_dynamics,"bfgs; sd; cg; cg_bfgs;"); //pengfei add 2013-08-15
	OUTP(ofs,"out_level",out_level,"ie(for electrons); i(for ions);");
	OUTP(ofs,"out_dm",out_dm,">0 output density matrix");

	// for deepks
	OUTP(ofs,"out_descriptor",out_descriptor,">0 compute descriptor for deepks");
	OUTP(ofs,"lmax_descriptor",lmax_descriptor,">0 lmax used in descriptor for deepks");

	ofs << "\n#Parameters (4.LCAO)" << endl;
	OUTP(ofs,"basis_type",basis_type,"PW; LCAO in pw; LCAO");
	OUTP(ofs,"new_dm",new_dm,"Type of density matrix; 0: old 1: new");
	if(ks_solver=="HPSEPS" || ks_solver=="genelpa" || ks_solver=="scalapack_gvx")
	{
		OUTP(ofs,"nb2d",nb2d,"2d distribution of atoms");
	}
	OUTP(ofs,"search_radius",search_radius,"input search radius (Bohr)");
	OUTP(ofs,"search_pbc",search_pbc,"input periodic boundary condition");
	OUTP(ofs,"lcao_ecut",lcao_ecut,"energy cutoff for LCAO");
	OUTP(ofs,"lcao_dk",lcao_dk,"delta k for 1D integration in LCAO");
	OUTP(ofs,"lcao_dr",lcao_dr,"delta r for 1D integration in LCAO");
	OUTP(ofs,"lcao_rmax",lcao_rmax,"max R for 1D two-center integration table");
	OUTP(ofs,"out_hs",out_hs,"output H and S matrix");
	OUTP(ofs,"out_hs2",out_hs2,"output H(R) and S(R) matrix");
	OUTP(ofs,"out_r",out_r_matrix,"output r(R) matrix");
	OUTP(ofs,"out_lowf",out_lowf,"ouput LCAO wave functions");
	OUTP(ofs,"bx",bx,"division of an element grid in FFT grid along x");
	OUTP(ofs,"by",by,"division of an element grid in FFT grid along y");
	OUTP(ofs,"bz",bz,"division of an element grid in FFT grid along z");

	ofs << "\n#Parameters (5.Smearing)" << endl;
	OUTP(ofs,"smearing",smearing,"type of smearing: gauss; fd; fixed; mp; mp2; mv");
	OUTP(ofs,"sigma",degauss,"energy range for smearing");
	
	ofs << "\n#Parameters (6.Charge Mixing)" << endl;
	OUTP(ofs,"mixing_type",mixing_mode,"plain; kerker; pulay; pulay-kerker; broyden");
	OUTP(ofs,"mixing_beta",mixing_beta,"mixing parameter: 0 means no new charge");
	OUTP(ofs,"mixing_ndim",mixing_ndim,"mixing dimension in pulay");
	OUTP(ofs,"mixing_gg0",mixing_gg0,"mixing parameter in kerker");

	ofs << "\n#Parameters (7.DOS)" << endl;
	OUTP(ofs,"dos_emin_ev",dos_emin_ev,"minimal range for dos");
	OUTP(ofs,"dos_emax_ev",dos_emax_ev,"maximal range for dos");
	OUTP(ofs,"dos_edelta_ev",dos_edelta_ev,"delta energy for dos");
	OUTP(ofs,"dos_scale",dos_scale,"scale dos range by");
	OUTP(ofs,"dos_sigma",b_coef,"gauss b coefficeinet(default=0.07)");

	ofs << "\n#Parameters (8.Technique)" << endl;
	OUTP(ofs,"gamma_only",gamma_only,"gamma only, only used in LCAO basis");
	OUTP(ofs,"diago_proc",GlobalV::DIAGO_PROC,"number of proc used to diago");
	OUTP(ofs,"npool",npool,"number of pools for k points, pw only");
	OUTP(ofs,"mem_saver",mem_saver,"memory saver for many k points used");
	OUTP(ofs,"printe",printe,"print band energy for selectively ionic steps");

	ofs << "\n#Parameters (9.SIAO)" << endl;
	OUTP(ofs,"selinv_npole",selinv_npole,"number of selected poles");
	OUTP(ofs,"selinv_temp",selinv_temp,"temperature for Fermi-Dirac distribution");
	OUTP(ofs,"selinv_gap",selinv_gap,"supposed gap in the calculation");
	OUTP(ofs,"selinv_deltae",selinv_deltae,"expected energy range");
	OUTP(ofs,"selinv_mu",selinv_mu,"chosen mu as Fermi energy");
	OUTP(ofs,"selinv_threshold",selinv_threshold,"threshold for calculated electron number");
	OUTP(ofs,"selinv_niter",selinv_niter,"max number of steps to update mu");

	ofs << "\n#Parameters (10.Molecular dynamics)" << endl;
	OUTP(ofs,"md_mdtype",mdp.mdtype,"choose ensemble");
	//OUTP(ofs,"md_potential",mdp.md_potential,"choose potential for md");
	OUTP(ofs,"md_dt",mdp.dt,"time step");
	OUTP(ofs,"mnhc",mdp.MNHC,"number of Nose-Hoover chains");
	OUTP(ofs,"md_qmass",mdp.Qmass,"mass of thermostat");
	OUTP(ofs,"md_tfirst",mdp.tfirst,"temperature first");
	OUTP(ofs,"md_tlast",mdp.tlast,"temperature last");
	OUTP(ofs,"md_dumpmdfred",mdp.recordFreq,"The period to dump MD information for monitoring and restarting MD");
	OUTP(ofs,"md_mdoutpath",mdp.mdoutputpath,"output path of md");
	OUTP(ofs,"md_rstmd",mdp.rstMD,"whether restart");
	OUTP(ofs,"md_fixtemperature",mdp.fixTemperature,"period to change temperature");
	OUTP(ofs,"md_ediff",mdp.ediff,"parameter for constraining total energy change");
	OUTP(ofs,"md_ediffg",mdp.ediffg,"parameter for constraining max force change");
	OUTP(ofs,"NVT_tau",mdp.NVT_tau,"parameter for adjust effect of thermostat");
	OUTP(ofs,"NVT_control",mdp.NVT_control,"choose which thermostat used in NVT ensemble");

	ofs << "\n#Parameters (11.Efield)" << endl;
	OUTP(ofs,"efield",efield,"add electric field");
	OUTP(ofs,"edir",edir,"add electric field");
	OUTP(ofs,"emaxpos",emaxpos,"maximal position of efield [0,1)");
	OUTP(ofs,"eopreg",eopreg,"where sawlike potential decrease");
	OUTP(ofs,"eamp",eamp,"amplitute of the efield, unit is a.u.");
	OUTP(ofs,"eamp_v",eamp*51.44,"amplitute of the efield, unit is V/A");

	ofs << "\n#Parameters (12.Test)" << endl;
	OUTP(ofs,"out_alllog",out_alllog,"output information for each processor, when parallel");
	OUTP(ofs,"nurse", nurse,"for coders");
	OUTP(ofs,"colour", colour,"for coders, make their live colourful");
	OUTP(ofs,"t_in_h", t_in_h,"calculate the kinetic energy or not");
	OUTP(ofs,"vl_in_h", vl_in_h,"calculate the local potential or not");
	OUTP(ofs,"vnl_in_h", vnl_in_h,"calculate the nonlocal potential or not");
	OUTP(ofs,"vh_in_h", vh_in_h,"calculate the hartree potential or not");
	OUTP(ofs,"vxc_in_h", vxc_in_h,"calculate the xc potential or not");
	OUTP(ofs,"vion_in_h", vion_in_h,"calculate the local ionic potential or not");
	OUTP(ofs,"test_force", test_force, "test the force");
	OUTP(ofs,"test_stress", test_stress, "test the force");
	
	ofs << "\n#Parameters (13.Other Methods)" << endl;
	OUTP(ofs,"mlwf_flag",mlwf_flag,"turn MLWF on or off");
	OUTP(ofs,"opt_epsilon2",opt_epsilon2,"calculate the dielectic function");
	OUTP(ofs,"opt_nbands",opt_nbands,"number of bands for optical calculation");
	
	ofs << "\n#Parameters (14.VdW Correction)" << endl;								
	OUTP(ofs,"vdw_method",vdw_method,"the method of calculating vdw (none ; d2 ; d3_0 ; d3_bj");
	OUTP(ofs,"vdw_s6",vdw_s6,"scale parameter of d2/d3_0/d3_bj");
    OUTP(ofs,"vdw_s8",vdw_s8,"scale parameter of d3_0/d3_bj");
    OUTP(ofs,"vdw_a1",vdw_a1,"damping parameter of d3_0/d3_bj");
    OUTP(ofs,"vdw_a2",vdw_a2,"damping parameter of d3_bj");
    OUTP(ofs,"vdw_d",vdw_d,"damping parameter of d2");		
    OUTP(ofs,"vdw_abc",vdw_abc,"third-order term?");
	OUTP(ofs,"vdw_C6_file",vdw_C6_file,"filename of C6");
    OUTP(ofs,"vdw_C6_unit",vdw_C6_unit,"unit of C6, Jnm6/mol or eVA6");
    OUTP(ofs,"vdw_R0_file",vdw_R0_file,"filename of R0");
    OUTP(ofs,"vdw_R0_unit",vdw_R0_unit,"unit of R0, A or Bohr");
	OUTP(ofs,"vdw_model",vdw_model,"expression model of periodic structure, radius or period");
    OUTP(ofs,"vdw_radius",vdw_radius,"radius cutoff for periodic structure");
    OUTP(ofs,"vdw_radius_unit",vdw_radius_unit,"unit of radius cutoff for periodic structure");
    OUTP(ofs,"vdw_cn_thr",vdw_cn_thr,"radius cutoff for cn");
    OUTP(ofs,"vdw_cn_thr_unit",vdw_cn_thr_unit,"unit of cn_thr, Bohr or Angstrom");
	ofs << setw(20) << "vdw_period" << vdw_period.x 
			<< " " << vdw_period.y << " " 
			<< vdw_period.z<< " #periods of periodic structure" << endl;
	
	
	ofs << "\n#Parameters (15.spectrum)" << endl;              // pengfei Li add 2016-11-23
	OUTP(ofs,"spectral_type",spectral_type,"the type of the calculated spectrum");
	OUTP(ofs,"spectral_method",spectral_method,"0: tddft(linear response)");
	OUTP(ofs,"kernel_type",kernel_type,"the kernel type: rpa, tdlda ...");
	OUTP(ofs,"eels_method",eels_method,"0: hilbert_transform method; 1: standard method");
	OUTP(ofs,"absorption_method",absorption_method,"0: vasp's method  1: pwscf's method");
	OUTP(ofs,"system",system_type,"the calculate system");
	OUTP(ofs,"eta",eta,"eta(Ry)");
	OUTP(ofs,"domega",domega,"domega(Ry)");
	OUTP(ofs,"nomega",nomega,"nomega");
	OUTP(ofs,"ecut_chi",ecut_chi,"the dimension of chi matrix");
	ofs << setw(20) <<"q_start"<<q_start[0]<<"   "
		<<q_start[1]<<"   "<<q_start[2]
		<<"  #the position of the first q point in direct coordinate" <<endl;
	ofs << setw(20) <<"q_direction"<<q_direct[0]<<"   "<<q_direct[1]<<"   "<<q_direct[2]<<"  #the q direction" <<endl;
	OUTP(ofs,"nq",nq,"the total number of qpoints for calculation");
	OUTP(ofs,"out_epsilon",out_epsilon,"output epsilon or not");
	OUTP(ofs,"out_chi",out_chi,"output chi or not");
	OUTP(ofs,"out_chi0",out_chi0,"output chi0 or not");
	OUTP(ofs,"fermi_level",fermi_level,"the change of the fermi_level(Ry)");
	OUTP(ofs,"coulomb_cutoff",coulomb_cutoff," turn on the coulomb_cutoff or not");
	OUTP(ofs,"kmesh_interpolation",kmesh_interpolation,"calculting <i,0|j,R>");
	for(int i=0; i<nq; i++)
	{
		ofs << setw(20) <<"qcar" << qcar[i][0] <<"   "<< qcar[i][1] <<"   "<<qcar[i][2]<<"  #(unit: 2PI/lat0)" << endl;
	}
	OUTP(ofs,"ocp",GlobalV::ocp,"change occupation or not");
	OUTP(ofs,"ocp_set",GlobalV::ocp_set,"set occupation");
	//OUTP(ofs,"ocp_n",ocp_n,"number of occupation");
	// for(int i=0; i<ocp_n; i++)
	// {
		// ofs << setw(20) <<"ocp_kb" << GlobalV::ocp_kb[i]<< endl;
	// }
	ofs << setw(20) <<"lcao_box"<<lcao_box[0]
		<<"   "<<lcao_box[1]<<"   "
		<<lcao_box[2]<<"  #the scale for searching the existence of the overlap <i,0|j,R>" <<endl;
	OUTP(ofs," mulliken", GlobalV::mulliken," mulliken  charge or not");//qifeng add 2019/9/10
	
	//OUTP(ofs,"epsilon0",epsilon0,"calculate the macroscopic dielectric constant or not");
	OUTP(ofs,"intrasmear",intrasmear,"Eta");
	OUTP(ofs,"shift",shift,"shift");
	OUTP(ofs,"metalcalc",metalcalc,"metal or not");
	OUTP(ofs,"eps_degauss",eps_degauss,"degauss in calculating epsilon0");
	OUTP(ofs,"noncolin",noncolin,"using non-collinear-spin");
	OUTP(ofs,"lspinorb",lspinorb,"consider the spin-orbit interaction");
	
	//OUTP(ofs,"epsilon0_choice",epsilon0_choice,"0: vasp's method  1: pwscf's method");
	
	ofs << "\n#Parameters (17.exx)" << endl;
	OUTP(ofs,"exx_hybrid_type",exx_hybrid_type,"no, hf, pbe0, hse or opt_orb");
	OUTP(ofs,"exx_hybrid_alpha",exx_hybrid_alpha,"");
	OUTP(ofs,"exx_hse_omega",exx_hse_omega,"");
	OUTP(ofs,"exx_separate_loop",exx_separate_loop,"0 or 1");
	OUTP(ofs,"exx_hybrid_step",exx_hybrid_step,"");
	OUTP(ofs,"exx_lambda",exx_lambda,"");
	OUTP(ofs,"exx_pca_threshold",exx_pca_threshold,"");
	OUTP(ofs,"exx_c_threshold",exx_c_threshold,"");
	OUTP(ofs,"exx_v_threshold",exx_v_threshold,"");
	OUTP(ofs,"exx_dm_threshold",exx_dm_threshold,"");
	OUTP(ofs,"exx_schwarz_threshold",exx_schwarz_threshold,"");
	OUTP(ofs,"exx_cauchy_threshold",exx_cauchy_threshold,"");
	OUTP(ofs,"exx_ccp_threshold",exx_ccp_threshold,"");
	OUTP(ofs,"exx_ccp_rmesh_times",exx_ccp_rmesh_times,"");
	OUTP(ofs,"exx_distribute_type",exx_distribute_type,"htime or kmeans1 or kmeans2");
	OUTP(ofs,"exx_opt_orb_lmax",exx_opt_orb_lmax,"");
	OUTP(ofs,"exx_opt_orb_ecut",exx_opt_orb_ecut,"");
	OUTP(ofs,"exx_opt_orb_tolerence",exx_opt_orb_tolerence,"");

	ofs << "\n#Parameters (17.tddft)" << endl;
	OUTP(ofs,"tddft",tddft,"calculate tddft or not");
	OUTP(ofs,"td_dr2",td_dr2,"threshold for electronic iteration of tddft");
	OUTP(ofs,"td_dt",td_dt,"time of ion step");
	OUTP(ofs,"td_force_dt",td_force_dt,"time of force change");
	OUTP(ofs,"td_val_elec_01",td_val_elec_01,"td_val_elec_01");
	OUTP(ofs,"td_val_elec_02",td_val_elec_02,"td_val_elec_02");
	OUTP(ofs,"td_val_elec_03",td_val_elec_03,"td_val_elec_03");
	OUTP(ofs,"td_vext",td_vext,"add extern potential or not");
	OUTP(ofs,"td_vext_dire",td_vext_dire,"extern potential direction");
	OUTP(ofs,"td_timescale",td_timescale,"extern potential td_timescale");
	OUTP(ofs,"td_vexttype",td_vexttype,"extern potential type");
	OUTP(ofs,"td_vextout",td_vextout,"output extern potential or not");
	OUTP(ofs,"td_dipoleout",td_dipoleout,"output dipole or not");

	
	ofs << "\n#Parameters (18.berry_wannier)" << endl;
	OUTP(ofs,"berry_phase",berry_phase,"calculate berry phase or not");
	OUTP(ofs,"gdir",gdir,"calculate the polarization in the direction of the lattice vector");
	OUTP(ofs,"towannier90",towannier90,"use wannier90 code interface or not");
	OUTP(ofs,"nnkpfile",NNKP,"the wannier90 code nnkp file name");
	OUTP(ofs,"wannier_spin",wannier_spin,"calculate spin in wannier90 code interface");
	
    ofs.close();
    return;
}
