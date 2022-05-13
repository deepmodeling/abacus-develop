#include "gtest/gtest.h"
#include "mpi.h"
#include <iostream>
#include <fstream>
#include <complex>

#include "gamma_rho_mock.h"

/***************************************************************
*  unit test of functions in src_lcao/gint_gamma_rho.cpp
****************************************************************/

/**
 * This unit test is designed to test the functions in src_lcao/
 * gint_gamma_rho.cpp: Gint_Gamma::cal_rho(), Gint_Gamma::gamma_charge()
 * and sum_up_rho. 
 * It can make two comparisons:
 * (1) compare density matrix calculated from wavefunction with that read
 * from file (LOWF_GAMMA_S1.dat).
 * (2) compare charge density calcualted from density matrix with that read
 * from file (SPIN1_DM_Gamma)
 * The input and output files were taken from an integrated test:
 * abacus-develop/tests/integrate/815_NO_LT_sc by adding the following
 * INPUT tags:
 *       chg_extrap dm
 *       init_wfc atomic
 *       gamma_only 1
 *       out_wfc_lcao 1
 *       out_chg 1
 *       out_dm 1
 */

namespace
{

/******************************
 * Prepare Input variables
 ******************************/
struct ENVPrepare
{
    std::string basis_type_;
    int         nelec_;
    std::string calculation_;
    std::string pseudo_dir_;
    std::string stru_file_;
    std::string running_log_;
    std::string latname_;
    int         ntype_;
    int         lmaxmax_;
    bool        init_vel_;
    std::string fixed_axes_;
    double      pseudo_rcut_;
    bool        symm_flag_;
    std::string kpoint_card_;
    int         nspin_;
    bool        gamma_only_;
    double      ecutwfc_;
    double      ecutrho_;  
    int         nx_;
    int         ny_;
    int         nz_;
    int         ncx_;
    int         ncy_;
    int         ncz_;
    int         bx_;
    int         by_;
    int         bz_;
    int         pw_seed_;
    int         nbspline_;
    int         test_pw_;
    int         nbands_;
    double      cell_factor_;
    std::string init_chg_;
    std::string init_wfc_;
    int         out_wfc_pw_;
    std::string ks_solver_;
    std::string pseudo_type_; // used in init_pseudo_reader(const std::string &fn)
    int         nb2d_;
    int         colour_;
    int         drank_;
    int         my_rank_;
    std::string orbital_dir_;
    double lcao_ecut_; // ecut of two center integral
    double lcao_dk_; // delta k used in two center integral
    double lcao_dr_; // dr used in two center integral
    double lcao_rmax_; // rmax(a.u.) to make table.
    std::string out_dir_;

    // default values
    ENVPrepare()
    {
	basis_type_     = "lcao";
	nelec_          = 2;
	calculation_    = "scf";
    	pseudo_dir_	= "./support/";		
    	stru_file_	= "./support/STRU";		
    	running_log_	= "./support/running.log";	
    	latname_	= "sc";		
    	ntype_		= 1;		
    	lmaxmax_	= 2;		
    	init_vel_	= false;		
    	fixed_axes_	= "None";		
    	pseudo_rcut_	= 15.0;		
    	symm_flag_	= false;		
    	kpoint_card_	= "./support/KPT_Gamma";		
    	nspin_		= 1;		
    	gamma_only_	= true;          
    	ecutwfc_	= 25.0;           
    	ecutrho_  	= 0.0;            
    	nx_		= 0;              
    	ny_		= 0;              
    	nz_		= 0;              
    	ncx_		= 0;              
    	ncy_		= 0;              
    	ncz_		= 0;              
    	bx_		= 2;              
    	by_		= 2;              
    	bz_		= 2;              
    	pw_seed_	= 1;              
    	nbspline_	= -1;             
    	test_pw_	= 2;		
    	nbands_		= 4;
	cell_factor_    = 1.2;
	init_chg_       = "atomic";
	init_wfc_       = "atomic";
	out_wfc_pw_     = 1;
	ks_solver_      = "genelpa"; // default in input.cpp
	pseudo_type_    = "auto"; // default in input.cpp, not in global_variable.cpp
	nb2d_           = 0; // default in input.cpp
	colour_         = 0; // default in input.cpp
	drank_          = 0;
	my_rank_        = 0;
	orbital_dir_    = "./support/";
	//needed in read_orb_file
        lcao_ecut_ = ecutwfc_; // (Ry)
        lcao_dk_ = 0.01;
        lcao_dr_ = 0.01;
        lcao_rmax_ = 30; // (a.u.)
	out_dir_   = "./support/";
    }
};

ENVPrepare ENVP;
std::string running_log;	// to replace GlobalV::ofs_running
Local_Orbital_Charge LOC;
Local_Orbital_wfc LOWF;
Parallel_Orbitals ParaV;
Gint_Gamma GG;


class ENVEnvironment : public ::testing::Environment
{
public:
  ~ENVEnvironment(){}

  // Here is the interface to ABACUS.
  // We need to modify here after reconstruction.
  void SetUp() override 
  {
	  GlobalV::BASIS_TYPE        = env->basis_type_;
	  GlobalC::CHR.nelec         = env->nelec_;
	  GlobalV::CALCULATION       = env->calculation_;
	  GlobalV::global_pseudo_dir = env->pseudo_dir_;
	  GlobalV::stru_file         = env->stru_file_;
	  running_log                = env->running_log_;
	  // INPUT is an object of Input, and declared in input.h
	  GlobalC::ucell.latName     = env->latname_;
	  GlobalC::ucell.ntype       = env->ntype_;
	  INPUT.lmaxmax              = env->lmaxmax_;
	  INPUT.init_vel              = env->init_vel_;
	  INPUT.fixed_axes           = env->fixed_axes_;
	  // important in pseudo_nc::set_pseudo_atom
	  GlobalV::PSEUDORCUT        = env->pseudo_rcut_;
	  ModuleSymmetry::Symmetry::symm_flag 
                                     = env->symm_flag_;
	  GlobalV::NSPIN             = env->nspin_;
	  GlobalV::global_kpoint_card 
                                     = env->kpoint_card_;
	  GlobalV::GAMMA_ONLY_LOCAL  = env->gamma_only_; // set in input.cpp
          INPUT.ecutwfc              = env->ecutwfc_;
          INPUT.ecutrho              = env->ecutrho_;    
          INPUT.nx                   = env->nx_;
          INPUT.ny                   = env->ny_;
          INPUT.nz                   = env->nz_;
          INPUT.ncx                  = env->ncx_;
          INPUT.ncy                  = env->ncy_;
          INPUT.ncz                  = env->ncz_;
          INPUT.bx                   = env->bx_;
          INPUT.by                   = env->by_;
          INPUT.bz                   = env->bz_;
          INPUT.pw_seed              = env->pw_seed_;
          INPUT.nbspline             = env->nbspline_;
	  GlobalV::test_pw           = env->test_pw_;
	  GlobalV::NBANDS            = env->nbands_;
	  INPUT.cell_factor          = env->cell_factor_;
	  GlobalC::pot.init_chg      = env->init_chg_;
	  GlobalC::wf.init_wfc       = env->init_wfc_;
	  GlobalC::wf.out_wfc_pw     = env->out_wfc_pw_;
	  GlobalV::KS_SOLVER         = env->ks_solver_;
	  GlobalV::global_pseudo_type= env->pseudo_type_;
	  GlobalV::NB2D              = env->nb2d_;
	  GlobalV::DCOLOR            = env->colour_;
	  GlobalV::DRANK             = env->drank_;
	  GlobalV::MY_RANK           = env->my_rank_;
	  GlobalV::global_orbital_dir= env->orbital_dir_;
	  INPUT.lcao_ecut            = env->lcao_ecut_;
	  INPUT.lcao_dk              = env->lcao_dk_;
	  INPUT.lcao_dr              = env->lcao_dr_;
	  INPUT.lcao_rmax            = env->lcao_rmax_;
          GlobalV::global_out_dir    = env->out_dir_;
  }

  void set_variables_in_set_up(ENVPrepare* envtmp)
  {
	  env = envtmp;
  }
private:
  ENVPrepare* env;
};

class LCAOTest : public testing::Test
{
protected:
	void SetUp()
	{
		// here we use GlobalV::ofs_running directly
		GlobalV::ofs_running.open(running_log.c_str());
	}
	void TearDown()
	{
		GlobalV::ofs_running.close();
	}
};

TEST_F(LCAOTest,GammaRho)
{
	// initialization
	ModuleESolver::ESolver *p_esolver;
	Run_lcao::lcao_line(p_esolver);
	// initialization happened in ORB_control::setup_2d_division()
	// and ORB_control::set_trace()
	ParaV.nloc = GlobalV::NLOCAL*GlobalV::NLOCAL;
	ParaV.trace_loc_row = new int[GlobalV::NLOCAL];
	ParaV.trace_loc_col = new int[GlobalV::NLOCAL];
	for (int i=0;i<GlobalV::NLOCAL;i++)
	{
		ParaV.trace_loc_row[i] = i;
		ParaV.trace_loc_col[i] = i;
	}
	ParaV.ncol = GlobalV::NLOCAL;
	ParaV.nrow = GlobalV::NLOCAL;
	LOC.ParaV = &ParaV;
	//
        GlobalC::wf.allocate_ekb_wg(GlobalC::kv.nks);
        GlobalC::UFFT.allocate();
	// Grid Technique
	set_matrix_grid();

	// read gamma wavefunction
	LOWF.wfc_gamma.resize(GlobalV::NSPIN);
	LOC.gamma_file(LOWF);
	// allocate space for DM
	GlobalC::GridT.lgd = GlobalV::NLOCAL;
	LOC.allocate_gamma(GlobalC::GridT.lgd);
	// dm_gamma is another way to save density matrix
	LOC.dm_gamma.resize(GlobalV::NSPIN); // originally inside allocate_gamma
	// calculate density matrix from wavefunction
	LOC.cal_dm(GlobalC::wf.wg,LOWF.wfc_gamma,LOC.dm_gamma);
	//LOC.cal_dk_gamma_from_2D_pub(); ignore this step because it requres mpi
	// from dm_gamma to DM
	// create a temprary quantity to make comparison.
        std::vector<ModuleBase::matrix> DM_calculated;	
	DM_calculated.resize(GlobalV::NSPIN);
        for(int is=0; is<GlobalV::NSPIN; ++is)
        {
		DM_calculated[is].create(ParaV.ncol,ParaV.nrow);
	}
        for(int is=0; is<GlobalV::NSPIN; ++is)
        {
            for (int i = 0; i < GlobalV::NLOCAL; ++i)
            {
                for (int j = 0; j < GlobalV::NLOCAL; ++j)
                {
                     DM_calculated[is](i,j)=LOC.dm_gamma[is](j,i);
		     //std::cout<<DM_calculated[is](i,j)<<" ";
                }
		//std::cout<<endl;
            }
	}
	// read density matrix
	for(int is=0; is<GlobalV::NSPIN; is++)
	{
		ModuleBase::GlobalFunc::ZEROS(GlobalC::CHR.rho[is], GlobalC::pw.nrxx);
		std::stringstream ssd;
		ssd << GlobalV::global_out_dir << "SPIN" << is + 1 << "_DM_Gamma" ;
		//std::cout<<"ssd "<<ssd.str()<<std::endl;
		// reading density matrix,
		LOC.read_dm(is, ssd.str() );
	}
	// compare density matrix read and calculated
        for(int is=0; is<GlobalV::NSPIN; ++is)
        {
            for (int i = 0; i < GlobalV::NLOCAL; ++i)
            {
                for (int j = 0; j < GlobalV::NLOCAL; ++j)
                {
		     EXPECT_NEAR(DM_calculated[is](i,j),LOC.DM[is][i][j],1e-4);
		     //std::cout<<LOC.DM[is][i][j]<<" ";
                }
		//std::cout<<endl;
            }
	}
	// calculate the charge density
	if(GlobalV::GAMMA_ONLY_LOCAL)
	{
		GG.cal_rho(LOC.DM, (Charge*)(&GlobalC::CHR));
		//std::cout<<"number of elec: "<<nelec<<std::endl;
	}
	//std::cout<<"rho in test "<<GlobalC::CHR.rho[0][0]<<std::endl;

	GlobalC::CHR.renormalize_rho();
	//====== read rho ===================================================
	double** rho_for_compare;
	rho_for_compare = new double*[GlobalV::NSPIN];
	double totale = 0.0;
	for (int is = 0; is < GlobalV::NSPIN; is++)
	{
	    rho_for_compare[is] = new double[GlobalC::pw.nrxx];
	    std::stringstream ssc;
	    ssc << GlobalV::global_out_dir<< "SPIN" << is + 1 << "_CHG_Gamma";
	    GlobalC::CHR.read_rho(is, ssc.str(), rho_for_compare[is]);
	    for (int ix = 0; ix < GlobalC::pw.nrxx; ix++)
	    //for (int ix = 0; ix < 5; ix++)
	    {
	        totale += rho_for_compare[is][ix];
	        // compare rho read and rho calculated from wavefunctions
	        //std::cout<<"read "<< rho_for_compare[is][ix]<<" calc "<<GlobalC::CHR.rho[is][ix]<<std::endl;
	        EXPECT_NEAR(rho_for_compare[is][ix], GlobalC::CHR.rho[is][ix], 1e-5);
	    }
	}
	// check total number of electrons
	totale = totale * GlobalC::ucell.omega / GlobalC::pw.nrxx;
	EXPECT_NEAR(totale, GlobalC::CHR.nelec, 1e-8);
}


void Check(bool condition, const char* msg)
{
  if (!condition)
  {
    printf("FAILED: %s\n", msg);
  }
}

int RunAllTests(ENVEnvironment* env, ENVPrepare* ENVP)
{
    //env->Reset();
    env->set_variables_in_set_up(ENVP);
    return RUN_ALL_TESTS();
}

}//namespace

int main(int argc, char **argv)
{

    MPI_Init(&argc, &argv);

    testing::InitGoogleTest(&argc, argv);

    ENVEnvironment* const env = new ENVEnvironment;
    testing::AddGlobalTestEnvironment(env);
    Check (RunAllTests(env,&ENVP)==0,"");

    MPI_Finalize();

    return 0;
}
