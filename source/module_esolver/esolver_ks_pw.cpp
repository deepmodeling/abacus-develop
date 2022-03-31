#include "esolver_ks_pw.h"

//--------------temporary----------------------------
#include "../src_pw/global.h"
#include "../module_base/global_function.h"
#include "../module_symmetry/symmetry.h"
#include "../src_pw/vdwd2.h"
#include "../src_pw/vdwd3.h"
#include "../src_pw/vdwd2_parameters.h"
#include "../src_pw/vdwd3_parameters.h"
#include "../src_pw/pw_complement.h"
#include "../src_pw/pw_basis.h"
#include "../src_io/print_info.h"
//-----force-------------------
#include "../src_pw/forces.h"
//-----stress------------------
#include "../src_pw/stress_pw.h"
//---------------------------------------------------

namespace ModuleESolver
{

void ESolver_KS_PW::Init(Input &inp, UnitCell_pseudo &ucell)
{
    // setup GlobalV::NBANDS 
	// Yu Liu add 2021-07-03
	GlobalC::CHR.cal_nelec();

	if(GlobalC::ucell.atoms[0].xc_func=="HSE"||GlobalC::ucell.atoms[0].xc_func=="PBE0")
	{
		XC_Functional::set_xc_type("pbe");
	}
	else
	{
		XC_Functional::set_xc_type(GlobalC::ucell.atoms[0].xc_func);
	}

    ModuleBase::GlobalFunc::DONE(GlobalV::ofs_running, "SETUP UNITCELL");

    //-----------------------------------------------------

    // symmetry analysis should be performed every time the cell is changed
    if (ModuleSymmetry::Symmetry::symm_flag)
    {
        GlobalC::symm.analy_sys(GlobalC::ucell, GlobalV::ofs_running);
        ModuleBase::GlobalFunc::DONE(GlobalV::ofs_running, "SYMMETRY");
    }

    // Setup the k points according to symmetry.
    GlobalC::kv.set( GlobalC::symm, GlobalV::global_kpoint_card, GlobalV::NSPIN, GlobalC::ucell.G, GlobalC::ucell.latvec );
    ModuleBase::GlobalFunc::DONE(GlobalV::ofs_running,"INIT K-POINTS");

    // print information
    // mohan add 2021-01-30
    Print_Info::setup_parameters(GlobalC::ucell, GlobalC::kv);

    // Initalize the plane wave basis set
    GlobalC::pw.gen_pw(GlobalV::ofs_running, GlobalC::ucell, GlobalC::kv);
    ModuleBase::GlobalFunc::DONE(GlobalV::ofs_running,"INIT PLANEWAVE");
    std::cout << " UNIFORM GRID DIM     : " << GlobalC::pw.nx <<" * " << GlobalC::pw.ny <<" * "<< GlobalC::pw.nz << std::endl;
    std::cout << " UNIFORM GRID DIM(BIG): " << GlobalC::pw.nbx <<" * " << GlobalC::pw.nby <<" * "<< GlobalC::pw.nbz << std::endl;

    // mohan add 2010-09-13
    // initialize the real-space uniform grid for FFT and parallel
    // distribution of plane waves
    GlobalC::Pgrid.init(GlobalC::pw.ncx, GlobalC::pw.ncy, GlobalC::pw.ncz, GlobalC::pw.nczp,
    GlobalC::pw.nrxx, GlobalC::pw.nbz, GlobalC::pw.bz); // mohan add 2010-07-22, update 2011-05-04
        

    // Calculate Structure factor
    GlobalC::pw.setup_structure_factor();
    // cout<<"after pgrid init nrxx = "<<GlobalC::pw.nrxx<<endl;
    
    //----------------------------------------------------------
    // 1 read in initial data:
    //   a lattice structure:atom_species,atom_positions,lattice vector
    //   b k_points
    //   c pseudopotential
    // 2 setup planeware basis, FFT,structure factor, ...
    // 3 initialize local and nonlocal pseudopotential in G_space
    // 4 initialize charge desity and warefunctios in G_space
    //----------------------------------------------------------

    //=====================================
    // init charge/potential/wave functions
    //=====================================
    GlobalC::CHR.allocate(GlobalV::NSPIN, GlobalC::pw.nrxx, GlobalC::pw.ngmc);
    GlobalC::pot.allocate(GlobalC::pw.nrxx);

    GlobalC::wf.allocate(GlobalC::kv.nks);

    // cout<<GlobalC::pw.nrxx<<endl;
    // cout<<"before ufft allocate"<<endl;
    GlobalC::UFFT.allocate();

    // cout<<"after ufft allocate"<<endl;

    //=======================
    // init pseudopotential
    //=======================
    GlobalC::ppcell.init(GlobalC::ucell.ntype);

    //=====================
    // init hamiltonian
    // only allocate in the beginning of ELEC LOOP!
    //=====================
    GlobalC::hm.hpw.allocate(GlobalC::wf.npwx, GlobalV::NPOL, GlobalC::ppcell.nkb, GlobalC::pw.nrxx);

    //=================================
    // initalize local pseudopotential
    //=================================
    GlobalC::ppcell.init_vloc(GlobalC::pw.nggm, GlobalC::ppcell.vloc);
    ModuleBase::GlobalFunc::DONE(GlobalV::ofs_running, "LOCAL POTENTIAL");

    //======================================
    // Initalize non local pseudopotential
    //======================================
    GlobalC::ppcell.init_vnl(GlobalC::ucell);
    ModuleBase::GlobalFunc::DONE(GlobalV::ofs_running, "NON-LOCAL POTENTIAL");

    //=========================================================
    // calculate the total local pseudopotential in real space
    //=========================================================
    GlobalC::pot.init_pot(0, GlobalC::pw.strucFac); //atomic_rho, v_of_rho, set_vrs

    GlobalC::pot.newd();

    ModuleBase::GlobalFunc::DONE(GlobalV::ofs_running, "INIT POTENTIAL");

    //==================================================
    // create GlobalC::ppcell.tab_at , for trial wave functions.
    //==================================================
    GlobalC::wf.init_at_1();

    //================================
    // Initial start wave functions
    //================================
    if (GlobalV::NBANDS != 0 || (GlobalV::CALCULATION != "scf-sto" && GlobalV::CALCULATION != "relax-sto" && GlobalV::CALCULATION != "md-sto")) //qianrui add
    {
        GlobalC::wf.wfcinit();
    }

#ifdef __LCAO
#ifdef __MPI
    switch (GlobalC::exx_global.info.hybrid_type) // Peize Lin add 2019-03-09
    {
    case Exx_Global::Hybrid_Type::HF:
    case Exx_Global::Hybrid_Type::PBE0:
    case Exx_Global::Hybrid_Type::HSE:
        GlobalC::exx_lip.init(&GlobalC::kv, &GlobalC::wf, &GlobalC::pw, &GlobalC::UFFT, &GlobalC::ucell);
        break;
    case Exx_Global::Hybrid_Type::No:
        break;
    case Exx_Global::Hybrid_Type::Generate_Matrix:
    default:
        throw std::invalid_argument(ModuleBase::GlobalFunc::TO_STRING(__FILE__) + ModuleBase::GlobalFunc::TO_STRING(__LINE__));
    }
#endif
#endif

    ModuleBase::GlobalFunc::DONE(GlobalV::ofs_running, "INIT BASIS");



}

void ESolver_KS_PW::Run(int istep,UnitCell_pseudo &ucell)
{
    elec.self_consistent(istep-1);
    return ;
}

void ESolver_KS_PW::cal_Energy(energy &en)
{

}

void ESolver_KS_PW::cal_Force(ModuleBase::matrix &force)
{
	Forces ff;
	ff.init(force);
}
void ESolver_KS_PW::cal_Stress(ModuleBase::matrix &stress)
{
	Stress_PW ss;
	ss.cal_stress(stress);
}

int ESolver_KS_PW:: getiter()
{
    return this->elec.iter;
}

}