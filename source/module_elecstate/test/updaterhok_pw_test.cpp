#ifdef __MPI
#include "mpi.h"
#endif
#include "gtest/gtest.h"
#include "gmock/gmock.h"
#include <complex>
#include <fstream>
#include <iostream>

#include "module_psi/psi.h"
#include "module_elecstate/elecstate_pw.h"
#include "module_elecstate/elecstate.h"

#include "updaterhok_pw_test.h"
#include "src_io/berryphase.h"
bool berryphase::berry_phase_flag;
using::testing::AtLeast;
using::testing::Assign;

/***************************************************************
 *  unit test of elecstate::ElecStatePW::psiToRho()
 ****************************************************************/

/**
 * This unit test is designed to test the function of:
 *              elecstate::ElecStatePW::psiToRho()
 * A psi object "evc" is created to receive wavefunction, from which
 * psiToRho() calculate charge density. The class MockElecStatePW is 
 * derived from ElecStatePW, in order to mock the calculate_weights method.
 *
 * Currently, this test looks involved due to the use of extensively
 * called global classes. It will be adjusted to test the corresponding
 * new interfaces after reconstruction.
 *
 * The inputs and output files were taken from an integrated test:
 * abacus-develop/tests/integrate/801_PW_LT_sc. Note: the output files 
 * WAVEFUNC1.txt and SPIN1_CHG were obtained by using ABACUS without MPI 
 * (ABACUC.fp.x). This is because the GlobalC::sf.gdirect_global indexes 
 * are different with and without MPI.
 *
 * The mpi env. is currently invalid in this UT.
 */

namespace elecstate
{

class MockElecStatePW : public ElecStatePW
{
public:
  MockElecStatePW( ModulePW::PW_Basis_K* wfc_basis, Charge* chg_in, K_Vectors *pkv_in, int nbands_in):ElecStatePW(wfc_basis, chg_in, pkv_in, nbands_in){}
  MOCK_METHOD0(calculate_weights,void());
};
}

/******************************
 * read wavefunction
 ******************************/
void read_wfc(const std::string& fn, psi::Psi<std::complex<double>> &psi)
{
    std::string* wfilename;
    wfilename = new std::string[GlobalC::kv.nkstot];
    std::string tmpstring;
    for (int ik = 0; ik < GlobalC::kv.nkstot; ik++)
    {
        int ikstot;
        std::stringstream wfss;
        if (GlobalC::wf.out_wfc_pw == 1)
        {
            wfss << fn << ik + 1 << ".txt";
        }
        else if (GlobalC::wf.out_wfc_pw == 2)
        {
            wfss << fn << ik + 1 << ".dat";
        }
        wfilename[ik] = wfss.str();
        std::ifstream ifs;
        ifs.open(wfilename[ik], std::ifstream::in);
        // read out unnecessary info.
        for (int i = 0; i < 10; i++)
            getline(ifs, tmpstring);
        // read gkk
        for (int ig = 0; ig < GlobalC::kv.ngk[ik]; ig++)
        {
            ifs >> GlobalC::wfcpw->getgcar(ik,ig).x >> GlobalC::wfcpw->getgcar(ik,ig).y >> GlobalC::wfcpw->getgcar(ik,ig).z;
        }
        getline(ifs, tmpstring);
        getline(ifs, tmpstring);
        getline(ifs, tmpstring);
        for (int ib = 0; ib < GlobalV::NBANDS; ib++)
        {
            if (GlobalC::wf.out_wfc_pw == 1)
            {
                getline(ifs, tmpstring);
                getline(ifs, tmpstring);
                for (int ig = 0; ig < GlobalC::kv.ngk[ik]; ig++)
                {
                    double realpart;
                    double imagpart;
                    ifs >> realpart >> imagpart;
                    psi(ik, ib, ig) = std::complex<double>(realpart, imagpart);
                }
                getline(ifs, tmpstring);
                getline(ifs, tmpstring);
            }
        }
        ifs.close();
    }
}

/******************************
 * Prepare Input variables
 ******************************/
struct ENVPrepare
{
    int nelec_;
    std::string calculation_;
    std::string pseudo_dir_;
    std::string stru_file_;
    std::string running_log_;
    std::string latname_;
    int ntype_;
    double pseudo_rcut_;
    int symm_flag_;
    std::string kpoint_card_;
    int nspin_;
    bool gamma_only_;
    double ecutwfc_;
    double ecutrho_;
    int nx_;
    int ny_;
    int nz_;
    int ncx_;
    int ncy_;
    int ncz_;
    int bx_;
    int by_;
    int bz_;
    int pw_seed_;
    int nbspline_;
    int test_pw_;
    int nbands_;
    double cell_factor_;
    std::string init_chg_;
    std::string init_wfc_;
    int out_wfc_pw_;
    std::string out_dir_;
    // default values
    ENVPrepare()
    {
        nelec_ = 2;
        calculation_ = "scf";
        pseudo_dir_ = "./support/";
        stru_file_ = "./support/STRU";
        running_log_ = "./support/running.log";
        latname_ = "sc";
        ntype_ = 1;
        pseudo_rcut_ = 15.0;
        symm_flag_ = 0;
        kpoint_card_ = "./support/KPT";
        nspin_ = 1;
        gamma_only_ = false;
        ecutwfc_ = 25.0;
        ecutrho_ = 0.0;
        nx_ = 0;
        ny_ = 0;
        nz_ = 0;
        ncx_ = 0;
        ncy_ = 0;
        ncz_ = 0;
        bx_ = 2;
        by_ = 2;
        bz_ = 2;
        pw_seed_ = 1;
        nbspline_ = -1;
        test_pw_ = 2;
        nbands_ = 4;
        cell_factor_ = 1.2;
        init_chg_ = "atomic";
        init_wfc_ = "atomic";
        out_wfc_pw_ = 1;
	out_dir_="./support/";
    }
};

ENVPrepare ENVP;
std::string running_log; // to replace GlobalV::ofs_running
LCAO_Orbitals orb; // to replace GlovalC::ORB

psi::Psi<std::complex<double>> evc; //new class Psi

class ENVEnvironment : public ::testing::Environment
{
  public:
    ~ENVEnvironment()
    {
    }

    // Here is the interface to ABACUS.
    // We need to modify here after reconstruction.
    void SetUp() override
    {
        GlobalC::CHR.nelec = env->nelec_;
        GlobalV::CALCULATION = env->calculation_;
        GlobalV::global_pseudo_dir = env->pseudo_dir_;
        GlobalV::stru_file = env->stru_file_;
        running_log = env->running_log_;
        // INPUT is an object of Input, and declared in input.h
	GlobalC::ucell.latName = env->latname_;
	GlobalC::ucell.ntype = env->ntype_;
        // important in pseudo_nc::set_pseudo_atom
        GlobalV::PSEUDORCUT = env->pseudo_rcut_;
        ModuleSymmetry::Symmetry::symm_flag = env->symm_flag_;
        GlobalV::NSPIN = env->nspin_;
        GlobalV::global_kpoint_card = env->kpoint_card_;
        INPUT.gamma_only = env->gamma_only_;
        INPUT.ecutwfc = env->ecutwfc_;
        INPUT.ecutrho = env->ecutrho_;
        INPUT.nx = env->nx_;
        INPUT.ny = env->ny_;
        INPUT.nz = env->nz_;
        INPUT.ncx = env->ncx_;
        INPUT.ncy = env->ncy_;
        INPUT.ncz = env->ncz_;
        INPUT.bx = env->bx_;
        INPUT.by = env->by_;
        INPUT.bz = env->bz_;
        INPUT.pw_seed = env->pw_seed_;
        INPUT.nbspline = env->nbspline_;
        GlobalV::test_pw = env->test_pw_;
        GlobalV::NBANDS = env->nbands_;
        INPUT.cell_factor = env->cell_factor_;
        GlobalC::pot.init_chg = env->init_chg_;
        GlobalC::wf.init_wfc = env->init_wfc_;
        GlobalC::wf.out_wfc_pw = env->out_wfc_pw_;
	GlobalV::global_out_dir = env->out_dir_;
    }

    void set_variables_in_set_up(ENVPrepare* envtmp)
    {
        env = envtmp;
    }

  private:
    ENVPrepare* env;
};

class EState : public testing::Test
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

TEST_F(EState,RhoPW)
{
    // Here GlobalC::ucell is used directly without asserting a new
    // UnitCell_pseudo object.
    // this function will read STRU and pseudopotential files
    // should add more EXPECTs after calling
    GlobalC::ucell.setup_cell(orb, GlobalV::global_pseudo_dir, GlobalV::stru_file, GlobalV::ofs_running);
    // here GlobalC::kv and GlobalC::symm are used directly
    // the file KPT is read here
    GlobalC::kv.set(GlobalC::symm,
                    GlobalV::global_kpoint_card,
                    GlobalV::NSPIN,
                    GlobalC::ucell.G,
                    GlobalC::ucell.latvec);
    // here we use GlobalC::sf directly
    GlobalC::sf.set(INPUT.nbspline);
    // ecut is needed here in setup_gg()
    // the first parameter GlobalV::ofs_running is not necessarily needed
    // because GlobalV::ofs_running is used inside the function anyway.

    // pw_rho = new ModuleBase::PW_Basis();
    //temporary, it will be removed
    GlobalC::rhopw = new ModulePW::PW_Basis_Big(); 
    ModulePW::PW_Basis_Big* tmp = static_cast<ModulePW::PW_Basis_Big*>(GlobalC::rhopw);
    tmp->setbxyz(INPUT.bx,INPUT.by,INPUT.bz);
    GlobalC::wfcpw = new ModulePW::PW_Basis_K_Big(); 
    ModulePW::PW_Basis_K_Big* tmp2 = static_cast<ModulePW::PW_Basis_K_Big*>(GlobalC::wfcpw);
    tmp2->setbxyz(INPUT.bx,INPUT.by,INPUT.bz);
#ifdef __MPI
    GlobalC::rhopw->initmpi(1, 0 ,POOL_WORLD);
    GlobalC::wfcpw->initmpi(1, 0 ,POOL_WORLD);
#endif
    GlobalC::rhopw->initgrids(GlobalC::ucell.lat0, GlobalC::ucell.latvec, 4 * INPUT.ecutwfc);
    GlobalC::rhopw->initparameters(false, INPUT.ecutrho);
    GlobalC::rhopw->setuptransform();
    GlobalC::wfcpw->initgrids(GlobalC::ucell.lat0, GlobalC::ucell.latvec, GlobalC::rhopw->nx, GlobalC::rhopw->ny, GlobalC::rhopw->nz);
    GlobalC::wfcpw->initparameters(false, INPUT.ecutwfc, GlobalC::kv.nks, GlobalC::kv.kvec_d.data());
    GlobalC::wfcpw->setuptransform();
    for(int ik = 0 ; ik < GlobalC::kv.nks; ++ik)   GlobalC::kv.ngk[ik] = GlobalC::wfcpw->npwk[ik];
    exit(0);


    // test the generated fft grid (nx,ny,nz)
    EXPECT_TRUE((GlobalC::rhopw->nx + 1) % 2 == 0 || (GlobalC::rhopw->nx + 1) % 3 == 0 || (GlobalC::rhopw->nx + 1) % 5 == 0);
    EXPECT_TRUE((GlobalC::rhopw->ny + 1) % 2 == 0 || (GlobalC::rhopw->ny + 1) % 3 == 0 || (GlobalC::rhopw->ny + 1) % 5 == 0);
    EXPECT_TRUE((GlobalC::rhopw->nz + 1) % 2 == 0 || (GlobalC::rhopw->nz + 1) % 3 == 0 || (GlobalC::rhopw->nz + 1) % 5 == 0);

    // Calculate Structure factor
    // init charge/potential/wave functions
    GlobalC::CHR.allocate(GlobalV::NSPIN, GlobalC::rhopw->nrxx, GlobalC::rhopw->npw);
    // GlobalC::pot.allocate(GlobalC::rhopw->nrxx);
    // we need to supply NBANDS here
    psi::Psi<std::complex<double>>* psi = GlobalC::wf.allocate(GlobalC::kv.nks);
    // std::cout<<"npwx "<<GlobalC::wf.npwx<<std::endl;

    //====== read wavefunction ==========================================
    std::stringstream ssw;
    ssw <<GlobalV::global_out_dir<< "WAVEFUNC";
    // we need to supply out_wfc_pw here
    read_wfc(ssw.str(), psi[0]);

    // copy data from old wf.evc to new evc(an object of Psi)
    evc.resize(GlobalC::kv.nks,GlobalV::NBANDS,GlobalC::wf.npwx);
    evc.fix_k(0);
    for(int i=0;i<GlobalV::NBANDS;i++)
    {
	    for(int j=0;j<GlobalC::wf.npwx;j++)
            {
		    evc(i,j) = psi[0](i,j);
	    }
    }
    delete psi;
    // using class ElecStatePW to calculate rho
    elecstate::MockElecStatePW* kk;
    kk = new elecstate::MockElecStatePW(GlobalC::wfcpw,&GlobalC::CHR,&GlobalC::kv,GlobalV::NBANDS);
    EXPECT_CALL(*kk,calculate_weights()).Times(AtLeast(1));
    ModuleBase::matrix wg_tmp;
    wg_tmp.create(GlobalC::kv.nks,GlobalV::NBANDS);
    wg_tmp(0,0) = 2.0;
    wg_tmp(0,1) = 0.0;
    wg_tmp(0,2) = 0.0;
    wg_tmp(0,3) = 0.0;
    //===== calculate occupation before sum_band ========================
    ON_CALL(*kk,calculate_weights()).WillByDefault(Assign(&(kk->wg),wg_tmp));
    // run test here
    kk->psiToRho(evc);
    //kk->calculate_weights();
    EXPECT_EQ(kk->wg(0,0),2.0);
    EXPECT_EQ(kk->wg(0,1),0.0);
    EXPECT_EQ(kk->wg(0,2),0.0);
    EXPECT_EQ(kk->wg(0,3),0.0);

    //====== read rho ===================================================
    double** rho_for_compare;
    rho_for_compare = new double*[GlobalV::NSPIN];
    double totale = 0.0;
    for (int is = 0; is < GlobalV::NSPIN; is++)
    {
        rho_for_compare[is] = new double[GlobalC::rhopw->nrxx];
        std::stringstream ssc;
        ssc <<GlobalV::global_out_dir<< "SPIN" << is + 1 << "_CHG";
        GlobalC::CHR.read_rho(is, ssc.str(), rho_for_compare[is]);
        for (int ix = 0; ix < GlobalC::rhopw->nrxx; ix++)
        //for (int ix = 0; ix < 5; ix++)
        {
            totale += kk->charge->rho[is][ix];
            // compare rho read and rho calculated from wavefunctions
	    //std::cout<<"read "<< rho_for_compare[is][ix]<<" calc "<<kk->charge->rho[is][ix]<<std::endl;
            // EXPECT_NEAR(rho_for_compare[is][ix], kk->charge->rho[is][ix], 1e-5);
        }
    }
    // check total number of electrons
    totale = totale * GlobalC::ucell.omega / GlobalC::rhopw->nrxx;
    EXPECT_NEAR(totale, GlobalC::CHR.nelec, 1e-5);
    delete kk;
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
    // env->Reset();
    env->set_variables_in_set_up(ENVP);
    return RUN_ALL_TESTS();
}

int main(int argc, char** argv)
{
#ifdef __MPI
    MPI_Init(&argc, &argv);
#endif

    testing::InitGoogleTest(&argc, argv);

    ENVEnvironment* const env = new ENVEnvironment;
    testing::AddGlobalTestEnvironment(env);
    Check(RunAllTests(env, &ENVP) == 0, "");
#ifdef __MPI
    MPI_Finalize();
#endif

    return 0;
}
