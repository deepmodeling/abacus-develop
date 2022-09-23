#include "module_elecstate/elecstate_lcao.h"
#include "src_parallel/parallel_reduce.h"
#include "src_pdiag/test/diago_elpa_utils.h"
#include "gtest/gtest.h"
#include "mpi.h"

#define THRESHOLD 1e-10

/************************************************
 *  unit test of ElecStateLCAO
 ***********************************************/

/**
 * Test function psiToRho(), which is used to calculate 
 * desity matrix from wave function.
 * Containing the tests for gamma_only and multiple k points
 * cases.
 * 
 */

#include "src_parallel/parallel_grid.h"
#include "src_pw/wavefunc.h"
#include "src_pw/VNL_in_pw.h"
#include "src_pw/energy.h"
#include "src_pw/structure_factor.h"
#include "module_neighbor/sltk_atom_arrange.h"
#include "module_pw/pw_basis_k.h"
#include "module_xc/xc_functional.h"
#include "module_xc/exx_global.h"

Magnetism::Magnetism(){}
Magnetism::~Magnetism(){}
K_Vectors::K_Vectors(){}
K_Vectors::~K_Vectors(){}
Structure_Factor::Structure_Factor(){}
Structure_Factor::~Structure_Factor(){}
LCAO_Hamilt::LCAO_Hamilt(){}
LCAO_Hamilt::~LCAO_Hamilt(){}
LCAO_gen_fixedH::LCAO_gen_fixedH(){}
LCAO_gen_fixedH::~LCAO_gen_fixedH(){}
wavefunc::wavefunc(){}
wavefunc::~wavefunc(){}
WF_atomic::WF_atomic(){}
WF_atomic::~WF_atomic(){}
WF_igk::WF_igk(){}
WF_igk::~WF_igk(){}
pseudopot_cell_vnl::pseudopot_cell_vnl(){}
pseudopot_cell_vnl::~pseudopot_cell_vnl(){}
pseudopot_cell_vl::pseudopot_cell_vl(){}
pseudopot_cell_vl::~pseudopot_cell_vl(){}
energy::energy(){}
energy::~energy(){}
namespace GlobalC
{
    energy en;
    Structure_Factor sf;
    K_Vectors kv;
    UnitCell_pseudo ucell;
    pseudopot_cell_vnl ppcell;
    ModulePW::PW_Basis* rhopw;
    ModulePW::PW_Basis_Big *bigpw = static_cast<ModulePW::PW_Basis_Big*>(rhopw);
    ModulePW::PW_Basis_K* wfcpw;
    wavefunc wf;
    Charge CHR;
    Grid_Driver GridD(GlobalV::test_deconstructor, GlobalV::test_grid_driver,GlobalV::test_grid);
}


XC_Functional::XC_Functional(){}
XC_Functional::~XC_Functional(){}
int XC_Functional::get_func_type(){return 0;}

#ifdef __MPI
#include "src_ri/exx_lcao.h"
Exx_Lcao::Exx_Info::Exx_Info( const Exx_Global::Exx_Info &info_global )
    :hybrid_type(info_global.hybrid_type),hse_omega(info_global.hse_omega){}
Exx_Lcao::Exx_Lcao(const Exx_Global::Exx_Info &info_global ):info(info_global){}
namespace GlobalC
{
    Exx_Global exx_global;
    Exx_Lcao exx_lcao(GlobalC::exx_global.info); 
}
#endif

namespace WF_Local
{
    int read_lowf(double** ctot, const int& is, const Parallel_Orbitals* ParaV, psi::Psi<double>*) {return 1;};
    int read_lowf_complex(std::complex<double>** ctot, const int& ik, const Parallel_Orbitals* ParaV, psi::Psi<std::complex<double> >*) {return 1;}
    void write_lowf(const std::string &name, double **ctot, const ModuleBase::matrix& ekb, const ModuleBase::matrix& wg) {}
    void write_lowf_complex(const std::string &name, std::complex<double>** ctot, const int &ik, const ModuleBase::matrix& ekb, const ModuleBase::matrix& wg) {}
}

//mock the unrelated functions in charge.cpp
#include "src_pw/use_fft.h"
#include "src_pw/occupy.h"
namespace GlobalC {Use_FFT UFFT;}
Use_FFT::Use_FFT(){}
Use_FFT::~Use_FFT(){}
void Use_FFT::ToRealSpace(const int &is, const ModuleBase::ComplexMatrix &vg, double *vr, ModulePW::PW_Basis* rho_basis) {return;}
void Use_FFT::ToRealSpace(const std::complex<double> *vg, double *vr, ModulePW::PW_Basis* rho_basis) {return;};
bool Occupy::use_gaussian_broadening = false;
bool Occupy::use_tetrahedron_method = false;
double Magnetism::get_nelup(void) {return 0;}
double Magnetism::get_neldw(void) {return 0;}

void set_pw()
{
    GlobalC::rhopw->nx = 36; //should be only divided by 2,3,5
    GlobalC::rhopw->ny = 36;
    GlobalC::rhopw->nz = 36;
    GlobalC::bigpw->bx = 2; 
    GlobalC::bigpw->by = 2; 
    GlobalC::bigpw->bz = 2;
	GlobalC::bigpw->nbx = GlobalC::rhopw->nx/GlobalC::bigpw->bx; 
    GlobalC::bigpw->nby = GlobalC::rhopw->nx/GlobalC::bigpw->bx; 
    GlobalC::bigpw->nbz = GlobalC::rhopw->nx/GlobalC::bigpw->bx;
    GlobalC::bigpw->bxyz = GlobalC::bigpw->bx * GlobalC::bigpw->by * GlobalC::bigpw->bz;
    GlobalC::rhopw->nxyz = GlobalC::rhopw->nx * GlobalC::rhopw->ny * GlobalC::rhopw->nz;
	 
    if (GlobalV::MY_RANK < (GlobalC::bigpw->nbz % GlobalV::DSIZE))
    {
        GlobalC::bigpw->nbzp = GlobalC::bigpw->nbz/GlobalV::DSIZE + 1;
        GlobalC::bigpw->nbzp_start = (GlobalC::bigpw->nbz/GlobalV::DSIZE + 1) * GlobalV::MY_RANK;
    }
    else
    {
        GlobalC::bigpw->nbzp = GlobalC::bigpw->nbz/GlobalV::DSIZE;
        GlobalC::bigpw->nbzp_start = (GlobalC::bigpw->nbz/GlobalV::DSIZE) * GlobalV::MY_RANK + (GlobalC::bigpw->nbz % GlobalV::DSIZE); 
    }
    GlobalC::rhopw->nplane = GlobalC::bigpw->nbzp * GlobalC::bigpw->bz;
    GlobalC::rhopw->startz_current = GlobalC::bigpw->nbzp_start * GlobalC::bigpw->bz;
    GlobalC::bigpw->nbxx = GlobalC::bigpw->nbzp*GlobalC::bigpw->nbx*GlobalC::bigpw->nby;
    GlobalC::rhopw->nrxx = GlobalC::rhopw->nplane*GlobalC::rhopw->nx*GlobalC::rhopw->ny;
}

void init()
{
#ifdef __MPI
    MPI_Comm_size(MPI_COMM_WORLD,&GlobalV::NPROC);
	MPI_Comm_rank(MPI_COMM_WORLD,&GlobalV::MY_RANK); 
    Parallel_Global::split_diag_world(GlobalV::NPROC);
    Parallel_Global::split_grid_world(GlobalV::NPROC);
    MPI_Comm_split(MPI_COMM_WORLD,0,1,&POOL_WORLD); //in LCAO kpar=1
#endif

    GlobalV::BASIS_TYPE = "lcao";
    GlobalV::stru_file = "./support/si2.STRU";
    GlobalV::CALCULATION = "scf";
    GlobalV::global_pseudo_dir = "./support/";
    GlobalV::global_orbital_dir = "./support/";
    GlobalV::PSEUDORCUT = 15.0;
    GlobalC::ucell.ntype = 1;
    GlobalC::ORB.ecutwfc = 50;
    GlobalC::ORB.dk = 0.01;
    GlobalC::ORB.dR = 0.01;
    GlobalC::ORB.Rmax = 30;
    GlobalC::ORB.dr_uniform = 0.001;
    GlobalV::KS_SOLVER = "genelpa";
    GlobalV::NSPIN = 1;
    GlobalC::wf.init_wfc="atomic";
    GlobalC::rhopw = new ModulePW::PW_Basis_Big();
    GlobalC::bigpw = static_cast<ModulePW::PW_Basis_Big*>(GlobalC::rhopw);
    GlobalC::wfcpw = new ModulePW::PW_Basis_K_Big(); 
    ModulePW::PW_Basis_K_Big* tmp2 = static_cast<ModulePW::PW_Basis_K_Big*>(GlobalC::wfcpw);
    tmp2->setbxyz(GlobalC::bigpw->bx,GlobalC::bigpw->by,GlobalC::bigpw->bz);

    //GlobalC::ucell.setup(INPUT.latname, INPUT.ntype, INPUT.lmaxmax, INPUT.init_vel, INPUT.fixed_axes);
    GlobalC::ucell.setup("test", 1, 2, false, "None");
    GlobalC::ucell.setup_cell(GlobalC::ORB, GlobalV::global_pseudo_dir, GlobalV::stru_file, GlobalV::ofs_running);
    GlobalC::CHR.cal_nelec();
    int out_mat_r = 0;
    GlobalC::ORB.Read_Orbitals(GlobalV::ofs_running,GlobalC::ucell.ntype,GlobalC::ucell.lmax,GlobalV::deepks_setorb,
                                out_mat_r,GlobalV::CAL_FORCE,GlobalV::MY_RANK);
    ModuleBase::Ylm::set_coefficients();   
    set_pw();
    GlobalC::CHR.allocate(GlobalV::NSPIN, GlobalC::rhopw->nrxx, 0);                             
}

namespace elecstate
{
    class MockElecStateLCAO : public ElecStateLCAO
    {
        public:
        MockElecStateLCAO(Charge* chg_in,
                      const K_Vectors* klist_in,
                      int nks_in,
                      int nbands_in,
                      Local_Orbital_Charge* loc_in,
                      LCAO_Hamilt* uhm_in,
                      Local_Orbital_wfc* lowf_in,
                      ModuleBase::matrix &wg_in)
                      :elecstate::ElecStateLCAO(chg_in,klist_in,nks_in,nbands_in,loc_in,uhm_in,lowf_in)
                      {
                          this->wg = wg_in;
                      }               
    };

    const double* ElecState::getRho(int spin) const
    {
        return &(this->charge->rho[spin][0]);
    } 

    void ElecState::calculate_weights(void) {}
    void ElecState::calEBand() {}
    void ElecState::allocateRho(const int &nspin_in, const int &nrxx_in, const int &ngmc_in){}
}

template<class T>
class ElecStateLCAOPrepare
{
    public:
    ElecStateLCAOPrepare(int nk, int nbands, int nbasis)
    :nk(nk),nbands(nbands), nbasis(nbasis)
    {
        if(std::is_same<T, double>::value)
        {
            GlobalV::GAMMA_ONLY_LOCAL = true;
            nk = 1;
        }
        else
        {
            GlobalV::GAMMA_ONLY_LOCAL = false;
        }
        this->wg.create(nk, nbands);
    }
    ~ElecStateLCAOPrepare()
    {
        for(int is = 0; is < GlobalV::NSPIN; is++) 
        {
            delete [] this->rho_ref[is];
            delete [] this->rho_cal[is];
        }
        delete [] this->rho_ref;
        delete [] this->rho_cal;
    }

    int nk, nbands, nbasis, nelec;
    ModuleBase::matrix wg;
    psi::Psi<T> psi;
    psi::Psi<T> psi_local;
    double** rho_ref;
    double** rho_cal;

    void read_psi()
    {
        psi.resize(nk,nbands,nbasis);
        for(int i=1;i<=nk;i++)
        {
            std::string fname;
            if(std::is_same<T, double>::value)
                fname = std::string("./support/LOWF_GAMMA_S");
            else
                fname = std::string("./support/LOWF_K_");
            fname += to_string(i) + std::string(".dat");
            std::ifstream inf(fname);

            if(! inf.is_open())
            {
                std::cout << "Error: open file " << fname << " failed, skip!" << std::endl;
                exit(1);
            }

            string stmp;
            if(std::is_same<T, std::complex<double>>::value)
            {
                getline(inf,stmp);
                getline(inf,stmp);
            }

            int f_nbasis,f_nbands;
            inf >> f_nbands; getline(inf,stmp);
            inf >> f_nbasis; getline(inf,stmp);
            if(f_nbasis != nbasis || f_nbands != nbands)
            {
                std::cout << "Error: nbasis in file " << fname << " is " << f_nbasis << ", should be " << nbasis << std::endl;
                std::cout << "Error: nbands in file " << fname << " is " << f_nbands << ", should be " << nbands << std::endl;
                exit(1);
            }

            int band_local=0;
            std::string bandline("(band)");
            while(getline(inf,stmp))
            {
                if (stmp.find(bandline) != std::string::npos)
                {
                    getline(inf,stmp);
                    inf >> wg(i-1,band_local) >> stmp;
                    for(int j=0;j<f_nbasis;j++) 
                    {
                        inf >> psi(i-1,band_local,j);
                        if(std::is_same<T, std::complex<double>>::value)
                            inf >> reinterpret_cast<double *>(&psi(i-1,band_local,j))[1];
                    }
                    band_local += 1;
                } 
            }
        }
    }

    void set_env()
    {
        GlobalV::NBANDS = nbands;
        GlobalC::wf.wg = this->wg;

        GlobalC::kv.nks = GlobalC::kv.nkstot = nk;
        GlobalC::kv.isk.resize(nk,0);
        GlobalC::kv.kvec_d.resize(nk);
        for(int i=0;i<nk;i++)
        {
            GlobalC::kv.kvec_d[i].x = 0.0 + double(i)/double(nk);
            GlobalC::kv.kvec_d[i].y = 0.0;
            GlobalC::kv.kvec_d[i].z = 0.0;
        }
    }

    void distribute_psi_2d(psi::Psi<T> &psi, psi::Psi<T> &psi_local,const int* desc)
    {
        int nb2d = desc[4];
        int icontxt = desc[1];
        int nprows, npcols, myprow, mypcol;
        Cblacs_gridinfo(icontxt, &nprows, &npcols, &myprow, &mypcol);
        int nk = psi.get_nk();
        int nbasis = psi.get_nbasis();
        int nbands = psi.get_nbands();
        int local_basis = LCAO_DIAGO_TEST::na_rc(nbasis,nb2d,nprows,myprow);
        int local_bands = LCAO_DIAGO_TEST::na_rc(nbands,nb2d,npcols,mypcol);
        psi_local.resize(nk,local_bands,local_basis);
        T* matrix_all = new T[nbasis*nbasis];
        T* matrix_local = new T[local_basis*local_basis];
        ModuleBase::GlobalFunc::ZEROS(matrix_all, nbasis*nbasis);

        for (int k=0;k<nk;k++)
        {
            for(int i=0;i<nbasis;i++)
            {
                for(int j=0;j<nbands;j++)
                {
                    matrix_all[i*nbasis + j] = psi(k,j,i);
                }
            }
            //matrix_all should be row-first, matrix_local will be column-first
            LCAO_DIAGO_TEST::distribute_data<T>(matrix_all,matrix_local,nbasis,nb2d,local_basis,local_basis,icontxt);
            for(int i=0;i<local_basis;i++)
            {
                for(int j=0;j<local_bands;j++)
                {
                    psi_local(k,j,i) = matrix_local[j*local_basis + i];
                }
            }
        }
    }

    void read_ref_rho()
    {
        //rho is parallel along ncz, but store is ncz_rank:ncy:ncx
        this->rho_ref = new double* [GlobalV::NSPIN];
        for(int is = 0; is < GlobalV::NSPIN; is++) this->rho_ref[is] = new double[GlobalC::rhopw->nxyz];
        if(GlobalV::MY_RANK == 0)
        {
            std::ifstream ifs;
            if(GlobalV::GAMMA_ONLY_LOCAL) ifs.open("./support/rho_gamma_ref.dat");
            else ifs.open("./support/rho_k_ref.dat");

            for(int is = 0; is < GlobalV::NSPIN; is++)
            {
                for(int i=0;i<GlobalC::rhopw->nxyz;i++)
                {
                    ifs >> this->rho_ref[is][i];
                }
            }
            ifs.close();
        }
    }

    void gather_rho(Charge* charge)
    {
        this->rho_cal = new double* [GlobalV::NSPIN];
        for(int is = 0; is < GlobalV::NSPIN; is++) this->rho_cal[is] = new double[GlobalC::rhopw->nxyz];

#ifdef __MPI
        MPI_Status ierror;
        int ncx = GlobalC::rhopw->nx;
        int ncy = GlobalC::rhopw->ny;
        int ncz = GlobalC::rhopw->nz;
        int nczp = GlobalC::rhopw->nplane;

        int* nz = new int[GlobalV::NPROC];
        for(int i=0;i<GlobalV::NPROC;i++)
        {
            nz[i] = GlobalC::bigpw->nbz/GlobalV::NPROC;
            if (i < (GlobalC::bigpw->nbz % GlobalV::NPROC)) nz[i] += 1;
            nz[i] *= GlobalC::bigpw->bz;
        }

        for(int is = 0; is < GlobalV::NSPIN; is++)
        {
            for(int i=0;i<ncx;i++)
            {
                for(int j=0;j<ncy; j++)
                {
                    if(GlobalV::MY_RANK != 0)
                    {
                        
                        MPI_Send(&charge->rho[is][i*ncy*nczp+j*nczp], nczp, MPI_DOUBLE, 0, 0, POOL_WORLD);
                    }
                    else
                    {
                        int position = nz[0];
                        for(int k=1;k<GlobalV::NPROC;k++)
                        {
                            MPI_Recv(&rho_cal[is][i*ncy*ncz+j*ncz+position], nz[k],MPI_DOUBLE,k,0,POOL_WORLD, &ierror);
                            position += nz[k];
                        }
                        for(int k=0;k<nczp;k++)
                            rho_cal[is][i*ncy*ncz+j*ncz+k] = charge->rho[is][i*ncy*nczp+j*nczp+k];
                    }
                }
            }
        }
        delete [] nz;
#else
        for(int is = 0; is < GlobalV::NSPIN; is++)
        {
            for(int i=0;i<GlobalC::rhopw->nxyz;i++)
                rho_cal[is][i] = charge->rho[is][i];            
        }
#endif
    }

    void out_rho(Charge* charge)
    {
        //use only one core to output the data.
        std::ofstream ofs("rho_ref.dat"); 
        if(GlobalV::MY_RANK == 0)
        {
            for (int is = 0; is < GlobalV::NSPIN; is++)
            {
                for(int i=0;i<charge->nrxx;i++) 
                    ofs << std::setiosflags(std::ios_base::scientific) << std::setprecision(15) << charge->rho[is][i] << " ";
                ofs << std::endl;
            }
        }
        ofs.close();
    }

    void run()
    {
        this->set_env();

        Local_Orbital_Charge loc;
        LCAO_Hamilt uhm;
        Local_Orbital_wfc lowf;

        ORB_control orb_con(GlobalV::GAMMA_ONLY_LOCAL, GlobalV::NLOCAL, GlobalV::NBANDS, 
                            GlobalV::NSPIN, GlobalV::DSIZE, GlobalV::NB2D, GlobalV::DCOLOR, 
                            GlobalV::DRANK, GlobalV::MY_RANK, GlobalV::CALCULATION, GlobalV::KS_SOLVER);
        orb_con.setup_2d_division(GlobalV::ofs_running, GlobalV::ofs_warning);

        loc.ParaV = lowf.ParaV = &(orb_con.ParaV);
#ifdef __MPI        
        this->distribute_psi_2d(this->psi,this->psi_local,loc.ParaV->desc_wfc);
#else
        this->psi_local = this->psi;
#endif        
        GlobalV::SEARCH_RADIUS = atom_arrange::set_sr_NL(
            GlobalV::ofs_running,
            GlobalV::OUT_LEVEL,
            GlobalC::ORB.get_rcutmax_Phi(),
            GlobalC::ucell.infoNL.get_rcutmax_Beta(),
            GlobalV::GAMMA_ONLY_LOCAL);

        atom_arrange::search(
            GlobalV::SEARCH_PBC,
            GlobalV::ofs_running,
            GlobalC::GridD,
            GlobalC::ucell,
            GlobalV::SEARCH_RADIUS,
            GlobalV::test_atom_input);

        GlobalC::GridT.set_pbc_grid(
			GlobalC::rhopw->nx, GlobalC::rhopw->ny, GlobalC::rhopw->nz,
			GlobalC::bigpw->bx, GlobalC::bigpw->by, GlobalC::bigpw->bz,
			GlobalC::bigpw->nbx, GlobalC::bigpw->nby, GlobalC::bigpw->nbz,
			GlobalC::bigpw->nbxx, GlobalC::bigpw->nbzp_start, GlobalC::bigpw->nbzp);
        if (!GlobalV::GAMMA_ONLY_LOCAL)
        {
            GlobalC::GridT.cal_nnrg(&(orb_con.ParaV));
        }
        psi::Psi<double>* psigo = nullptr;
        psi::Psi<std::complex<double>>* psik = nullptr;
        if(GlobalV::GAMMA_ONLY_LOCAL)
        {
            psigo = (psi::Psi<double>*)(&(this->psi));
        }
        else
        {
            psik = (psi::Psi<std::complex<double>>*)(&(this->psi));
        }
        loc.allocate_dm_wfc(GlobalC::GridT.lgd, lowf, psigo, psik);

        elecstate::MockElecStateLCAO mesl(&GlobalC::CHR,&GlobalC::kv,nk,nbands,&loc,&uhm,&lowf,this->wg);

        mesl.psiToRho(psi_local);  
    }

    void compare()
    {
        //this->out_rho(&GlobalC::CHR);       
        this->gather_rho(&GlobalC::CHR);          
        this->read_ref_rho();
        if(GlobalV::MY_RANK == 0)
        {
            bool pass = true;
            double error; 
            for (int is = 0; is < GlobalV::NSPIN; is++)
            {
                for(int i=0;i<GlobalC::rhopw->nxyz;i++)
                {
                    error = abs(this->rho_cal[is][i] - this->rho_ref[is][i]);
                    if ( error > THRESHOLD)
                    {
                        pass = false;
                        break;
                    }
                } 
            }
            EXPECT_TRUE(pass) << "error is " << error;
        }
    }
};

TEST(LCAOPsiToRho,GammaOnly)
{
    ElecStateLCAOPrepare<double> eslp(1,6,26);
    eslp.read_psi();
    eslp.run();
    eslp.compare();
}

TEST(LCAOPsiToRho,MultipleK)
{
    ElecStateLCAOPrepare<std::complex<double>> eslp(3,6,26);
    eslp.read_psi();
    eslp.run();
    eslp.compare();
}

int main(int argc, char **argv)
{
    MPI_Init(&argc, &argv);
    testing::InitGoogleTest(&argc, argv);

    init();

    ::testing::TestEventListeners &listeners = ::testing::UnitTest::GetInstance()->listeners();
    if (GlobalV::MY_RANK != 0)
    {
        delete listeners.Release(listeners.default_result_printer());
    }

    int result = RUN_ALL_TESTS();
    delete GlobalC::rhopw;
    delete GlobalC::wfcpw;
    if (GlobalV::MY_RANK == 0 && result != 0)
    {
        std::cout << "ERROR:some tests are not passed" << std::endl;
        return result;
    }
    else
    {
        MPI_Finalize();
        return 0;
    }
}
