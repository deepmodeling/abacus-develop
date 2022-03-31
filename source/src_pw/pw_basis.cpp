// MPI BUG FIX STEPS:
// 1. add GlobalC::pw.nrxx_mpi ,used in sum_band_k , rho's Allgather
// 2. change ./src_parallel/ft.cpp 's some group size
// 3. delete condition justification in pw_and_columns distributins 2

#include "../module_base/global_function.h"
#include "../module_base/global_variable.h"
#include "../module_base/memory.h"
#include "../src_parallel/parallel_reduce.h"
#include "pw_basis.h"
#include "pw_complement.h"
#include "../module_base/timer.h"

#ifdef _OPENMP
#include <omp.h>
#endif

PW_Basis::PW_Basis()
{
//	std::cout << "\n PW_Basis " << std::endl;
    ig2fftw = nullptr;
    ig2fftc = nullptr;

	// if not parallel, gg_global == gg, gdirect_global == gdirect
	// gcar_global == gcar
    gg_global = nullptr;
    gdirect_global = nullptr;
    gcar_global = nullptr;

#ifdef __MPI
    gg = nullptr;
    gdirect = nullptr;
    gcar = nullptr;
#endif

    ggs = nullptr;
    ig2ngg = nullptr;

    ig1 = nullptr;
    ig2 = nullptr;
    ig3 = nullptr;

    this->nczp_start = 0;
    gg_global0 = nullptr; //LiuXh 20180515
    cutgg_num_table = nullptr; //LiuXh 0180515
    ggchg_time_global = 0; //LiuXh 20180515

}

PW_Basis::~PW_Basis()
{
	if(GlobalV::test_deconstructor)
	{
		std::cout << " ~PW_Basis()" << std::endl;
	}
    delete [] ig2fftw;
    delete [] ig2fftc;

    delete [] gcar_global;
    delete [] gdirect_global;
    delete [] gg_global;
    delete [] gg_global0; //LiuXh 20180515
	delete [] cutgg_num_table;

#ifdef __MPI
#if ((!defined __CUDA) && (!defined __ROCM))
    delete [] gcar;
    delete [] gdirect;
    delete [] gg;
#endif
#endif

    delete [] ggs;
    delete [] ig2ngg;

    delete [] ig1;
    delete [] ig2;
    delete [] ig3;
}

// called in input.cpp
void PW_Basis::set
(
    const bool &gamma_only_in,
    const double &ecutwfc_in,
    const double &ecutrho_in,
    const int &nx_in,
    const int &ny_in,
    const int &nz_in,
    const int &ncx_in,
    const int &ncy_in,
    const int &ncz_in,
	const int &bx_in,
	const int &by_in,
	const int &bz_in,
    const int &seed_in,
    const int &nbspline_in
)
{
    ModuleBase::TITLE("PW_Basis","set");
    this->gamma_only = gamma_only_in;
    this->ecutwfc = ecutwfc_in;
    this->ecutrho = ecutrho_in,
    this->nx = nx_in;
    this->ny = ny_in;
    this->nz = nz_in;
    this->ncx = ncx_in;
    this->ncy = ncy_in;
    this->ncz = ncz_in;
	this->bx = bx_in;
	this->by = by_in;
	this->bz = bz_in;
    this->pw_seed = seed_in;
    this->nbspline = nbspline_in;

    if (ecutwfc <= 0.00)
    {
        ModuleBase::WARNING_QUIT("PW_Basis::set","ecutwfc < 0 is not allowed !");
    }

    if (ecutrho <= 0.00)
    {
        this->wfac = 4.0;
    }
    else
    {
        this->wfac = ecutrho/ecutwfc;
        if (wfac <= 1.0)
        {
            ModuleBase::WARNING_QUIT("input","GlobalC::pw.wfac <= 1.0 is not allowed !");
        }
    }
    return;
}


// initialize of plane wave basis.
void PW_Basis::gen_pw(std::ofstream &runlog, const UnitCell &Ucell_in, const K_Vectors &Klist_in)
{
    ModuleBase::TITLE("PW_Basis","gen_pw");
    ModuleBase::timer::tick("PW_Basis","gen_pw");


	GlobalV::ofs_running << "\n\n\n\n";
	GlobalV::ofs_running << " >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>" << std::endl;
	GlobalV::ofs_running << " |                                                                    |" << std::endl;
	GlobalV::ofs_running << " | Setup plane waves:                                                 |" << std::endl;
	GlobalV::ofs_running << " | Use the energy cutoff and the lattice vectors to generate the      |" << std::endl;
	GlobalV::ofs_running << " | dimensions of FFT grid. The number of FFT grid on each processor   |" << std::endl;
	GlobalV::ofs_running << " | is 'nrxx'. The number of plane wave basis in reciprocal space is   |" << std::endl;
	GlobalV::ofs_running << " | different for charege/potential and wave functions. We also set    |" << std::endl;
	GlobalV::ofs_running << " | the 'sticks' for the parallel of FFT.                              |" << std::endl;
	GlobalV::ofs_running << " |                                                                    |" << std::endl;
	GlobalV::ofs_running << " <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<" << std::endl;
	GlobalV::ofs_running << "\n\n\n\n";



	GlobalV::ofs_running << "\n SETUP THE PLANE WAVE BASIS" << std::endl;

    // Now set up a few parameters that related to Ecutwfc
    // read Ecutwfc from the parameter table in input_parametes class

	this->Ucell = &Ucell_in;
	this->Klist = &Klist_in;

    //mohan modify 2008-3-25
	// setup for ggpsi, ggwfc, ggwfc2, ggchg
    this->setup_gg();
    this->setup_FFT_dimension();

	//----------------------------------------
	// if calculation is 'test', we don't need
	// to allocate the arrays, claculate the
	// structure factors, etc.
	// just return here.
	//----------------------------------------
	if(GlobalV::CALCULATION=="test")
	{
		return;
	}

#ifdef __MPI
#if ((defined __CUDA) || (defined __ROCM))
    bool cutgg_flag = false;
#else
    this->divide_fft_grid();
    bool cutgg_flag = true;

	// mohan add 2011-07-23
	//if(winput::out_spillage)
	//{
	//	cutgg_flag = false;
	//	GlobalV::ofs_running << " Turn off the cutgg function." << std::endl;
	//}
	//else
	//{
	//	GlobalV::ofs_running << " Turn on the cutgg function." << std::endl;
	//}
#endif
#else
	// mohan update 2011-09-21
	this->nbzp=nbz; //nbz shoud equal nz for single proc.
	this->nczp=nbzp*bz; 
	this->nbxx=nbz*nbx*nby; //mohan fix 2021-02-17
	this->nbzp_start=0;
    bool cutgg_flag = false;
#endif

    if (cutgg_flag)
    {
#ifdef __MPI
        GlobalV::ofs_running << "\n SETUP COORDINATES OF PLANE WAVES" << std::endl;
		
		// get the number of total plane waves within a spheri.
        this->ngmc_g = PW_complement::get_total_pw_number(0.0, ggchg, ncx, ncy, ncz, Ucell->GGT);
		ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running,"number of total plane waves",ngmc_g);

        double cutgg_pieces = 10;

		// ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running,"cutgg_pieces",cutgg_pieces);
        const double cutgg_delta = ggchg / std::pow( (double)cutgg_pieces, 2.0/3.0 );

        // get cutgg_delta from input.
		// ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running,"cutgg_delta",cutgg_delta);

        //int cutgg_num_start = 0;
        double ggchg_start = 0.0;
        double ggchg_end = 0.0;

        int ggchg_time = 1;
        do
        {
            // get the cutoff now.
            // 2.0/3.0 factor: because gg = |g|^{2},
            // we ensure increasing the same volumn each time
            // to make sure to get the nearly same plane wave number.
            ggchg_end = cutgg_delta * std::pow( (double)ggchg_time, 2.0/3.0 ) ;
            if ( abs(ggchg_end - ggchg)<1.0e-8 ) ggchg_end = ggchg;
            //std::cout << "\n ggchg_start = " << ggchg_start;
            //std::cout << "\n ggchg_end = " << ggchg_end << std::endl;

            // get the number of plane waves between two shells.
            const int cutgg_num_now = PW_complement::get_total_pw_number( ggchg_start, ggchg_end, ncx, ncy, ncz, Ucell->GGT);
            //std::cout << "\n cutgg_num_now = " << cutgg_num_now << std::endl;

            delete[] gg_global;
            delete[] gdirect_global;
            delete[] gcar_global;
            gg_global = new double[cutgg_num_now];
            gdirect_global = new ModuleBase::Vector3<double>[cutgg_num_now];
            gcar_global = new ModuleBase::Vector3<double>[cutgg_num_now];

			//GlobalV::ofs_running << " setup |g|^2" << std::endl;
            PW_complement::get_total_pw(gg_global, gdirect_global, ggchg_start, ggchg_end,
                                        ncx, ncy, ncz, Ucell->GGT, ngmc_g);
            
			PW_complement::setup_GVectors(Ucell->G, cutgg_num_now, gg_global, gdirect_global, gcar_global);

            //FFT_chg.fft_map(this->ig2fftc, this->ngmc, cutgg_num_now);
            //FFT_wfc.fft_map(this->ig2fftw, this->ngmw, cutgg_num_now);
            //FFT_chg.fft_map(this->ig2fftc, this->ngmc, cutgg_num_now, ggchg_time);
            //FFT_wfc.fft_map(this->ig2fftw, this->ngmw, cutgg_num_now, ggchg_time);
            if(!GlobalV::FINAL_SCF) //LiuXh add 20180619
            {
                FFT_chg.fft_map(this->ig2fftc, this->ngmc, cutgg_num_now, ggchg_time);
                FFT_wfc.fft_map(this->ig2fftw, this->ngmw, cutgg_num_now, ggchg_time);
            }
            else
            {
                FFT_chg.fft_map_final_scf(this->ig2fftc, this->ngmc, cutgg_num_now);
                FFT_wfc.fft_map_final_scf(this->ig2fftw, this->ngmw, cutgg_num_now);
            }

            ggchg_start = ggchg_end;
            ++ggchg_time;

//			std::cout << "\n ggchg_end = " << ggchg_end;
//			std::cout << "\n ggchg = " << ggchg;
        } while ( ggchg_end < ggchg );

        if(!GlobalV::FINAL_SCF) //LiuXh add 20180619
        {
            //LiuXh add 20180515, begin
            ggchg_time_global = ggchg_time;
			delete [] cutgg_num_table;
            cutgg_num_table = new int[ggchg_time_global];
            ModuleBase::GlobalFunc::ZEROS(cutgg_num_table, ggchg_time_global);
            GlobalV::ofs_running << "\n SETUP COORDINATES OF PLANE WAVES" << std::endl;

            double cutgg_pieces2 = 10;

            const double cutgg_delta2 = ggchg / std::pow( (double)cutgg_pieces2, 2.0/3.0 );

            //int cutgg_num_start2 = 0;
            double ggchg_start2 = 0.0;
            double ggchg_end2 = 0.0;

            int ggchg_time2 = 1;
            do
            {
                ggchg_end2 = cutgg_delta2 * std::pow( (double)ggchg_time2, 2.0/3.0 ) ;
                if ( ggchg_end2 > ggchg ) ggchg_end2 = ggchg;

                const int cutgg_num_now2 = PW_complement::get_total_pw_number( ggchg_start2, ggchg_end2, ncx, ncy, ncz, Ucell->GGT);
                cutgg_num_table[ggchg_time2-1] = cutgg_num_now2;

                ggchg_start2 = ggchg_end2;
                ++ggchg_time2;
            } while ( ggchg_end2 < ggchg ); //LiuXh add 20180515, end
        }

#endif
        /*
        std::stringstream ss;
        ss << "./out_data/ig2fftc.dat" << cutgg_delta;
        std::ofstream ofsc( ss.str().c_str() );

        ofsc << "\n ig2fftc" << std::endl;
        for(int ig=0; ig<ngmc; ig++)
        {
        	ofsc << std::setw(10) << ig << std::setw(20) << ig2fftc[ig] << std::endl;
        }
        */

#ifdef __MPI
        this->get_MPI_GVectors();
        //std::cout<<" UNIFORM GRID DIM     : "<<this->nx <<" * "<<this->ny <<" * "<<this->nz  << "," << this->nrxx << std::endl;
        //std::cout<<" UNIFORM GRID DIM     : "<<this->ncx<<" * "<<this->ncy<<" * "<<this->ncz << "," << this->nrxx  << std::endl;

        FFT_wfc.setup_MPI_FFT3D(this->nx, this->ny, this->nz,this->nrxx,1);
        FFT_chg.setup_MPI_FFT3D(this->ncx, this->ncy, this->ncz,this->nrxx,1);
#endif
    }
    else
    {
        this->ngmc_g = PW_complement::get_total_pw_number(0.0, ggchg, ncx, ncy, ncz, Ucell->GGT);
		ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running,"ngmc_g",ngmc_g);

        delete[] gg_global;
        gg_global = new double[ngmc_g];// store the |G|^2 of the 1d array
        ModuleBase::Memory::record("PW_complement","gg_global",ngmc_g,"double");

        delete[] gdirect_global;
        gdirect_global = new ModuleBase::Vector3<double>[ngmc_g];// indices of G vectors
        ModuleBase::Memory::record("PW_complement","gdirect_global",ngmc_g,"ModuleBase::Vector3<double>");

        delete[] gcar_global;
        gcar_global = new ModuleBase::Vector3<double>[ngmc_g];
        ModuleBase::Memory::record("PW_complement","gcar",ngmc_g,"ModuleBase::Vector3<double>");

        PW_complement::get_total_pw(gg_global, gdirect_global, 0.0, ggchg, ncx, ncy, ncz, Ucell->GGT, ngmc_g);
        PW_complement::setup_GVectors(Ucell->G, ngmc_g, gg_global, gdirect_global, gcar_global);

#ifdef __MPI
#if ((defined __CUDA) || (defined __ROCM))
        this->get_GVectors();
        FFT_wfc.setupFFT3D(this->nx, this->ny,this->nz);
        FFT_chg.setupFFT3D(this->ncx, this->ncy,this->ncz);
#else
        //FFT_chg.fft_map(this->ig2fftc, this->ngmc, ngmc_g);
        //FFT_wfc.fft_map(this->ig2fftw, this->ngmw, ngmc_g);
        FFT_chg.fft_map(this->ig2fftc, this->ngmc, ngmc_g, 0); //LiuXh add 20180619
        FFT_wfc.fft_map(this->ig2fftw, this->ngmw, ngmc_g, 0); //LiuXh add 20180619

        this->get_MPI_GVectors();

        FFT_wfc.setup_MPI_FFT3D(this->nx, this->ny, this->nz,this->nrxx,1);
        FFT_chg.setup_MPI_FFT3D(this->ncx, this->ncy, this->ncz,this->nrxx,1);
#endif
#else
        this->get_GVectors();
        FFT_wfc.setupFFT3D(this->nx, this->ny,this->nz);
        FFT_chg.setupFFT3D(this->ncx, this->ncy,this->ncz);
#endif
		ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running,"fft charge/potential grid",ncx,ncy,ncz);
		ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running,"fft wavefunction grid",nx,ny,nz);
    }

    this->get_nggm(this->ngmc);

//	this->printPW("src_check/check_pw.txt");
    ModuleBase::timer::tick("PW_Basis","gen_pw");
    return;
}

void PW_Basis::setup_gg(void)
{
    ModuleBase::TITLE("PW_Basis","setup_gg");

    if (Ucell->tpiba2 <= 0)
    {
        ModuleBase::WARNING_QUIT("PW_Basis::setup_gg","tpiba2 <= 0");
    }
    this->ggpsi = this->ecutwfc / Ucell->tpiba2;
    //=================================
    // ecutwfc is the energy cutoff of
    // wave function (in Ry),
    // ggpsi wave function cut off for
    // Harmiltonian in the unit of
    // Ry*a0^2/(2*PI)^2
    //=================================

    // FFT cut off for wave function
    this->ggwfc = 4 * this->ggpsi;

    // FFT cut off for charge/potential
    this->ggchg = wfac * this->ggpsi;

    this->ggwfc2 = 0;

    if (gamma_only)
    {
        ggwfc2 = ggwfc;
    }
    else
    {
        for (int ik = 0;ik < this->Klist->nks;ik++)
        {
            const double k_mod = sqrt( Klist->kvec_c[ik] * Klist->kvec_c[ik]);
            const double tmp = sqrt(this->ggpsi) + k_mod;
            const double tmp2 = tmp * tmp ;
            if (this->ggwfc2 < tmp2) this->ggwfc2 = tmp2;
        }
        //qianrui add 2021-8-13 to make different kpar parameters can get the same results
#ifdef __MPI
        if(pw_seed > 0)    MPI_Allreduce(MPI_IN_PLACE, &ggwfc2, 1, MPI_DOUBLE, MPI_MAX , MPI_COMM_WORLD);
#endif
    }

	ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running,"energy cutoff for wavefunc (unit:Ry)",ecutwfc);

    return;
}


//  First stage Basis initialization.
//  Set up crystal structure parameters.
void PW_Basis::setup_FFT_dimension(void)
{
    if (GlobalV::test_pw) ModuleBase::TITLE("PW_Basis","setup_FFT_dimension");

    this->nxyz = nx * ny * nz;
    this->ncxyz = ncx * ncy * ncz;
	this->bxyz = bx * by * bz;

    if(GlobalV::FINAL_SCF) //LiuXh add 20180619
    {
        this->nxyz = 0;
        this->ncxyz = 0;
    }

    if (this->nxyz == 0)
    {
        PW_complement::get_FFT_dimension(Ucell->latvec, ggwfc, nx, ny, nz, bx, by, bz);
        this->nxyz = nx * ny * nz;
    }
	else
	{
		GlobalV::ofs_running << " use input fft dimensions for wave functions." << std::endl;
	}

    if (this->ncxyz == 0)
    {
        PW_complement::get_FFT_dimension(Ucell->latvec, ggchg, ncx, ncy, ncz, bx, by, bz);
        this->ncxyz = ncx * ncy * ncz;
    }
	else
	{
		GlobalV::ofs_running << " use input fft dimensions for charge/potential. " << std::endl;
	}

	assert(ncxyz >= nxyz);

	this->nbx = this->ncx / this->bx;
	this->nby = this->ncy / this->by;
	this->nbz = this->ncz / this->bz;
	this->nbxyz = this->nbx * this->nby * this->nbz;//mohan add 2011-06-02

	ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running,"fft grid for wave functions",nx,ny,nz);
	ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running,"fft grid for charge/potential",ncx,ncy,ncz);
	ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running,"fft grid division",bx,by,bz);
	ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running,"big fft grid for charge/potential",nbx,nby,nbz);

	// used in Grid_Base_Beta.init() (Grid integration)
	// mohan add 2009-11-09
	this->nczp = ncz;

    return;
}


#ifdef __MPI
void PW_Basis::divide_fft_grid(void)
{
    ModuleBase::TITLE("PW_Basis","divide_fft_grid");

    //----------------------------------------------
    // set charge/potential grid : nrxx
    //----------------------------------------------
    const int remain_planes = this->nbz%GlobalV::NPROC_IN_POOL;
    this->nbzp = this->nbz/GlobalV::NPROC_IN_POOL;
    
	if (GlobalV::RANK_IN_POOL < remain_planes)
    {
        nbzp++;
    }

	nczp = nbzp * bz;

	this->nbxx = nbzp*this->nbx*this->nby;
    this->nrxx = nczp*this->ncx*this->ncy;

    //=====================================
    // set nrxx_start
    //=====================================
    if (GlobalV::RANK_IN_POOL < remain_planes)
    {
        this->nrxx_start = GlobalV::RANK_IN_POOL * nrxx;
		this->nczp_start = GlobalV::RANK_IN_POOL * nczp;
		this->nbzp_start = GlobalV::RANK_IN_POOL * nbzp;
    }
    else
    {
        this->nrxx_start = GlobalV::RANK_IN_POOL * nrxx + remain_planes * bz * ncx * ncy;
		this->nczp_start = GlobalV::RANK_IN_POOL * nczp + remain_planes * bz;
		this->nbzp_start = GlobalV::RANK_IN_POOL * nbzp + remain_planes;
    }

	if(GlobalV::test_pw)ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running,"remain planes",remain_planes);
	if(GlobalV::test_pw)ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running,"small planes in this processor(charge)",nczp);
	if(GlobalV::test_pw)ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running,"big   planes in this processor(charge)",nbzp);
	if(GlobalV::test_pw)ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running,"nczp_start",nczp);
	if(GlobalV::test_pw)ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running,"nbzp_start",nbzp);
	ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running,"nbxx",nbxx);
	ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running,"nrxx",nrxx);
	if(GlobalV::test_pw)ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running,"nrxx_start",nrxx_start);

    //=====================================
    // generate nst,st_i,st_j,st_k,npps
    //=====================================
	
	GlobalV::ofs_running << "\n SETUP PLANE WAVES FOR CHARGE/POTENTIAL" << std::endl;
    FFT_chg.init(this->ggchg, this->ncx, this->ncy, this->ncz, this->bz, GlobalV::NPROC_IN_POOL,GlobalV::RANK_IN_POOL);
    FFT_chg.columns_map();
    FFT_chg.restore_st();
    FFT_chg.columns_and_pw_distribution();

	GlobalV::ofs_running << "\n SETUP PLANE WAVES FOR WAVE FUNCTIONS" << std::endl;
    FFT_wfc.init(this->ggwfc2, this->nx, this->ny, this->nz, this->bz, GlobalV::NPROC_IN_POOL,GlobalV::RANK_IN_POOL);
    FFT_wfc.columns_map();
    FFT_wfc.restore_st();
    FFT_wfc.columns_and_pw_distribution();
    //=====================================
    this->columns_and_pw_distribution_2();
    this->ngmc = FFT_chg.npw_per[GlobalV::RANK_IN_POOL];
	this->ngmw = FFT_wfc.npw_per[GlobalV::RANK_IN_POOL];
	//=====================================
	// generate isind,ismap,st_start
    // isind is used in fft_map
    // ismap is used in FFT.
    // st_start is used in FFT.
    FFT_chg.fft_dlay_set();
    FFT_wfc.fft_dlay_set();
    //=====================================
    //generate ig2fftc and ig2fftw
    delete[] ig2fftc;
    delete[] ig2fftw;
    delete[] gdirect;
    delete[] gcar;
    delete[] gg;
    this->ig2fftc = new int[ngmc];
    this->ig2fftw = new int[ngmw];
    this->gdirect = new ModuleBase::Vector3<double>[ngmc];
    this->gcar  = new ModuleBase::Vector3<double>[ngmc];
    this->gg = new double[ngmc];
    return;
}


//////////////////////////////////  EXPLAIN    //////////////////////////////////
//  M. Gong has made mistakes during calculating the effective structure factor
//  ig1[i], ig2[i], ig3[i] store the $G$ points using Cartesian coordinate,
//  where  -ncx <= ig1[] <= ncx, -ncy <= ig2[] <= ncy and -ncz <= ig3[] <= ncz
//  ngmc > ngmw. ig1, ig2 and ig3 is a mapping between (k+G) <-> G
//////////////////////////////////  EXPLAIN    //////////////////////////////////


void PW_Basis::get_MPI_GVectors(void)
{
    if (GlobalV::test_pw) ModuleBase::TITLE("PW_Basis","get_MPI_GVectors");

    delete[] ig1;
    delete[] ig2;
    delete[] ig3;
    this->ig1 = new int[ngmc];
    this->ig2 = new int[ngmc];
    this->ig3 = new int[ngmc];
	ModuleBase::GlobalFunc::ZEROS(ig1, ngmc);
	ModuleBase::GlobalFunc::ZEROS(ig2, ngmc);
	ModuleBase::GlobalFunc::ZEROS(ig3, ngmc);

    for (int i = 0; i < ngmc;i++)
    {
        this->ig1[i] = int(this->gdirect[i].x) + ncx;
        this->ig2[i] = int(this->gdirect[i].y) + ncy;
        this->ig3[i] = int(this->gdirect[i].z) + ncz;
    }

//	for( int i = 0; i < ngmc; ++i)
//	{
//		std::cout << gcar[i].x << " " << gcar[i].y << " " << gcar[i].z << std::endl;
//	}

    //=====================================
    if (GlobalV::MY_RANK==0)
    {
//        std::stringstream ssc,ssw;
//        ssc << GlobalV::global_out_dir << "test.FFT_chg";
//        ssw << GlobalV::global_out_dir << "test.FFT_wfc";
//        std::ofstream charge_data( ssc.str().c_str() );
//        std::ofstream wavefun_data( ssw.str().c_str() );
//        FFT_chg.print_data(charge_data);
//        FFT_wfc.print_data(wavefun_data);
//        charge_data.close();
//        wavefun_data.close();
    }
    //=====================================
}//end setup_mpi_GVectors
// #else
#endif
void PW_Basis::get_GVectors(void)
{
    if (GlobalV::test_pw) ModuleBase::TITLE("PW_Basis","get_GVectors");
    ModuleBase::timer::tick("PW_Basis","get_GVectors");

    this->nrxx = this->ncxyz;
    this->ngmc=this->ngmc_g;

    //************************************************************
    // g  : Store the G vectors in 1d array (Cartian coordinate)
    // ig : Store the G vectors in 1d array (Direct coordinate)
    // gg : store the |G|^2 of the 1d array
    //************************************************************

    //----------------------------------------------------------
    // EXPLAIN : if not parallel case, we use pointer
    // g and g_global are pointers to the same array
    //----------------------------------------------------------
    this->gcar = this->gcar_global;
    this->gdirect = this->gdirect_global;
    this->gg = this->gg_global;

    // (2) calculate ig2fftc
    assert(ngmc>0);
    delete[] ig2fftc;
    delete[] ig1;
    delete[] ig2;
    delete[] ig3;
    this->ig2fftc = new int[ngmc];
    this->ig1 = new int[ngmc];
    this->ig2 = new int[ngmc];
    this->ig3 = new int[ngmc];
	ModuleBase::GlobalFunc::ZEROS(ig2fftc, ngmc);
	ModuleBase::GlobalFunc::ZEROS(ig1, ngmc);
	ModuleBase::GlobalFunc::ZEROS(ig2, ngmc);
	ModuleBase::GlobalFunc::ZEROS(ig3, ngmc);
    PW_complement::get_ig2fftc(ngmc, ncx, ncy, ncz, gdirect, ig1, ig2, ig3, ig2fftc);

    // (3) calculate ngmw: number of plane wave to describe wave functions.
    PW_complement::get_ngmw(ngmc, ggwfc2, gg_global, ngmw);

    assert(ngmw>0);
    delete[] ig2fftw;
    ig2fftw = new int[ngmw];
	ModuleBase::GlobalFunc::ZEROS(ig2fftw, ngmw);
    ModuleBase::Memory::record("PW_complement","ig2fftw",ngmw,"int");
    ModuleBase::GlobalFunc::ZEROS(ig2fftw, ngmw);

    PW_complement::get_ig2fftw(ngmw, nx, ny, nz, gdirect, ig2fftw);

    ModuleBase::timer::tick("PW_Basis","get_GVectors");
    return;
}//end get_GVectors;

void PW_Basis::get_nggm(const int ngmc_local)
{
    ModuleBase::TITLE("PW_Basis","get_nggm");
    ModuleBase::timer::tick("PW_Basis","get_nggm");

//	GlobalV::ofs_running << " calculate the norm of G vectors." << std::endl;
    //*********************************************
    // Group the g vectors that have the same norm
    //*********************************************
	assert(ngmc_local>0);
    double *tmp = new double[ngmc_local];
	ModuleBase::GlobalFunc::ZEROS(tmp, ngmc_local);


    delete[] ig2ngg;
    this->ig2ngg = new int[ngmc_local];
	ModuleBase::GlobalFunc::ZEROS(ig2ngg, ngmc_local);

    tmp[0] = this->gg[0];
    ig2ngg[0] = 0;
    int ng = 0;

    for (int ig = 1; ig < ngmc_local; ig++)
    {
        if (abs(this->gg[ig] - tmp[ng]) > 1.0e-8)
        {
            ng++;
            tmp[ng] = this->gg[ig];
        }
        this->ig2ngg[ig] = ng;
    }
    ng++;
	

    //********************************
    // number of different |G| values
    // *******************************
    this->nggm = ng;

	ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running,"number of |g|",nggm);

    delete[] ggs;
    this->ggs = new double[this->nggm];
	ModuleBase::GlobalFunc::ZEROS(ggs, this->nggm);

	// mohan update 2011-06-12
	for(int ig=0; ig<nggm; ig++)
	{
		this->ggs[ig] = tmp[ig];
	}

	ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running,"max |g|",ggs[nggm-1]);
	ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running,"min |g|",ggs[0]);

    if (this->ggs[0]>1.0e-4) // mohan modified 2009-07-07
    {
        gstart = 0;
    }
    else
    {
        gstart = 1;
    }

	if(GlobalV::test_pw)ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running,"gstart",gstart);

    delete[] tmp;
    ModuleBase::timer::tick("PW_Basis","get_nggm");
    return;
}

#include "../module_base/constants.h"
//  Calculate structure factor
void PW_Basis::setup_structure_factor(void)			// Peize Lin optimize and add OpenMP 2021.04.01
{
    ModuleBase::TITLE("PW_Basis","setup_structure_factor");
    ModuleBase::timer::tick("PW_Basis","setup_struc_factor");
    const std::complex<double> ci_tpi = ModuleBase::NEG_IMAG_UNIT * ModuleBase::TWO_PI;

    this->strucFac.create(Ucell->ntype, this->ngmc);
    ModuleBase::Memory::record("PW_Basis","struc_fac", Ucell->ntype*this->ngmc,"complexmatrix");

//	std::string outstr;
//	outstr = GlobalV::global_out_dir + "strucFac.dat"; 
//	std::ofstream ofs( outstr.c_str() ) ;
    bool usebspline;
    if(nbspline > 0)   usebspline = true;
    else    usebspline = false;
    
    if(usebspline)
    {
        nbspline = int((nbspline+1)/2)*2; // nbspline must be a positive even number.
        this->bspline_sf(nbspline);
    }
    else
    {
        for (int it=0; it<Ucell->ntype; it++)
        {
	    	const int na = Ucell->atoms[it].na;
	    	const ModuleBase::Vector3<double> * const tau = Ucell->atoms[it].tau;
#ifdef _OPENMP
		    #pragma omp parallel for schedule(static)
#endif
            for (int ig=0; ig<this->ngmc; ig++)
            {
		    	const ModuleBase::Vector3<double> gcar_ig = gcar[ig];
                std::complex<double> sum_phase = ModuleBase::ZERO;
                for (int ia=0; ia<na; ia++)
                {
                    // e^{-i G*tau}
                    sum_phase += exp( ci_tpi * (gcar_ig * tau[ia]) );
                }
                this->strucFac(it,ig) = sum_phase;
            }
        }
    }

//	ofs.close();

    int i,j; //ng;
    this->eigts1.create(Ucell->nat, 2*this->ncx + 1);
    this->eigts2.create(Ucell->nat, 2*this->ncy + 1);
    this->eigts3.create(Ucell->nat, 2*this->ncz + 1);

    ModuleBase::Memory::record("PW_Basis","eigts1",Ucell->nat*2*this->ncx + 1,"complexmatrix");
    ModuleBase::Memory::record("PW_Basis","eigts2",Ucell->nat*2*this->ncy + 1,"complexmatrix");
    ModuleBase::Memory::record("PW_Basis","eigts3",Ucell->nat*2*this->ncz + 1,"complexmatrix");

    ModuleBase::Vector3<double> gtau;
    int inat = 0;
    for (i = 0; i < Ucell->ntype; i++)
    {
        if (GlobalV::test_pw > 1)
        {
            ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running,"eigts",i);
        }
        for (j = 0; j < Ucell->atoms[i].na;j++)
        {
            gtau = Ucell->G * Ucell->atoms[i].tau[j];  //HLX: fixed on 10/13/2006
            for (int n1 = -ncx; n1 <= ncx;n1++)
            {
                double arg = n1 * gtau.x;
                this->eigts1(inat, n1 + ncx) = exp( ci_tpi*arg  );
            }
            for (int n2 = -ncy; n2 <= ncy;n2++)
            {
                double arg = n2 * gtau.y;
                this->eigts2(inat, n2 + ncy) = exp( ci_tpi*arg );
            }
            for (int n3 = -ncz; n3 <= ncz;n3++)
            {
                double arg = n3 * gtau.z;
                this->eigts3(inat, n3 + ncz) = exp( ci_tpi*arg );
            }
            inat++;
        }
    }
    ModuleBase::timer::tick("PW_Basis","setup_struc_factor"); 
    return;
}

#ifdef __MPI
void PW_Basis::columns_and_pw_distribution_2(void)
{
    ModuleBase::TITLE("PW_Basis","columns_and_pw_distribution_2");

    // time count the number of sticks in charge grid.
    int time=0;
    //int sum_ngmc_g=0;
    // circle of charge sticks.
    while (time<FFT_chg.nst)
    {
        time++;
        int ip1=0;
        bool find=false;
        int max_npw=0;
        int max_i = 0;
        int max_j = 0;
        int pw_wf=0;
        
		// which column now has max plane wave, select it.
		FFT_chg.max_pw_column(max_npw,max_i,max_j);

        for (int j=0;j<FFT_wfc.nst;j++)
        {
            if (FFT_wfc.st_i[j] == max_i && FFT_wfc.st_j[j] == max_j )
            {
                //------------------------------------------------------
                // if this charge column is also one of the wave function column,
                // Then we distribute this column according to each process's
                // planewave number or column numbers.
                //------------------------------------------------------
                if (FFT_wfc.st_n[j]>0)
                {
                    find=true;
                    // number of planewaves for wave function grid.
                    pw_wf=FFT_wfc.st_n[j];
                    for (int ip2=0;ip2<GlobalV::NPROC_IN_POOL;ip2++)
                    {
						const int ngrid = this->nx*this->ny*FFT_wfc.npps[ip2];
						const int non_zero_grid = FFT_wfc.nst_per[ip2]*this->nz;
						
						const int npw1 = FFT_wfc.npw_per[ip1];
						const int npw2 = FFT_wfc.npw_per[ip2];

						const int nst1 = FFT_wfc.nst_per[ip1];
						const int nst2 = FFT_wfc.nst_per[ip2];

                        //------------------------------------------------------
                        // compare planewave numbers in each processor
                        // npw_per : number of planewaves in processor 'ip'
                        //------------------------------------------------------
                        if (npw2 < npw1)
                        {
							if (non_zero_grid + nz <= ngrid) //qianrui fix a bug 2021-5-20 to make sure non_zero_grid < ngrid after distributing pw
							{
								// ip1 save the process number which has smallest number of plane wave
								// in this pool.
								// P.S. There's one more condition, that the index should not exceed
								// the number of FFT grid in each processor.
								ip1=ip2;
							}
						}
                        //----------------------------------------------------
                        // if two processor has the same planewave number,
                        // then we compare column number in each processor
                        //----------------------------------------------------
                        else if ( npw2 == npw1)
                        {
                            if (nst2 < nst1)
                            {
                                if (non_zero_grid + nz <= ngrid) //qianrui add nz
                                {
                                    ip1=ip2;
                                }
                            }
                        }
                    }//end ip2
                }//end FFT_wfc.st_n
            }//end FFT_wfc.st__i ; FFT_wfc.st__j
        }//end j
        // if this stick only belong to charge grid: find = true.
        if (!find)
        {
            // search from '0' processor.
            int pw_tmp=FFT_chg.npw_per[0];
            //std::cout<<"pw_tmp="<<pw_tmp<<std::endl;
            for (int ip2=0;ip2<GlobalV::NPROC_IN_POOL;ip2++)
            {

				const int ngrid = this->nx*this->ny*FFT_chg.npps[ip2];
				const int non_zero_grid = FFT_chg.nst_per[ip2]*this->nz;

				//const int npw1 = FFT_chg.npw_per[ip1];
				const int npw2 = FFT_chg.npw_per[ip2];

				//const int nst1 = FFT_chg.nst_per[ip1];
				//const int nst2 = FFT_chg.nst_per[ip2];

                if (npw2 < pw_tmp)
                {
                    pw_tmp = npw2;
                    if (non_zero_grid + nz <= ngrid) //qianrui add nz
                    {
                        // ip1 is the index of processor which
                        // has smallest number of plane wave.
                        ip1=ip2;
                    }
                }
            }
        }

        // update the information:
        FFT_chg.nst_per[ip1]++;
        FFT_chg.npw_per[ip1]+=max_npw;
        FFT_chg.index_ip[max_j+FFT_chg.n2*max_i] = ip1;//mohan 2009-11-09

        if (find)
        {
            FFT_wfc.nst_per[ip1]++;
            FFT_wfc.npw_per[ip1]+=pw_wf;
            //std::cout<<std::setw(12)<<"time"<<std::setw(12)<<"Cols"<<std::setw(12)<<"PWs"<<std::endl;
            //for(int i=0;i<GlobalV::NPROC_IN_POOL;i++)
            //{
            //	std::cout<<std::setw(12)<<time<<std::setw(12)<<FFT_wfc.nst_per[i]<<std::setw(12)<<FFT_wfc.npw_per[i]<<std::endl;
            //}
            
			//----------------------------------------------------------
			// this may be the reason why x and y is not correct!!!!!!!!
			// mohan 2011-04-21
			//----------------------------------------------------------
			FFT_wfc.index_ip[max_j+FFT_wfc.n2*max_i] = ip1;// mohan 2009-11-09
        }
    }//end i
    int sum_pw=0;
    int i=0;

	//GlobalV::ofs_running << " distribution of plane waves in pool " << GlobalV::MY_POOL << std::endl;
	//GlobalV::ofs_running << " number of processors in this pool is " << GlobalV::NPROC_IN_POOL << std::endl;
    if (1)
    {
		GlobalV::ofs_running << "\n PARALLEL PW FOR CHARGE/POTENTIAL" << std::endl;
        GlobalV::ofs_running << " "
        << std::setw(8)  << "PROC"
        << std::setw(15) << "COLUMNS(POT)"
        << std::setw(15) << "PW" << std::endl;
// charge
        for (i = 0;i < GlobalV::NPROC_IN_POOL;i++)
        {
            GlobalV::ofs_running << " "
            << std::setw(8)  << i+1
            << std::setw(15) << FFT_chg.nst_per[i]
            << std::setw(15) << FFT_chg.npw_per[i] << std::endl;
            sum_pw += FFT_chg.npw_per[i];
        }
        GlobalV::ofs_running << " --------------- sum -------------------" << std::endl;
        GlobalV::ofs_running << " "
        << std::setw(8)  << GlobalV::NPROC_IN_POOL
        << std::setw(15) << FFT_chg.nst
        << std::setw(15) << sum_pw << std::endl;
// wavefunction


		GlobalV::ofs_running << "\n PARALLEL PW FOR WAVE FUNCTIONS" << std::endl;
        sum_pw=0;
        GlobalV::ofs_running << " "
        << std::setw(8)  << "PROC"
        << std::setw(15) << "COLUMNS(W)"
        << std::setw(15) << "PW" << std::endl;

        for (i = 0;i < GlobalV::NPROC_IN_POOL;i++)
        {
            GlobalV::ofs_running << " "
            << std::setw(8)  << i+1
            << std::setw(15) << FFT_wfc.nst_per[i]
            << std::setw(15) << FFT_wfc.npw_per[i] << std::endl;
            sum_pw += FFT_wfc.npw_per[i];
        }
        GlobalV::ofs_running << " --------------- sum -------------------" << std::endl;
        GlobalV::ofs_running << " "
        << std::setw(8)  << GlobalV::NPROC_IN_POOL
        << std::setw(15) << FFT_wfc.nst
        << std::setw(15) << sum_pw << std::endl;
    }

    for (int ip=0;ip<GlobalV::NPROC_IN_POOL;ip++)
    {
        //	GlobalV::ofs_running<<"\n FFT_chg.npps = "<<FFT_chg.npps[ip]<<std::endl;
        //	GlobalV::ofs_running<<"\n FFT_chg.nst_per = "<<FFT_chg.nst_per[ip]<<std::endl;
        int ngrid = this->ncx*this->ncy*FFT_chg.npps[ip];
        int non_zero_grid = FFT_chg.nst_per[ip]*this->ncz;
        if (non_zero_grid > ngrid)
        {
            GlobalV::ofs_running<<" too many sticks for cpu = "<<ip<<std::endl;
            GlobalV::ofs_running<<" ngrid is = "<< ngrid << std::endl;
            GlobalV::ofs_running<<" In fact , non_zero_grid = "<< non_zero_grid << std::endl;
            ModuleBase::WARNING_QUIT("PW_Basis::columns_and_pw_distribution_2","conflict about pw distribution.");
        }
    }


	// check. mohan add 2011-04-25
	int no_pw = 0;
	for(int i=0; i<GlobalV::NPROC_IN_POOL; i++)
	{
		if(FFT_wfc.npw_per[i]==0)
		{
			++no_pw;
			break;
		}
	}
	Parallel_Reduce::reduce_int_all(no_pw);

	if(no_pw>0)
	{
		ModuleBase::WARNING_QUIT("distribution of pw","some processor has no plane waves!");
	}


    return;
}
#endif

//LiuXh add a new function here,
//20180515
void PW_Basis::update_gvectors(std::ofstream &runlog, const UnitCell &Ucell_in)
{
    ModuleBase::TITLE("PW_Basis","update_gvectors");
    ModuleBase::timer::tick("PW_Basis","update_gvectors");

#ifdef __MPI
    bool cutgg_flag = true;

#else
    bool cutgg_flag = false;
#endif

    if (cutgg_flag)
    {
#ifdef __MPI
        ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running,"number of total plane waves",ngmc_g);

        double cutgg_pieces = 10;

        const double cutgg_delta = ggchg / std::pow( (double)cutgg_pieces, 2.0/3.0 );

        //int cutgg_num_start = 0;
        double ggchg_start = 0.0;
        double ggchg_end = 0.0;

        int ggchg_time = 1;
        do
        {
            ggchg_end = cutgg_delta * std::pow( (double)ggchg_time, 2.0/3.0 ) ;
            if ( abs(ggchg_end - ggchg)<1.0e-8 ) ggchg_end = ggchg;

            delete[] gg_global0;
            delete[] gg_global;
            delete[] gdirect_global;
            delete[] gcar_global;

            int ii = ggchg_time-1;
            int cutgg_num_now = cutgg_num_table[ii];
            int cutgg_num_now2 = cutgg_num_table[ii];
            gg_global0 = new double[cutgg_num_now2];
            gg_global = new double[cutgg_num_now2];
            gdirect_global = new ModuleBase::Vector3<double>[cutgg_num_now2];
            gcar_global = new ModuleBase::Vector3<double>[cutgg_num_now2];

            PW_complement::get_total_pw_after_vc(gg_global0, gg_global, gdirect_global, ggchg_start, ggchg_end,
                    ncx, ncy, ncz, Ucell->GGT, Ucell->GGT0, ngmc_g);

            PW_complement::setup_GVectors(Ucell->G, cutgg_num_now, gg_global, gdirect_global, gcar_global);

            FFT_chg.fft_map_after_vc(this->ig2fftc, this->ngmc, cutgg_num_now, ggchg_time);
            FFT_wfc.fft_map_after_vc(this->ig2fftw, this->ngmw, cutgg_num_now, ggchg_time);

            ggchg_start = ggchg_end;
            ++ggchg_time;
        } while ( ggchg_end < ggchg );
#endif
    }



#ifdef __MPI
    this->get_MPI_GVectors();
#endif

    this->get_nggm(this->ngmc);

    this->setup_structure_factor();

    ModuleBase::timer::tick("PW_Basis","update_gvectors");


    return;
}
