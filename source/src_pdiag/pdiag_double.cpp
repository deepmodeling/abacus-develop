/*test generalized Stardand double precision symmetric eigenproblem*/
#include "pdiag_double.h"
#include "../module_base/lapack_connector.h"
#include "../src_pw/occupy.h"
#include "../src_pw/global.h"
#include "../src_lcao/local_orbital_charge.h"
#include "../src_io/wf_local.h"
#include "../module_base/memory.h"
#include "../module_base/timer.h"


#ifdef __MPI
extern "C"
{
    #include "../module_base/blacs_connector.h"
    #include "my_elpa.h"
	#include "../module_base/scalapack_connector.h"
}
#include "pdgseps.h"
#include "pzgseps.h"
#include "../module_base/lapack_connector.h"
#endif

#include "../src_external/src_test/test_function.h"


#ifdef __MPI
inline int set_elpahandle(elpa_t &handle, const int *desc,const int local_nrows,const int local_ncols, const int nbands)
{
  int error;
  int nprows, npcols, myprow, mypcol;
  Cblacs_gridinfo(desc[1], &nprows, &npcols, &myprow, &mypcol);
  elpa_init(20210430);
  handle = elpa_allocate(&error);
  elpa_set_integer(handle, "na", desc[2], &error);
  elpa_set_integer(handle, "nev", nbands, &error);

  elpa_set_integer(handle, "local_nrows", local_nrows, &error);

  elpa_set_integer(handle, "local_ncols", local_ncols, &error);

  elpa_set_integer(handle, "nblk", desc[4], &error);

  elpa_set_integer(handle, "mpi_comm_parent", MPI_Comm_c2f(MPI_COMM_WORLD), &error);

  elpa_set_integer(handle, "process_row", myprow, &error);

  elpa_set_integer(handle, "process_col", mypcol, &error);

  elpa_set_integer(handle, "blacs_context", desc[1], &error);

  elpa_set_integer(handle, "cannon_for_generalized", 0, &error);
   /* Setup */
  elpa_setup(handle);   /* Set tunables */
  return 0;
}
#endif


inline bool ifElpaHandle(const bool& newIteration, const bool& ifNSCF)
{
    int doHandle = false;
	if(newIteration) doHandle = true;
	if(ifNSCF) doHandle = true;
	return doHandle;
}

int Pdiag_Double::out_mat_hs = 0;
int Pdiag_Double::out_mat_hsR = 0;
int Pdiag_Double::out_wfc_lcao = 0;

Pdiag_Double::Pdiag_Double() {}

Pdiag_Double::~Pdiag_Double(){}


void Pdiag_Double::diago_double_begin(
	const int &ik, // k-point index
	Local_Orbital_wfc &lowf,
	double* h_mat, // hamiltonian matrix
    double* s_mat, // overlap matrix
    double* Stmp, 	// because the output Stmp will be different from Sloc2, so we need to copy that.
    double* ekb) // eigenvalues for each k-point and band
{
    const Parallel_Orbitals* pv = this->ParaV = lowf.ParaV;
    
#ifdef TEST_DIAG
	{
		static int istep = 0;
		auto print_matrix_C = [&](const std::string &file_name, double*m)
		{
			std::ofstream ofs(file_name+"-C_"+ModuleBase::GlobalFunc::TO_STRING(istep)+"_"+ModuleBase::GlobalFunc::TO_STRING(GlobalV::MY_RANK));
			for(int ic=0; ic<pv->ncol; ++ic)
			{
				for(int ir=0; ir<pv->nrow; ++ir)
				{
					const int index=ic*pv->nrow+ir;
					if(abs(m[index])>1E-10)
						ofs<<m[index]<<"\t";
					else
						ofs<<0<<"\t";
				}
				ofs<<std::endl;
			}
		};
		auto print_matrix_F = [&](const std::string &file_name, double*m)
		{
			std::ofstream ofs(file_name+"-F_"+ModuleBase::GlobalFunc::TO_STRING(istep)+"_"+ModuleBase::GlobalFunc::TO_STRING(GlobalV::MY_RANK));
			for(int ir=0; ir<pv->nrow; ++ir)
			{
				for(int ic=0; ic<pv->ncol; ++ic)
				{
					const int index=ic*pv->nrow+ir;
					if(abs(m[index])>1E-10)
						ofs<<m[index]<<"\t";
					else
						ofs<<0<<"\t";
				}
				ofs<<std::endl;
			}
		};
		print_matrix_F("H_gamma", h_mat);
		print_matrix_F("S_gamma", s_mat);
		print_matrix_C("H_gamma", h_mat);
		print_matrix_C("S_gamma", s_mat);
		++istep;
	}
	#endif

#ifdef __MPI
	ModuleBase::TITLE("Pdiag_Double","diago_begin");
	assert(pv->loc_size > 0);
	assert(GlobalV::NLOCAL > 0);

	char uplo='U';
	const int inc=1;

    ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running,"start solver, ks_solver",GlobalV::KS_SOLVER);
    if(GlobalV::KS_SOLVER=="hpseps")
    {
        // parallel diagonalize the
        // H | psi > = S | psi >
        // problem.
        int loc_pos;
        
        int  myid;
        MPI_Comm_rank(pv->comm_2D, &myid);

        double* eigen = new double[GlobalV::NLOCAL];
        ModuleBase::GlobalFunc::ZEROS(eigen, GlobalV::NLOCAL);

        double* Z = new double[pv->loc_size * GlobalV::NLOCAL];
        ModuleBase::GlobalFunc::ZEROS(Z, pv->loc_size * GlobalV::NLOCAL);

        ModuleBase::Memory::record("Pdiag_Double","Z", pv->loc_size * GlobalV::NLOCAL,"double");
        ModuleBase::timer::tick("Diago_LCAO_Matrix","pdgseps");
        BlasConnector::copy(pv->nloc, s_mat, inc, Stmp, inc);
        int loc_size = pv->loc_size;
        pdgseps(pv->comm_2D, GlobalV::NLOCAL, pv->nb, h_mat, Stmp, Z, eigen, pv->MatrixInfo, uplo, loc_size, loc_pos);
        ModuleBase::timer::tick("Diago_LCAO_Matrix","pdgseps");

        if(myid <= pv->lastband_in_proc)
        {
            for(int i=0; i<pv->loc_sizes[myid]; i++)
            {
                for(int n=0; n<GlobalV::NLOCAL; n++)
                {
                    pv->Z_LOC[ik][n*pv->loc_sizes[myid] + i] = Z[n*pv->loc_sizes[myid] + i];
                }
            }
        }

        // the eigenvalues.
        //xiaohui modify 2014-06-15, move to the top
        BlasConnector::copy(GlobalV::NBANDS, eigen, inc, ekb, inc);
        delete[] eigen;
        //=====================================
        // gather the eigenvectors and
        // distribute them to each processor
        // Z is delete in gath_eig
        //=====================================

        //xiaohui modify 2014-06-18
        ModuleBase::timer::tick("Diago_LCAO_Matrix","gath_eig");
        this->gath_eig(DIAG_HPSEPS_WORLD, GlobalV::NLOCAL, Z);
        ModuleBase::timer::tick("Diago_LCAO_Matrix","gath_eig");
        //delete[] Z; //LiuXh 20171109
	}// HPSEPS method
    else if(GlobalV::KS_SOLVER=="genelpa")
    {
        double *eigen = new double[GlobalV::NLOCAL];
        ModuleBase::GlobalFunc::ZEROS(eigen, GlobalV::NLOCAL);

        long maxnloc; // maximum number of elements in local matrix
        MPI_Reduce(&pv->nloc_wfc, &maxnloc, 1, MPI_LONG, MPI_MAX, 0, pv->comm_2D);
        MPI_Bcast(&maxnloc, 1, MPI_LONG, 0, pv->comm_2D);
        lowf.wfc_gamma[ik].create(pv->ncol, pv->nrow);			// Fortran order
        
        static elpa_t handle;
        static bool has_set_elpa_handle = false;
        if(! has_set_elpa_handle)
        {
            set_elpahandle(handle, pv->desc, pv->nrow, pv->ncol, GlobalV::NBANDS);
            has_set_elpa_handle = true;
        }

        int is_already_decomposed;
        if(ifElpaHandle(GlobalC::CHR.get_new_e_iteration(), (GlobalV::CALCULATION=="nscf")))
        {
            ModuleBase::timer::tick("Diago_LCAO_Matrix","decompose_S");
            BlasConnector::copy(pv->nloc, s_mat, inc, Stmp, inc);
            is_already_decomposed=0;
            ModuleBase::timer::tick("Diago_LCAO_Matrix","decompose_S");
        }
        else
        {
            is_already_decomposed=1;
        }

        ModuleBase::timer::tick("Diago_LCAO_Matrix","elpa_solve");
        int elpa_error;
        elpa_generalized_eigenvectors_d(handle, h_mat, Stmp, eigen, lowf.wfc_gamma[ik].c, is_already_decomposed, &elpa_error);
        ModuleBase::timer::tick("Diago_LCAO_Matrix", "elpa_solve");

    	ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running,"K-S equation was solved by genelpa2");
        BlasConnector::copy(GlobalV::NBANDS, eigen, inc, ekb, inc);
        delete[] eigen;
	    ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running,"eigenvalues were copied to ekb");

        double** wfc_grid = nullptr;    //output but not do "2d-to-grid" conversion
        lowf.wfc_2d_to_grid(this->out_wfc_lcao, lowf.wfc_gamma[ik].c, wfc_grid);
        
    } // GenELPA method
	else if(GlobalV::KS_SOLVER=="lapack_gv")
	{
		lowf.wfc_gamma[ik].create(pv->ncol_bands, pv->nrow, false);
		memcpy( lowf.wfc_gamma[ik].c, h_mat, sizeof(double)*pv->ncol*pv->nrow );
		ModuleBase::matrix s_tmp(pv->ncol, pv->nrow, false);
		memcpy( s_tmp.c, s_mat, sizeof(double)*pv->ncol*pv->nrow );
		std::vector<double> ekb_tmp(GlobalV::NLOCAL,0);

		const char jobz='V', uplo='U';
		const int itype=1;
		int lwork=-1, info=0;
		std::vector<double> work(1,0);
		dsygv_(&itype, &jobz, &uplo, &GlobalV::NLOCAL, lowf.wfc_gamma[ik].c, &GlobalV::NLOCAL,
			s_tmp.c, &GlobalV::NLOCAL, ekb_tmp.data(), work.data(), &lwork, &info);

		if(info)
		{
			throw std::runtime_error("info="+ModuleBase::GlobalFunc::TO_STRING(info)+". "+ModuleBase::GlobalFunc::TO_STRING(__FILE__)+" line "+ModuleBase::GlobalFunc::TO_STRING(__LINE__));
		}

		lwork = work[0];
		work.resize(lwork);

		dsygv_(&itype, &jobz, &uplo, &GlobalV::NLOCAL, lowf.wfc_gamma[ik].c, &GlobalV::NLOCAL,
			s_tmp.c, &GlobalV::NLOCAL, ekb_tmp.data(), work.data(), &lwork, &info);

		if(info)
		{
			throw std::runtime_error("info="+ModuleBase::GlobalFunc::TO_STRING(info)+". "+ModuleBase::GlobalFunc::TO_STRING(__FILE__)+" line "+ModuleBase::GlobalFunc::TO_STRING(__LINE__));
		}
		memcpy( ekb, ekb_tmp.data(), sizeof(double)*GlobalV::NBANDS );

	}
	else if(GlobalV::KS_SOLVER=="lapack_gvx")
	{
		ModuleBase::matrix h_tmp(pv->ncol, pv->nrow, false);
		memcpy( h_tmp.c, h_mat, sizeof(double)*pv->ncol*pv->nrow );
		ModuleBase::matrix s_tmp(pv->ncol, pv->nrow, false);
		memcpy( s_tmp.c, s_mat, sizeof(double)*pv->ncol*pv->nrow );
		lowf.wfc_gamma[ik].create(pv->ncol, pv->nrow, false);

		const char jobz='V', range='I', uplo='U';
		const int itype=1, il=1, iu=GlobalV::NBANDS;
		int M=0, lwork=-1, info=0;
		const double abstol=0;
		std::vector<double> work(1,0);
		std::vector<int> iwork(5*GlobalV::NLOCAL,0);
		std::vector<int> ifail(GlobalV::NLOCAL,0);

		dsygvx_(&itype, &jobz, &range, &uplo,
			&GlobalV::NLOCAL, h_tmp.c, &GlobalV::NLOCAL, s_tmp.c, &GlobalV::NLOCAL, NULL, NULL, &il, &iu, &abstol,
			&M, ekb, lowf.wfc_gamma[ik].c, &GlobalV::NLOCAL, work.data(), &lwork, iwork.data(), ifail.data(), &info);

		if(info)
		{
			throw std::runtime_error("info="+ModuleBase::GlobalFunc::TO_STRING(info)+". "+ModuleBase::GlobalFunc::TO_STRING(__FILE__)+" line "+ModuleBase::GlobalFunc::TO_STRING(__LINE__));
		}

		lwork = work[0];
		work.resize(lwork);
		dsygvx_(&itype, &jobz, &range, &uplo,
			&GlobalV::NLOCAL, h_tmp.c, &GlobalV::NLOCAL, s_tmp.c, &GlobalV::NLOCAL, NULL, NULL, &il, &iu, &abstol,
			&M, ekb, lowf.wfc_gamma[ik].c, &GlobalV::NLOCAL, work.data(), &lwork, iwork.data(), ifail.data(), &info);

		if(info)
		{
			throw std::runtime_error("info="+ModuleBase::GlobalFunc::TO_STRING(info)+". "+ModuleBase::GlobalFunc::TO_STRING(__FILE__)+" line "+ModuleBase::GlobalFunc::TO_STRING(__LINE__));
		}
		if(M!=GlobalV::NBANDS)
		{
			throw std::runtime_error("M="+ModuleBase::GlobalFunc::TO_STRING(M)+". GlobalV::NBANDS="+ModuleBase::GlobalFunc::TO_STRING(GlobalV::NBANDS)+". "+ModuleBase::GlobalFunc::TO_STRING(__FILE__)+" line "+ModuleBase::GlobalFunc::TO_STRING(__LINE__));
		}

	}
	else if(GlobalV::KS_SOLVER=="scalapack_gvx")
	{
		diag_scalapack_gvx.pdsygvx_diag(pv->desc, pv->ncol, pv->nrow, h_mat, s_mat, ekb, lowf.wfc_gamma[ik]);		// Peize Lin add 2021.11.02
	}
    //delete[] Stmp; //LiuXh 20171109
#endif

#ifdef TEST_DIAG
	{
		static int istep = 0;
		{
			std::ofstream ofs("ekb_"+ModuleBase::GlobalFunc::TO_STRING(istep)+"_"+ModuleBase::GlobalFunc::TO_STRING(GlobalV::MY_RANK));
			for(int ib=0; ib<GlobalV::NBANDS; ++ib)
			{
				ofs<<ekb[ib]<<std::endl;
			}
		}
		{
			std::ofstream ofs("wfc-C_"+ModuleBase::GlobalFunc::TO_STRING(istep)+"_"+ModuleBase::GlobalFunc::TO_STRING(GlobalV::MY_RANK));
			ofs<<lowf.wfc_gamma[ik]<<std::endl;
		}
		{
			std::ofstream ofs("wfc-F_"+ModuleBase::GlobalFunc::TO_STRING(istep)+"_"+ModuleBase::GlobalFunc::TO_STRING(GlobalV::MY_RANK));
			ofs<<transpose(lowf.wfc_gamma[ik])<<std::endl;
		}
		++istep;
	}
#endif

	return;
}

void Pdiag_Double::diago_complex_begin(
	const int &ik,
    Local_Orbital_wfc &lowf,
	std::complex<double>* ch_mat,
    std::complex<double>* cs_mat,
    std::complex<double>* Stmp, 	// because the output Stmp will be different from Sloc2, so we need to copy that.
    double* ekb)
{
    const Parallel_Orbitals* pv = this->ParaV = lowf.ParaV;
    
#ifdef TEST_DIAG
   	{
		static int istep = 0;
		auto print_matrix_C = [&](const std::string &file_name, std::complex<double>*m)
		{
			std::ofstream ofs(file_name+"-C_"+ModuleBase::GlobalFunc::TO_STRING(istep)+"_"+ModuleBase::GlobalFunc::TO_STRING(GlobalV::MY_RANK));
			for(int ic=0; ic<pv->ncol; ++ic)
			{
				for(int ir=0; ir<pv->nrow; ++ir)
				{
					const int index=ic*pv->nrow+ir;
					if(std::norm(m[index])>1E-10)
                    {
                        if(std::imag(m[index])>1E-10)
						{
                            ofs<<m[index]<<"\t";
						}
                        else
						{
                            ofs<<std::real(m[index])<<"\t";
						}
                    }
					else
					{
						ofs<<0<<"\t";
					}
				}
				ofs<<std::endl;
			}
		};
		auto print_matrix_F = [&](const std::string &file_name, std::complex<double>*m)
		{
			std::ofstream ofs(file_name+"-F_"+ModuleBase::GlobalFunc::TO_STRING(istep)+"_"+ModuleBase::GlobalFunc::TO_STRING(GlobalV::MY_RANK));
			for(int ir=0; ir<pv->nrow; ++ir)
			{
				for(int ic=0; ic<pv->ncol; ++ic)
				{
					const int index=ic*pv->nrow+ir;
					if(std::norm(m[index])>1E-10)
                    {
                        if(std::imag(m[index])>1E-10)
						{
                            ofs<<m[index]<<"\t";
						}
                        else
						{
                            ofs<<std::real(m[index])<<"\t";
						}
                    }
					else
					{
						ofs<<0<<"\t";
					}
				}
				ofs<<std::endl;
			}
		};
		print_matrix_F("H_gamma", ch_mat);
		print_matrix_F("S_gamma", cs_mat);
		print_matrix_C("H_gamma", ch_mat);
		print_matrix_C("S_gamma", cs_mat);
		++istep;
	}
    #endif

#ifdef __MPI
	ModuleBase::TITLE("Pdiag_Double","diago_complex_begin");

	char uplo='U';
	const int inc=1;

	if(GlobalV::KS_SOLVER=="hpseps")
    {
        // parallel diagonalize the
        // H | psi > = S | psi >
        // problem.
        int loc_pos;
        
        double* eigen = new double[GlobalV::NLOCAL];
        ModuleBase::GlobalFunc::ZEROS(eigen, GlobalV::NLOCAL);

        assert(pv->loc_size > 0);
        std::complex<double>* Z = new std::complex<double>[pv->loc_size * GlobalV::NLOCAL];
        ModuleBase::GlobalFunc::ZEROS(Z, pv->loc_size * GlobalV::NLOCAL);

        ModuleBase::Memory::record("Pdiag_Double","Z",pv->loc_size * GlobalV::NLOCAL,"cdouble");
		int nbands_tmp = GlobalV::NBANDS;
        ModuleBase::timer::tick("Diago_LCAO_Matrix","pzgseps");
        BlasConnector::copy(pv->nloc, cs_mat, inc, Stmp, inc);
        int loc_size = pv->loc_size;
        pzgseps(pv->comm_2D, GlobalV::NLOCAL, pv->nb, nbands_tmp, ch_mat, Stmp, Z, eigen, pv->MatrixInfo, uplo, loc_size, loc_pos);
        ModuleBase::timer::tick("Diago_LCAO_Matrix","pzgseps");
        // the eigenvalues.
        BlasConnector::copy(GlobalV::NBANDS, eigen, inc, ekb, inc);
        delete[] eigen;

        // Z is delete in gath_eig
        ModuleBase::timer::tick("Diago_LCAO_Matrix","gath_eig_complex");
        this->gath_eig_complex(DIAG_HPSEPS_WORLD, GlobalV::NLOCAL, lowf.wfc_k_grid[ik], Z, ik);
        ModuleBase::timer::tick("Diago_LCAO_Matrix","gath_eig_complex");
        //delete[] Z; //LiuXh 20180329, fix bug of 'double free()'
        //this->gath_full_eig_complex(DIAG_WORLD, GlobalV::NLOCAL, c, Z);
	} // HPSEPS method
    else if(GlobalV::KS_SOLVER=="genelpa")
    {
        double *eigen = new double[GlobalV::NLOCAL];
        ModuleBase::GlobalFunc::ZEROS(eigen, GlobalV::NLOCAL);
        lowf.wfc_k[ik].create(pv->ncol_bands, pv->nrow);            // Fortran order

        static elpa_t handle;
        static bool has_set_elpa_handle = false;
        if(! has_set_elpa_handle)
        {
            set_elpahandle(handle, pv->desc, pv->nrow, pv->ncol, GlobalV::NBANDS);
            has_set_elpa_handle = true;
        }

        BlasConnector::copy(pv->nloc, cs_mat, inc, Stmp, inc);

        ModuleBase::timer::tick("Diago_LCAO_Matrix","elpa_solve");
        int elpa_derror;
        elpa_generalized_eigenvectors_dc(handle, reinterpret_cast<double _Complex*>(ch_mat),
                                         reinterpret_cast<double _Complex*>(Stmp),
                                         eigen, reinterpret_cast<double _Complex*>(lowf.wfc_k[ik].c), 0, &elpa_derror);
        ModuleBase::timer::tick("Diago_LCAO_Matrix","elpa_solve");

        // the eigenvalues.
        BlasConnector::copy(GlobalV::NBANDS, eigen, inc, ekb, inc);
        delete[] eigen;

        lowf.wfc_2d_to_grid(this->out_wfc_lcao, lowf.wfc_k[ik].c, lowf.wfc_k_grid[ik], ik);
        
    } // GenELPA method
	else if(GlobalV::KS_SOLVER=="scalapack_gvx")
	{
		diag_scalapack_gvx.pzhegvx_diag(pv->desc, pv->ncol, pv->nrow, ch_mat, cs_mat, ekb, lowf.wfc_k[ik]);		// Peize Lin add 2021.11.02

        lowf.wfc_2d_to_grid(this->out_wfc_lcao, lowf.wfc_k[ik].c, lowf.wfc_k_grid[ik], ik);

    }

#endif
	return;
}



#ifdef __MPI
#include "../src_pw/occupy.h"
void Pdiag_Double::gath_eig_complex(MPI_Comm comm,int n,std::complex<double> **cc,std::complex<double> *Z, const int &ik)
{
    ModuleBase::TITLE("Pdiag_Double","gath_eig_complex");
    time_t time_start = time(NULL);
    //GlobalV::ofs_running << " Start gath_eig_complex Time : " << ctime(&time_start);
    const Parallel_Orbitals* pv = this->ParaV;

    int i,j,k;
    int nprocs,myid;
    MPI_Status status;
    MPI_Comm_size(comm,&nprocs);
    MPI_Comm_rank(comm,&myid);

    std::complex<double> **ctot;

	// mohan add 2010-07-03
	// the occupied bands are useless
	// for calculating charge density.
	if(GlobalV::DRANK> pv->lastband_in_proc)
	{
		delete[] Z;
	}

	// first we need to collect all
	// the occupied bands.
	// GlobalV::NBANDS * GlobalV::NLOCAL	
	if(GlobalV::DRANK==0)
	{
		ctot = new std::complex<double>*[GlobalV::NBANDS];
    	for (int i=0; i<GlobalV::NBANDS; i++)
    	{
        	ctot[i] = new std::complex<double>[GlobalV::NLOCAL];
        	ModuleBase::GlobalFunc::ZEROS(ctot[i], GlobalV::NLOCAL);
    	}
    	ModuleBase::Memory::record("Pdiag_Double","ctot",GlobalV::NBANDS*GlobalV::NLOCAL,"cdouble");
	}

	k=0;
    if (myid==0)
    {
        // mohan add nbnd0 2010-07-02
        int nbnd0 = -1;
        if (GlobalV::NBANDS < pv->loc_sizes[0])
        {
			// means all bands in this processor
			// is needed ( is occupied)
            nbnd0 = GlobalV::NBANDS;
        }
        else
        {
			// means this processor only save
			// part of GlobalV::NBANDS.
            nbnd0 = pv->loc_sizes[0];
        }
        if(pv->testpb)GlobalV::ofs_running << " nbnd in processor 0 is " << nbnd0 << std::endl;

        for (i=0; i<nbnd0; i++)
        {
            for (j=0; j<GlobalV::NLOCAL; j++)
            {
				// change the order in processor 0.
				// the contribution from processor 0.
                ctot[k][j]=Z[j*pv->loc_sizes[0]+i];
            }
            k++;
        }
		// Z is useless in processor 0 now.
		delete[] Z;
    }
    MPI_Barrier(comm);

	for (i=1; i<= pv->lastband_in_proc; i++)
    {
        // mohan fix bug 2010-07-02
        // rows indicates the data structure of Z.
        // mpi_times indicates the data distribution
        // time, each time send a band.
        int rows = pv->loc_sizes[i];
        int mpi_times;
        if (i==pv->lastband_in_proc)
        {
            mpi_times = pv->lastband_number;
        }
        else
        {
            mpi_times = pv->loc_sizes[i];
        }
        if(pv->testpb)GlobalV::ofs_running << " nbnd in processor " << i << " is " << mpi_times << std::endl;
        if (myid==i)
        {
            for (j=0; j<mpi_times; j++)
            {
                int tag = j;
                std::complex<double> *send = new std::complex<double>[n];
                int count = 0;

                for (int m=0; m<rows*n; m++)
                {
                    if (m%rows==j)
                    {
                        send[count] = Z[m];
                        ++count;
                    }
                }

				// send the data to processor 0.
                MPI_Send(send,n,mpicomplex,0,tag,comm);

                delete[] send;
            }
			// third part to delete Z;
			delete[] Z;
        }
        else if (myid==0)
        {
            int col=0;
            for (j=0; j<mpi_times; j++)
            {
                std::complex<double> *ctmp = new std::complex<double>[GlobalV::NLOCAL];
                ModuleBase::GlobalFunc::ZEROS(ctmp, GlobalV::NLOCAL);
                int tag = j;
                
				// Processor 0 receive the data from other processors.
				MPI_Recv(ctmp,n,mpicomplex,i,tag,comm,&status);

                for (int m=0; m<GlobalV::NLOCAL; m++)
                {
                    ctot[k][m]=ctmp[m];
//					GlobalV::ofs_running << " receive Z=" << ctmp[m] << std::endl;
                }
                k++;

                delete[] ctmp;
            }
        }
        //MPI_Barrier(comm);
    }
    if(pv->testpb)ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running,"Final k",k);

	// output the wave function if required.
	// this is a bad position to output wave functions.
	// but it works!
	std::stringstream ss;
	ss << GlobalV::global_out_dir << "LOWF_K_" << ik+1 << ".dat";
    if(this->out_wfc_lcao)
	{
//		std::cout << " write the wave functions" << std::endl;
		WF_Local::write_lowf_complex( ss.str(), ctot, ik );//mohan add 2010-09-09        
	}

	// mohan add 2010-09-10
	// distribution of local wave functions 
	// to each processor.
	WF_Local::distri_lowf_complex( ctot, cc);
	
	// clean staff.
	if(GlobalV::DRANK==0)
	{
    	for (int i=0; i<GlobalV::NBANDS; i++)
    	{
        	delete[] ctot[i];
    	}
    	delete[] ctot;
	}

    time_t time_end = time(NULL);
//    GlobalV::ofs_running << " End   gath_eig_complex Time : " << ctime(&time_end);

	ModuleBase::GlobalFunc::OUT_TIME("gather eigenvalues",time_start,time_end);
    return;
}
#endif

#ifdef __MPI
void Pdiag_Double::gath_full_eig(MPI_Comm comm,int n,double **c,double *Z)
{
    ModuleBase::TITLE("Pdiag_Double","gath_full_eig");

    time_t time_start = time(NULL);
    //GlobalV::ofs_running << " Start gath_full_eig Time : " << ctime(&time_start);

    int i,j,k,incx=1;
    int *loc_sizes,loc_size,nprocs,myid;
    MPI_Status status;
    MPI_Comm_size(comm,&nprocs);
    MPI_Comm_rank(comm,&myid);
    loc_sizes =(int*)malloc(sizeof(int)*nprocs);
    loc_size=n/nprocs;
    for (i=0; i<nprocs; i++)
    {
        if (i<n%nprocs)
        {
            loc_sizes[i]=loc_size+1;
        }
        else
        {
            loc_sizes[i]=loc_size;
        }
    }
    if (myid==0)
    {
        k=0;
        for (i=0; i<loc_sizes[0]; i++)
        {
            for (j=0; j<n; j++)
            {
                // i : column index;
                // j : row index.
//              c[k*n+j]=Z[i*n+j];
                c[k][j]=Z[j*loc_sizes[0]+i];
                //			GlobalV::ofs_running << " Z=" << Z[i*n+j] << std::endl;
            }
            k++;
        }
    }
    MPI_Barrier(comm);

    for (i=1; i<nprocs; i++)
    {
        const int rows = loc_sizes[i];
        const int mpi_times = rows;
        if (myid==i)
        {
            for (j=0; j<mpi_times; j++)
            {
                int tag = j;

                double *send = new double[n];

                int count = 0;
                for (int m=0; m<rows*n; m++)
                {
                    if (m%rows==j)
                    {
                        send[count] = Z[m];
                        ++count;
                    }
                }

                MPI_Send(send,n,MPI_DOUBLE,0,tag,comm);

                delete[] send;
            }
        }
        else if (myid==0)
        {
            int col=0;
            for (j=0; j<mpi_times; j++)
            {
                double *ctmp = new double[n];
                ModuleBase::GlobalFunc::ZEROS(ctmp, n);
                int tag = j;
                MPI_Recv(ctmp,n,MPI_DOUBLE,i,tag,comm,&status);

                for (int m=0; m<n; m++)
                {
                    c[k][m]=ctmp[m];
//					GlobalV::ofs_running << " receive Z=" << ctmp[m] << std::endl;
                }
                k++;

                delete[] ctmp;
            }
        }
        MPI_Barrier(comm);
    }

    for (int i=0; i<GlobalV::NLOCAL; i++)
    {
        Parallel_Common::bcast_double(c[i],GlobalV::NLOCAL);
    }

    time_t time_end = time(NULL);
//    GlobalV::ofs_running << " End   gath_full_eig Time : " << ctime(&time_end);

    return;
}

void Pdiag_Double::gath_full_eig_complex(MPI_Comm comm,int n,std::complex<double> **c,std::complex<double> *Z)
{
    ModuleBase::TITLE("Pdiag_Double","gath_full_eig_complex");

    time_t time_start = time(NULL);
    //GlobalV::ofs_running << " Start gath_full_eig_complex Time : " << ctime(&time_start);

    int i,j,k,incx=1;
    int *loc_sizes,loc_size,nprocs,myid;
    MPI_Status status;
    MPI_Comm_size(comm,&nprocs);
    MPI_Comm_rank(comm,&myid);
    loc_sizes =(int*)malloc(sizeof(int)*nprocs);
    loc_size=n/nprocs;
    for (i=0; i<nprocs; i++)
    {
        if (i<n%nprocs)
        {
            loc_sizes[i]=loc_size+1;
        }
        else
        {
            loc_sizes[i]=loc_size;
        }
    }
    if (myid==0)
    {
        k=0;
        for (i=0; i<loc_sizes[0]; i++)
        {
            for (j=0; j<n; j++)
            {
                // i : column index;
                // j : row index.
//              c[k*n+j]=Z[i*n+j];
                c[k][j]=Z[j*loc_sizes[0]+i];
                //			GlobalV::ofs_running << " Z=" << Z[i*n+j] << std::endl;
            }
            k++;
        }
    }
    MPI_Barrier(comm);

	for (i=1; i<nprocs; i++)
	{
		const int rows = loc_sizes[i];
		const int mpi_times = rows;
		if (myid==i)
		{
			for (j=0; j<mpi_times; j++)
			{
				int tag = j;

				std::complex<double> *send = new std::complex<double>[n];

				int count = 0;
				for (int m=0; m<rows*n; m++)
				{
					if (m%rows==j)
					{
						send[count] = Z[m];
						++count;
					}
				}

				MPI_Send(send,n,mpicomplex,0,tag,comm);

				delete[] send;
			}
		}
		else if (myid==0)
		{
			int col=0;
			for (j=0; j<mpi_times; j++)
			{
				std::complex<double> *ctmp = new std::complex<double>[n];
				ModuleBase::GlobalFunc::ZEROS(ctmp, n);
				int tag = j;
				MPI_Recv(ctmp,n,mpicomplex,i,tag,comm,&status);

				for (int m=0; m<n; m++)
				{
					c[k][m]=ctmp[m];
					//					GlobalV::ofs_running << " receive Z=" << ctmp[m] << std::endl;
				}
				k++;

				delete[] ctmp;
			}
		}
		MPI_Barrier(comm);
	}

	for (int i=0; i<GlobalV::NLOCAL; i++)
	{
		Parallel_Common::bcast_complex_double(c[i],GlobalV::NLOCAL);
	}

    time_t time_end = time(NULL);
    //GlobalV::ofs_running << " End   gath_full_eig Time : " << ctime(&time_end);
	ModuleBase::GlobalFunc::OUT_TIME("gather full eigenvalues",time_start,time_end);

    return;
}
#endif

//LiuXh add 2021-09-06, clear memory, totwfc and WFC_GAMMA_aug not used now
#ifdef __MPI
#include "../src_pw/occupy.h"
void Pdiag_Double::gath_eig(MPI_Comm comm,int n,double *Z)
{
    ModuleBase::TITLE("Pdiag_Double","gath_eig");
    time_t time_start = time(NULL);
//  GlobalV::ofs_running << " Start gath_eig Time : " << ctime(&time_start);

    const Parallel_Orbitals* pv = this->ParaV;
    
    int i, j, k;
    int nprocs,myid;
    MPI_Status status;
    MPI_Comm_size(comm,&nprocs);
    MPI_Comm_rank(comm,&myid);

    double **ctot;

	// mohan add 2010-07-03
	// the occupied bands are useless
	// for calculating charge density.
	if(GlobalV::DRANK > pv->lastband_in_proc)
	{
		delete[] Z;
	}

	// first we need to collect all
	// the occupied bands.
	// GlobalV::NBANDS * GlobalV::NLOCAL	
	if(myid==0)
	{
		ctot = new double*[GlobalV::NBANDS];
    	for (int i=0; i<GlobalV::NBANDS; i++)
    	{
        	ctot[i] = new double[GlobalV::NLOCAL];
        	ModuleBase::GlobalFunc::ZEROS(ctot[i], GlobalV::NLOCAL);
    	}
    	ModuleBase::Memory::record("Pdiag_Double","ctot",GlobalV::NBANDS*GlobalV::NLOCAL,"double");
	}

    k=0;
    if (myid==0)
    {
        // mohan add nbnd0 2010-07-02
        int nbnd0 = -1;
        if (GlobalV::NBANDS < pv->loc_sizes[0])
        {
			// means all bands in this processor
			// is needed ( is occupied)
            nbnd0 = GlobalV::NBANDS;
        }
        else
        {
			// means this processor only save
			// part of GlobalV::NBANDS.
            nbnd0 = pv->loc_sizes[0];
        }
        if(pv->testpb) GlobalV::ofs_running << " nbnd in processor 0 is " << nbnd0 << std::endl;

//printf("from 0 to %d\n",nbnd0-1);
        for (i=0; i<nbnd0; i++)
        {
            for (j=0; j<GlobalV::NLOCAL; j++)
            {
				// change the order in processor 0.
				// the contribution from processor 0.
                ctot[k][j]=Z[j*pv->loc_sizes[0]+i];
            }
            k++;
        }
		// Z is useless in processor 0 now.
		delete[] Z;
    }
    MPI_Barrier(comm);


    for (i=1; i<= pv->lastband_in_proc; i++)
    {
        // mohan fix bug 2010-07-02
        // rows indicates the data structure of Z.
        // mpi_times indicates the data distribution
        // time, each time send a band.
        int rows = pv->loc_sizes[i];
        int mpi_times;
        if (i==pv->lastband_in_proc)
        {
            mpi_times = pv->lastband_number;
        }
        else
        {
            mpi_times = pv->loc_sizes[i];
        }
        if(pv->testpb)GlobalV::ofs_running << " nbnd in processor " << i << " is " << mpi_times << std::endl;
        if (myid==i)
        {
            for (j=0; j<mpi_times; j++)
            {
                int tag = j;
                double *send = new double[n];
                int count = 0;

                for (int m=0; m<rows*n; m++)
                {
                    if (m%rows==j)
                    {
                        send[count] = Z[m];
                        ++count;
                    }
                }

				// send the data to processor 0.
                MPI_Send(send,n,MPI_DOUBLE,0,tag,comm);

                delete[] send;
            }
			// third part to delete Z;
			delete[] Z;
        }
        else if (myid==0)
        {
            int col=0;
            for (j=0; j<mpi_times; j++)
            {
                double *ctmp = new double[GlobalV::NLOCAL];
                ModuleBase::GlobalFunc::ZEROS(ctmp, GlobalV::NLOCAL);
                int tag = j;
                
				// Processor 0 receive the data from other processors.
				MPI_Recv(ctmp,n,MPI_DOUBLE,i,tag,comm,&status);

                for (int m=0; m<GlobalV::NLOCAL; m++)
                {
                    ctot[k][m]=ctmp[m];
//					GlobalV::ofs_running << " receive Z=" << ctmp[m] << std::endl;
                }
                k++;

                delete[] ctmp;
            }
        }
        //MPI_Barrier(comm);
    }
    if(pv->testpb)ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running,"Final k",k);
/*
if(myid==0){
	double *vect=new double[GlobalV::NLOCAL*GlobalV::NBANDS];
	double *form=new double[GlobalV::NBANDS*GlobalV::NBANDS];
	int x,y;
	for(x=0;x<GlobalV::NBANDS;x++){
		for(y=0;y<GlobalV::NLOCAL;y++){
			vect[x*GlobalV::NLOCAL+y]=ctot[x][y];
		}
	}
		char chT='T';
		char chN='N';
		int	ne = GlobalV::NBANDS;
		int m1, n1, k1;
		double ONE=1.0,ZERO=0.0;
		m1 = ne;
		n1 = ne;
		k1 = GlobalV::NLOCAL;
		dgemm_(&chT, &chN, &m1, &n1, &k1, &ONE, vect, &k1, vect, &k1, &ZERO, form, &m1);
		double di=0.0,oth=0.0;
		for(x=0;x<ne;x++){
			for(y=0;y<ne;y++){
				if(x==y){
					di+=fabs(form[x*ne+y]);
				}else{
					oth+=fabs(form[x*ne+y]);
				}
			}
		}
		di-=ne;
		printf("\n\n\ndi=%.16lf\n\n\nother=%.16lf\n\n\n",di,oth);
		//assert(0>1);
}
*/
MPI_Barrier(comm);

	// mohan add 2010-09-10
	// output the wave function if required.
	// this is a bad position to output wave functions.
	// but it works!
    if(this->out_wfc_lcao)
	{
		// read is in ../src_algorithms/wf_local.cpp
		std::stringstream ss;
		ss << GlobalV::global_out_dir << "LOWF_GAMMA_S" << GlobalV::CURRENT_SPIN+1 << ".dat";
		// mohan add 2012-04-03, because we need the occupations for the
		// first iteration. 
		Occupy::calculate_weights();
		WF_Local::write_lowf( ss.str(), ctot );//mohan add 2010-09-09        
	}

	// mohan add 2010-09-10
	// distribution of local wave functions 
	// to each processor.
	// only used for GlobalV::GAMMA_ONLY_LOCAL
	//WF_Local::distri_lowf( ctot, wfc);


	// clean staff.
	if(myid==0)
	{
    	for (int i=0; i<GlobalV::NBANDS; i++)
    	{
        	delete[] ctot[i];
    	}
    	delete[] ctot;
	}

    time_t time_end = time(NULL);
    //GlobalV::ofs_running << " End   gath_eig Time : " << ctime(&time_end);
	ModuleBase::GlobalFunc::OUT_TIME("gather eigenvalues",time_start,time_end);
    return;
}
#endif
