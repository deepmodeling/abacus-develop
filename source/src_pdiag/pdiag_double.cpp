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
	#include "../module_base/scalapack_connector.h"
}
#include "../module_base/lapack_connector.h"
#endif

#ifdef __CUSOLVER_LCAO
#include "diag_cusolver.cuh"
#endif

#ifdef __CUSOLVER_LCAO
template <typename T>
void cusolver_helper_gather(const T* mat_loc, T* mat_glb, const Parallel_Orbitals* pv){
	int rank = 0;
    MPI_Comm_rank(pv->comm_2D, &rank);

	int maxncol;
	MPI_Allreduce(&pv->ncol, &maxncol, 1, MPI_INT, MPI_MAX, pv->comm_2D);

    MPI_Datatype datatype = std::is_same<T, double>::value ? MPI_DOUBLE : MPI_DOUBLE_COMPLEX;
	for(int i = 0; i < maxncol - 1; i++)
    {
		MPI_Gather(mat_loc + i*pv->nrow, pv->nrow, datatype,
                   mat_glb + i * pv->dim1 * pv->nrow, pv->nrow,
                   datatype, 0, pv->comm_2D);
	}

	int displs[pv->dim1], rcounts[pv->dim1];
	for (int j = 0; j < pv->dim1; j++ )
	{
        if (GlobalV::NLOCAL % pv->dim1 && j >= GlobalV::NLOCAL % pv->dim1)
        {
            rcounts[j] = 0;
            displs[j] = (GlobalV::NLOCAL % pv->dim1 ) * pv->nrow;
        } else {
            rcounts[j] = pv->nrow;
            displs[j] = j * pv->nrow ;
        }
	}
	MPI_Gatherv(mat_loc + (maxncol - 1) * pv->nrow, rcounts[rank], datatype,
                mat_glb + (maxncol - 1) * pv->dim1 * pv->nrow, rcounts, displs,
                datatype, 0, pv->comm_2D);
}

template<typename T>
void cusolver_helper_scatter(const T* mat_glb, T* mat_loc, const Parallel_Orbitals* pv){
    int rank = 0;
    MPI_Comm_rank(pv->comm_2D, &rank);
    MPI_Status status;
    MPI_Datatype datatype = std::is_same<T, double>::value ? MPI_DOUBLE : MPI_DOUBLE_COMPLEX;

	if (rank == 0){
		for (int i =0; i < GlobalV::NLOCAL; i++){
			if ((i % pv->dim1) == 0) continue;
			MPI_Send(mat_glb + i*pv->nrow, pv->nrow, datatype, i%pv->dim1, i/pv->dim1, pv->comm_2D);
		}
		for (int i =0; i < GlobalV::NLOCAL; i+=pv->dim1)
			memcpy(mat_loc + i/pv->dim1*pv->nrow, mat_glb + i*pv->nrow, pv->nrow*sizeof(T));
	} else {
		for (int i = 0; i < pv->ncol; i++)
			MPI_Recv(mat_loc + i*pv->nrow, pv->nrow, datatype, 0, i, pv->comm_2D, &status);
	}
}
#endif


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
//wfc_gamma has been replaced by psi, this part needs rewriting
/*
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
*/
    //delete[] Stmp; //LiuXh 20171109
#ifdef __CUSOLVER_LCAO
	if(GlobalV::KS_SOLVER=="cusolver")
	{
		int rank;
		MPI_Comm_rank(pv->comm_2D, &rank);

        lowf.wfc_gamma[ik].create(pv->ncol, pv->nrow);			// Fortran order
		std::vector<double> ekb_tmp(GlobalV::NLOCAL,0);

		ModuleBase::timer::tick("Diago_LCAO_Matrix","gath_HS");
		double *htot, *stot, *vtot;
		vtot = new double[GlobalV::NLOCAL * GlobalV::NLOCAL];
		ModuleBase::Memory::record("Pdiag_Basic","vtot",GlobalV::NLOCAL*GlobalV::NLOCAL,"double");
		if(rank==0)    // htot  NLOCAL * NLOCAL
		{
			htot = new double[GlobalV::NLOCAL * GlobalV::NLOCAL];
			ModuleBase::Memory::record("Pdiag_Basic","htot",GlobalV::NLOCAL*GlobalV::NLOCAL,"double");

			stot = new double[GlobalV::NLOCAL * GlobalV::NLOCAL];
			ModuleBase::Memory::record("Pdiag_Basic","stot",GlobalV::NLOCAL*GlobalV::NLOCAL,"double");
		}
		cusolver_helper_gather<double>(h_mat, htot, pv);
		cusolver_helper_gather<double>(s_mat, stot, pv);
		ModuleBase::timer::tick("Diago_LCAO_Matrix","gath_HS");

		ModuleBase::timer::tick("Diago_LCAO_Matrix","cusolver_gvd_solve");
		if (rank == 0){
            static Diag_Cusolver_gvd* diag_cusolver_gvd = new Diag_Cusolver_gvd();
			diag_cusolver_gvd->Dngvd_double(GlobalV::NLOCAL, GlobalV::NLOCAL, htot, stot, ekb_tmp.data(), vtot);
		}
		ModuleBase::timer::tick("Diago_LCAO_Matrix","cusolver_gvd_solve");
		ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running,"K-S equation was solved by cusolver");

		if (rank == 0)	memcpy( ekb, ekb_tmp.data(), sizeof(double)*GlobalV::NBANDS );
		MPI_Bcast(ekb, GlobalV::NBANDS, MPI_DOUBLE, 0, pv->comm_2D);
		ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running,"eigenvalues were copied to ekb");

        ModuleBase::timer::tick("Diago_LCAO_Matrix","DIVIDE_EIG");
		cusolver_helper_scatter<double>(vtot, lowf.wfc_gamma[ik].c, this->ParaV);
        ModuleBase::timer::tick("Diago_LCAO_Matrix","DIVIDE_EIG");

		double** wfc_grid = nullptr;    //output but not do "2d-to-grid" conversion
        lowf.wfc_2d_to_grid(this->out_wfc_lcao, lowf.wfc_gamma[ik].c, wfc_grid);

		if (rank == 0){
			delete[] htot;
			delete[] stot;
		}
		delete[] vtot;
	}
#endif	//__CUSOLVER_LCAO

#endif

#ifdef TEST_DIAG
	{
		//wfc_gamma has been replaced by psi
		//this part nees rewriting
		/*
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
		*/
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
#ifdef __CUSOLVER_LCAO
	if(GlobalV::KS_SOLVER=="cusolver")
	{
		int rank;
		MPI_Comm_rank(pv->comm_2D, &rank);
        lowf.wfc_k[ik].create(pv->ncol, pv->nrow);            // Fortran order
		std::vector<double> ekb_tmp(GlobalV::NLOCAL,0);

		ModuleBase::timer::tick("Diago_LCAO_Matrix","gath_HS");
		std::complex<double> *htot, *stot, *vtot;
		vtot = new std::complex<double>[GlobalV::NLOCAL * GlobalV::NLOCAL];
		ModuleBase::Memory::record("Pdiag_Basic","vtot",GlobalV::NLOCAL*GlobalV::NLOCAL,"cdouble");
		if(rank==0)    // htot  NLOCAL * NLOCAL
		{
			htot = new std::complex<double>[GlobalV::NLOCAL * GlobalV::NLOCAL];
			ModuleBase::Memory::record("Pdiag_Basic","htot",GlobalV::NLOCAL*GlobalV::NLOCAL,"cdouble");

			stot = new std::complex<double>[GlobalV::NLOCAL * GlobalV::NLOCAL];
			ModuleBase::Memory::record("Pdiag_Basic","stot",GlobalV::NLOCAL*GlobalV::NLOCAL,"cdouble");
		}
		cusolver_helper_gather<std::complex<double> >(ch_mat, htot, pv);
		cusolver_helper_gather<std::complex<double> >(cs_mat, stot, pv);
		ModuleBase::timer::tick("Diago_LCAO_Matrix","gath_HS");

		ModuleBase::timer::tick("Diago_LCAO_Matrix","cusolver_gvd_solve");
		if (rank == 0){
            static Diag_Cusolver_gvd* diag_cusolver_gvd = new Diag_Cusolver_gvd();
			diag_cusolver_gvd->Dngvd_complex(GlobalV::NLOCAL, GlobalV::NLOCAL, htot, stot, ekb_tmp.data(), vtot);
		}
		ModuleBase::timer::tick("Diago_LCAO_Matrix","cusolver_gvd_solve");
		ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running,"K-S equation was solved by cusolver");

		if (rank == 0)	memcpy( ekb, ekb_tmp.data(), sizeof(double)*GlobalV::NBANDS );
		MPI_Bcast(ekb, GlobalV::NBANDS, MPI_DOUBLE, 0, pv->comm_2D);
		ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running,"eigenvalues were copied to ekb");

		ModuleBase::timer::tick("Diago_LCAO_Matrix","DIVIDE_EIG");
		//cusolver_helper_scatter<std::complex<double> >(vtot, lowf.wfc_k[ik].c, pv); //lowf.wfc_k has been replaced by psi
		ModuleBase::timer::tick("Diago_LCAO_Matrix","DIVIDE_EIG");

		lowf.wfc_2d_to_grid(this->out_wfc_lcao, lowf.wfc_k[ik].c, lowf.wfc_k_grid[ik], ik);

		if (rank == 0){
			delete[] htot;
			delete[] stot;
		}
		delete[] vtot;
	}
#endif	//__CUSOLVER_LCAO

#endif
	return;
}

