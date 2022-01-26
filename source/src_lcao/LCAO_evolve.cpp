#include "LCAO_evolve.h"
#include "../src_pw/global.h"
#include "../src_pdiag/pdiag_double.h"
#include "../src_io/write_HS.h"
#include"../input.h"
#include <complex>
#include "../module_base/lapack_connector.h"
#include "../module_base/scalapack_connector.h"
#include "../module_base/blas_connector.h"
//fuxiang add 2016-10-28

Evolve_LCAO_Matrix::Evolve_LCAO_Matrix(){}
Evolve_LCAO_Matrix::~Evolve_LCAO_Matrix(){}

void Evolve_LCAO_Matrix::evolve_complex_matrix(const int &ik, std::complex<double>** WFC_K, ModuleBase::ComplexMatrix &wfc_2d)const
{
	ModuleBase::TITLE("Evolve_LCAO_Matrix","evolve_complex_matrix");
	time_t time_start = time(NULL);
	GlobalV::ofs_running << " Start Time : " << ctime(&time_start);

	if (INPUT.tddft==1)
	{
/*
#ifdef __MPI
		this->using_ScaLAPACK_complex(ik, cc, cc_init);
#else
		this->using_LAPACK_complex(ik, cc, cc_init);
#endif
*/
		//this->using_LAPACK_complex(ik, cc, cc_init);
		this->using_ScaLAPACK_complex(ik, WFC_K, wfc_2d);
	}
	else
	{
		ModuleBase::WARNING_QUIT("Evolve_LCAO_Matrix::evolve_complex_matrix","only tddft==1 cando evolve");
	}

	time_t time_end = time(NULL);
	ModuleBase::GlobalFunc::OUT_TIME("evolve(std::complex)", time_start, time_end);
	
	return;
}

inline int globalIndex(int localIndex, int nblk, int nprocs, int myproc)
{
    int iblock, gIndex;
    iblock=localIndex/nblk;
    gIndex=(iblock*nprocs+myproc)*nblk+localIndex%nblk;
    return gIndex;
}

inline int work2WFC(
	int naroc[2],
	int nb,
	int dim0,
	int dim1,
	int iprow,
	int ipcol,
	std::complex<double>* work,
	std::complex<double>** WFC,
	std::complex<double>** WFCAUG)
{
    for(int j=0; j<naroc[1]; ++j)
    {
        int igcol=globalIndex(j, nb, dim1, ipcol);
        if(igcol>=GlobalV::NBANDS) continue;
        for(int i=0; i<naroc[0]; ++i)
        {
            int igrow=globalIndex(i, nb, dim0, iprow);
	        int mu_local=GlobalC::GridT.trace_lo[igrow];
            if(mu_local>=0)
            {
                WFC[igcol][mu_local]=work[j*naroc[0]+i];
            }
	        int mu_aug=GlobalC::LOWF.trace_aug[igrow];
            if(mu_aug>=0)
            {
                WFCAUG[igcol][mu_aug]=work[j*naroc[0]+i];
            }
        }
    }
    return 0;
}

inline int WFC2work(
	int naroc[2],
	int nb,
	int dim0,
	int dim1,
	int iprow,
	int ipcol,
	std::complex<double>* work,
	std::complex<double>** WFC,
	std::complex<double>** WFCAUG)
{
    for(int j=0; j<naroc[1]; ++j)
    {
        int igcol=globalIndex(j, nb, dim1, ipcol);
        if(igcol>=GlobalV::NBANDS) continue;
        for(int i=0; i<naroc[0]; ++i)
        {
            int igrow=globalIndex(i, nb, dim0, iprow);
	        int mu_local=GlobalC::GridT.trace_lo[igrow];
            if(mu_local>=0)
            {
                work[j*naroc[0]+i]=WFC[igcol][mu_local];
            }
	        int mu_aug=GlobalC::LOWF.trace_aug[igrow];
            if(mu_aug>=0)
            {
                work[j*naroc[0]+i]=WFCAUG[igcol][mu_aug];
            }
        }
    }
    return 0;
}

void Evolve_LCAO_Matrix::using_LAPACK_complex(const int &ik, std::complex<double>** c, std::complex<double>** c_init)const
{
	                                                                                                                     ModuleBase::TITLE("Evolve_LCAO_Matrix","using_LAPACK_complex");

//	Calculate the U operator

	ModuleBase::ComplexMatrix Htmp(GlobalV::NLOCAL,GlobalV::NLOCAL);
	ModuleBase::ComplexMatrix Stmp(GlobalV::NLOCAL,GlobalV::NLOCAL);
	for(int i=0; i<GlobalV::NLOCAL; i++)
        {
        	for(int j=0; j<GlobalV::NLOCAL; j++)
                {
                	Htmp(i,j) = GlobalC::LM.Hloc2[i*GlobalV::NLOCAL+j];
                        Stmp(i,j) = GlobalC::LM.Sloc2[i*GlobalV::NLOCAL+j];
                }
        }

/*
	std::cout << " Htmp: " <<std::endl;
        for(int i=0; i<GlobalV::NLOCAL; i++)
        {
                for(int j=0; j<GlobalV::NLOCAL; j++)
                {
                        std::cout << Htmp(i,j) <<"\t";
                }
                std::cout <<std::endl;
        }
        std::cout <<std::endl;

	std::cout << " Stmp: " <<std::endl;
        for(int i=0; i<GlobalV::NLOCAL; i++)
        {
                for(int j=0; j<GlobalV::NLOCAL; j++)
                {
                        std::cout << Stmp(i,j) <<"\t";
                }
                std::cout <<std::endl;
        }
        std::cout <<std::endl;
*/
				
	int INFO;

        int LWORK=3*GlobalV::NLOCAL-1; //tmp
        std::complex<double> * WORK = new std::complex<double>[LWORK];
        ModuleBase::GlobalFunc::ZEROS(WORK, LWORK);
        int IPIV[GlobalV::NLOCAL];

        LapackConnector::zgetrf( GlobalV::NLOCAL, GlobalV::NLOCAL, Stmp, GlobalV::NLOCAL, IPIV, &INFO);
        LapackConnector::zgetri( GlobalV::NLOCAL, Stmp, GlobalV::NLOCAL, IPIV, WORK, LWORK, &INFO);

/*
        std::cout << " S^-1: " <<std::endl;
        for(int i=0; i<GlobalV::NLOCAL; i++)
        {
                for(int j=0; j<GlobalV::NLOCAL; j++)
                {
                        std::cout << Stmp(i,j) <<"\t";
                }
                std::cout <<std::endl;
        }
        std::cout <<std::endl;
*/

	ModuleBase::ComplexMatrix S_plus_H(GlobalV::NLOCAL,GlobalV::NLOCAL);
	S_plus_H = Stmp*Htmp;

/*
        std::cout << " S^-1*H: " <<std::endl;
        for(int i=0; i<GlobalV::NLOCAL; i++)
        {
                for(int j=0; j<GlobalV::NLOCAL; j++)
                {
                        std::cout << S_plus_H(i,j) <<"\t";
                }
                std::cout <<std::endl;
        }
        std::cout <<std::endl;
*/

	ModuleBase::ComplexMatrix Denominator(GlobalV::NLOCAL,GlobalV::NLOCAL);
	for (int i=0; i<GlobalV::NLOCAL; i++)
       	{
               	for (int j=0; j<GlobalV::NLOCAL; j++)
                {
                     /*   real(Denominator(i,j)) = -imag(S_plus_H(i,j));
                        imag(Denominator(i,j)) = real(S_plus_H(i,j));*/
                          Denominator(i,j) = std::complex<double>( -S_plus_H(i,j).imag(), S_plus_H(i,j).real() );
                }
        }
        
        ModuleBase::ComplexMatrix Idmat(GlobalV::NLOCAL,GlobalV::NLOCAL);
        for(int i=0; i<GlobalV::NLOCAL; i++)
        {
                for(int j=0; j<GlobalV::NLOCAL; j++)
                {
                        if(i==j) Idmat(i,j) = std::complex<double>(1.0, 0.0);
                       	else Idmat(i,j) = std::complex<double>(0.0, 0.0);
                }
        }
        //double delta_t;
//      delta_t = 0.2;	//identity: fs;
        ModuleBase::ComplexMatrix Numerator(GlobalV::NLOCAL,GlobalV::NLOCAL);
        Numerator = Idmat - 0.5*INPUT.mdp.dt*41.34*Denominator;
        Denominator = Idmat + 0.5*INPUT.mdp.dt*41.34*Denominator;

	int info;
        int lwork=3*GlobalV::NLOCAL-1; //tmp
        std::complex<double> * work = new std::complex<double>[lwork];
        ModuleBase::GlobalFunc::ZEROS(work, lwork);
        int ipiv[GlobalV::NLOCAL];
        
        LapackConnector::zgetrf( GlobalV::NLOCAL, GlobalV::NLOCAL, Denominator, GlobalV::NLOCAL, ipiv, &info);
        LapackConnector::zgetri( GlobalV::NLOCAL, Denominator, GlobalV::NLOCAL, ipiv, work, lwork, &info);

        ModuleBase::ComplexMatrix U_operator(GlobalV::NLOCAL,GlobalV::NLOCAL);
/*
        std::cout << " Numerator: " <<std::endl;
        for(int i=0; i<GlobalV::NLOCAL; i++)
        {
                for(int j=0; j<GlobalV::NLOCAL; j++)
                {
                        std::cout << Numerator(i,j) <<"\t";
                }
                std::cout <<std::endl;
        }
        std::cout <<std::endl;

        std::cout << " Denominator: " <<std::endl;
        for(int i=0; i<GlobalV::NLOCAL; i++)
        {
                for(int j=0; j<GlobalV::NLOCAL; j++)
                {
                        std::cout << Denominator(i,j) <<"\t";
                }
                std::cout <<std::endl;
        }
        std::cout <<std::endl;
*/

        U_operator = Numerator*Denominator;
/*
	std::cout << "U_operator Success!!!" <<std::endl;
        for(int i=0; i<GlobalV::NLOCAL; i++)
        {
                for(int j=0; j<GlobalV::NLOCAL; j++)
                {
                        std::cout << U_operator(i,j) <<"\t";
                }
                std::cout <<std::endl;
        }
	std::cout <<std::endl;
*/

// Calculate wave function at t+delta t
				
//	std::cout << "wave function coe at t+delta t !" << std::endl;

/*	for(int i=0; i<GlobalV::NBANDS; i++)
	{
		for(int j=0; j<GlobalV::NLOCAL; j++)
		{
			std::cout << c[i][j] << "\t";
		}
		std::cout <<std::endl;
	}
	std::cout << std::endl;
*/

/*	for(int i=0; i<GlobalV::NBANDS; i++)
	{
                for(int j=0; j<GlobalV::NLOCAL; j++)
		{
			std::cout << GlobalC::LOWF.WFC_K[ik][i][j] << "\t";
		}
		std::cout <<std::endl;
	}
	std::cout <<std::endl;
*/

	for(int i=0; i<GlobalV::NBANDS; i++)
	{
		std::complex<double> ccc[GlobalV::NLOCAL];
		for(int j=0; j<GlobalV::NLOCAL; j++)
		{	
			ccc[j] = std::complex<double>(0.0,0.0);
			for(int k=0; k<GlobalV::NLOCAL; k++)
			{
				 ccc[j] += U_operator(j,k)*c_init[i][k];
			}
		}
		for(int j=0; j<GlobalV::NLOCAL; j++)
		{
			c[i][j] = ccc[j];
			GlobalC::LOWF.WFC_K[ik][i][j] = ccc[j];
		}	
	}

/*	for(int i=0; i<GlobalV::NBANDS; i++)
	{
                for(int j=0; j<GlobalV::NLOCAL; j++)
		{
			std::cout << GlobalC::LOWF.WFC_K[ik][i][j] << "\t";
		}
		std::cout <<std::endl;
	}
*/

/*      std::cout << " c: " <<std::endl;
        for(int i=0; i<GlobalV::NLOCAL; i++)
        {
                for(int j=0; j<GlobalV::NLOCAL; j++)
                {
                        std::cout << c[i][j] <<"\t";
                }
                std::cout <<std::endl;
        }
        std::cout <<std::endl;
*/
/*
	for(int i=0; i<GlobalV::NBANDS; i++)
        {
                for(int j=0; j<GlobalV::NLOCAL; j++)
                {
                        std::cout << c[i][j] << "\t";
                }
                std::cout <<std::endl;
        }
        std::cout << std::endl;
*/

//	delete[] work;
//	delete[] ipiv;

	return;
}

int Evolve_LCAO_Matrix::using_ScaLAPACK_complex(const int &ik, complex<double>** WFC_K, ModuleBase::ComplexMatrix &wfc_2d)const
{
	ModuleBase::TITLE("Evolve_LCAO_Matrix","using_ScaLAPACK_complex");

	//inverse of matrix
	//pzgetrf (int *m, int *n, Complex16 *a, int ia, int ja, int *desca, int *ipiv, int info);
	//pzgetri (int *n, Complex16 *a, int *ia, int ja, int *desca, int *ipiv, Complex16 *Work, int *lwork, int *iwork, int *liwork, int *info);

	//product of vector and matrix
	//pzgemv(const char *trans, const int *m, const int *n, const Complex16 *alpha, const Complex16 *a, 
	//	const int *ia, const int *ja, const int *desca, const Complex16*x, const int *ix, const int *jx,
	//	const int *descx, const int *incx, const Complex16 *beta, Complex16 *y, const int *iy, 
	//	const int *jy, const int *descy, const int *incy);

	//matrix-matrix sum
	//pzgeadd (const char *trans, const int *m, const int *n, const Complex16 *alpha, const Complex16 *a, const int *ia,
	//	const int *ja, const int *desca, const Complex16 *beta, Complex16 *c, const int *ic, const int *jc, const int *descc);

	//matrix-matrix product
	//pzgemm

	char uplo = 'U';
	const int inc = 1;

	const int  one_int = 1;
	const int zero_int = 0;
	
	int ncol,nrow;
	ncol = GlobalC::ParaO.ncol;
	nrow = GlobalC::ParaO.nrow;

	int loc_pos;

	//complex<double>* Stmp = new complex<double> [GlobalC::ParaO.nloc];
	//complex<double>* Htmp1 = new complex<double> [GlobalC::ParaO.nloc];
	//complex<double>* Htmp2 = new complex<double> [GlobalC::ParaO.nloc];
	//complex<double>* Htmp3 = new complex<double> [GlobalC::ParaO.nloc];
	std::vector<std::complex<double>> Stmp;
	std::vector<std::complex<double>> Htmp1;
	std::vector<std::complex<double>> Htmp2;
	std::vector<std::complex<double>> Htmp3;
	ModuleBase::GlobalFunc::ZEROS(Stmp.data(),GlobalC::ParaO.nloc);
	ModuleBase::GlobalFunc::ZEROS(Htmp1.data(),GlobalC::ParaO.nloc);
	ModuleBase::GlobalFunc::ZEROS(Htmp2.data(),GlobalC::ParaO.nloc);
	ModuleBase::GlobalFunc::ZEROS(Htmp3.data(),GlobalC::ParaO.nloc);

	zcopy_(&GlobalC::ParaO.nloc, GlobalC::LM.Sloc2.data(), &inc, Stmp.data(), &inc);
	zcopy_(&GlobalC::ParaO.nloc, GlobalC::LM.Hloc2.data(), &inc, Htmp1.data(), &inc);
	zcopy_(&GlobalC::ParaO.nloc, GlobalC::LM.Hloc2.data(), &inc, Htmp2.data(), &inc);
	zcopy_(&GlobalC::ParaO.nloc, GlobalC::LM.Hloc2.data(), &inc, Htmp3.data(), &inc);
	
	double tmp_dt;
	//tmp_dt = INPUT.mdp.dt;
	tmp_dt = 0.01;
	complex<double> alpha = {1.0, 0.0};
	char transa = 'N';
	int desca = 0; 
	complex<double> beta = {0.0, -0.5*tmp_dt*41.34};  // this need modify
	int descc = 0;

        pzgeadd_(
		&transa,
		&GlobalV::NLOCAL, &GlobalV::NLOCAL,
		&alpha,
		Stmp.data(), &one_int, &one_int, GlobalC::ParaO.desc,
                &beta,
		Htmp1.data(), &one_int, &one_int, GlobalC::ParaO.desc);	

	//beta = (0.0, 0.5)*INPUT.md_dt;
	beta = {0.0, 0.5*tmp_dt*41.34}; // this need modify

	pzgeadd_(
		&transa,
		&GlobalV::NLOCAL, &GlobalV::NLOCAL,
		&alpha,
		Stmp.data(), &one_int, &one_int, GlobalC::ParaO.desc, 
		&beta,
		Htmp2.data(), &one_int, &one_int, GlobalC::ParaO.desc);

	//Next, invert the denominator
	int *ipiv = new int[ GlobalC::ParaO.nloc ];
	ModuleBase::GlobalFunc::ZEROS(ipiv,GlobalC::ParaO.nloc);
	int info;

	pzgetrf_(
		&GlobalV::NLOCAL, &GlobalV::NLOCAL, 
		Htmp2.data(), &one_int, &one_int, GlobalC::ParaO.desc,
		ipiv,  &info);

	int LWORK=-1, liWORK=-1;
	std::vector<std::complex<double>> WORK(1,0);
	std::vector<int> iWORK(1,0);

	pzgetri_(
		&GlobalV::NLOCAL, 
		Htmp2.data(), &one_int, &one_int, GlobalC::ParaO.desc,
		ipiv,  WORK.data(),  &LWORK, iWORK.data(), &liWORK, &info);

	LWORK = WORK[0].real();
	WORK.resize(LWORK, 0);
	liWORK = iWORK[0];
	iWORK.resize(liWORK, 0);

	pzgetri_(
		&GlobalV::NLOCAL, 
		Htmp2.data(), &one_int, &one_int, GlobalC::ParaO.desc,
		ipiv,  WORK.data(),  &LWORK, iWORK.data(), &liWORK, &info);
	
	char transb = 'T';
	int descb = 0; 

	double alpha_1 = 1.0;
	double beta_1 = 0.0;
	std::complex<double> im(0,1);

	transa='T';
	transb='T';
	pzgemm_(
		&transa, &transb,
		&GlobalV::NLOCAL, &GlobalV::NLOCAL, &GlobalV::NLOCAL,
		&alpha_1,
		Htmp2.data(), &one_int, &one_int, GlobalC::ParaO.desc,
		Htmp1.data(), &one_int, &one_int, GlobalC::ParaO.desc, 
		&beta_1,
		Htmp3.data(), &one_int, &one_int, GlobalC::ParaO.desc);
	
	int nprocs, myid;
	long maxnloc,nloc; // maximum number of elements in local matrix
	nloc=GlobalC::ParaO.nloc;
	MPI_Comm comm_2D_0;
	MPI_Comm comm_2D;
	comm_2D_0 = MPI_COMM_WORLD;
	int dim0,dim1;
	dim0 = (int)sqrt((double)GlobalV::DSIZE); 
	while (GlobalV::DSIZE%dim0!=0)
	{
		dim0 = dim0 - 1;
	}
	assert(dim0 > 0);
	dim1=GlobalV::DSIZE/dim0;
	int dim[2];
	dim[0]=dim0;
	dim[1]=dim1;
	int period[2]={1,1};
	int reorder=0;
	MPI_Cart_create(comm_2D_0,2,dim,period,reorder,&comm_2D);
    MPI_Status status;
    MPI_Comm_size(comm_2D, &nprocs);
    MPI_Comm_rank(comm_2D, &myid);
    MPI_Reduce(&nloc, &maxnloc, 1, MPI_LONG, MPI_MAX, 0, comm_2D);
    MPI_Bcast(&maxnloc, 1, MPI_LONG, 0, comm_2D);
	int naroc[2]; // maximum number of row or column
	int nb;
	nb = GlobalC::ParaO.nb;
	std::complex<double> *work1=new std::complex<double>[maxnloc];
	ModuleBase::GlobalFunc::ZEROS(work1,maxnloc);
	std::complex<double> *work2=new std::complex<double>[maxnloc];
	ModuleBase::GlobalFunc::ZEROS(work2,maxnloc);
    for(int iprow=0; iprow<dim0; ++iprow)
    {
        for(int ipcol=0; ipcol<dim1; ++ipcol)
        {
            const int coord[2]={iprow, ipcol};
            int src_rank;
            MPI_Cart_rank(comm_2D, coord, &src_rank);
            if(myid==src_rank)
        	{
                naroc[0]=nrow;
                naroc[1]=ncol;
				info=WFC2work(naroc, nb, dim0, dim1, iprow, ipcol, work1, GlobalC::LOWF.WFC_K_backup[ik], GlobalC::LOWF.WFC_K_aug_backup[ik]);
			}
		}
	}
	
	alpha_1=1.0;
	beta_1=0.0;
	transa = 'T';
	transb = 'T';
	pzgemm_(
		&transa, &transb,
		&GlobalV::NLOCAL, &GlobalV::NLOCAL, &GlobalV::NLOCAL,
		&alpha_1,
		work1, &one_int, &one_int, GlobalC::ParaO.desc,
		Htmp3.data(), &one_int, &one_int, GlobalC::ParaO.desc, 
		&beta_1,
		work2, &one_int, &one_int, GlobalC::ParaO.desc);
	
	alpha=1.0;
	beta=0.0;
	pzgeadd_(
		&transa,
		&GlobalV::NLOCAL, &GlobalV::NLOCAL,
		&alpha,
		work2, &one_int, &one_int, GlobalC::ParaO.desc,
        &beta,
		work1, &one_int, &one_int, GlobalC::ParaO.desc);
	
	for(int iprow=0; iprow<dim0; ++iprow)
    {
        for(int ipcol=0; ipcol<dim1; ++ipcol)
        {
            const int coord[2]={iprow, ipcol};
            int src_rank;
            MPI_Cart_rank(comm_2D, coord, &src_rank);
            if(myid==src_rank)
        	{
				zcopy_(&nloc, work1, &inc, work2, &inc);
                naroc[0]=nrow;
                naroc[1]=ncol;
				}
			info=MPI_Bcast(naroc, 2, MPI_INT, src_rank, comm_2D);
            info=MPI_Bcast(work2, maxnloc, MPI_DOUBLE_COMPLEX, src_rank, comm_2D);
			info=work2WFC(naroc, nb,dim0, dim1, iprow, ipcol, work2, GlobalC::LOWF.WFC_K[ik], GlobalC::LOWF.WFC_K_aug[ik]);
		}
	}
	delete work1;
	delete work2;
	delete ipiv;
	return 0;

}
