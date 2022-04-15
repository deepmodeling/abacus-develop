#include "LCAO_evolve.h"
#include "../src_pw/global.h"
#include "../src_pdiag/pdiag_double.h"
#include "../src_io/write_HS.h"
#include"../input.h"
#include <complex>
#include "../module_base/scalapack_connector.h"
//fuxiang add 2016-10-28

Evolve_LCAO_Matrix::Evolve_LCAO_Matrix(LCAO_Matrix* lm) :
    LM(lm)
{}
Evolve_LCAO_Matrix::~Evolve_LCAO_Matrix() {}

void Evolve_LCAO_Matrix::evolve_complex_matrix(const int &ik, Local_Orbital_wfc &lowf)const
{
	ModuleBase::TITLE("Evolve_LCAO_Matrix","evolve_complex_matrix");
	time_t time_start = time(NULL);
	GlobalV::ofs_running << " Start Time : " << ctime(&time_start);

	if (INPUT.tddft==1)
	{
///*
#ifdef __MPI
		this->using_ScaLAPACK_complex(ik, lowf.wfc_k[ik]);
#else
		//this->using_LAPACK_complex(ik, lowf.wfc_k_grid, lowf.wfc_k[ik]);
#endif
//*/
		//this->using_ScaLAPACK_complex(ik, lowf.wfc_k_grid, lowf.wfc_k[ik]);
	}
	else
	{
		ModuleBase::WARNING_QUIT("Evolve_LCAO_Matrix::evolve_complex_matrix","only tddft==1 cando evolve");
	}

	time_t time_end = time(NULL);
	ModuleBase::GlobalFunc::OUT_TIME("evolve(std::complex)", time_start, time_end);
	
	return;
}

void Evolve_LCAO_Matrix::using_LAPACK_complex(const int& ik, std::complex<double>*** wfc_k_grid, std::complex<double>** c_init)const
{
    ModuleBase::TITLE("Evolve_LCAO_Matrix","using_LAPACK_complex");

//	Calculate the U operator

	ModuleBase::ComplexMatrix Htmp(GlobalV::NLOCAL,GlobalV::NLOCAL);
	ModuleBase::ComplexMatrix Stmp(GlobalV::NLOCAL,GlobalV::NLOCAL);
	for(int i=0; i<GlobalV::NLOCAL; i++)
        {
        	for(int j=0; j<GlobalV::NLOCAL; j++)
                {
                	Htmp(i,j) = this->LM->Hloc2[i*GlobalV::NLOCAL+j];
                        Stmp(i,j) = this->LM->Sloc2[i*GlobalV::NLOCAL+j];
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
        Numerator = Idmat - 0.5*INPUT.mdp.md_dt*41.34*Denominator;
        Denominator = Idmat + 0.5*INPUT.mdp.md_dt*41.34*Denominator;

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
			std::cout << wfc_k_grid[ik][i][j] << "\t";
		}
		std::cout <<std::endl;
	}
	std::cout << std::endl;
*/

/*	for(int i=0; i<GlobalV::NBANDS; i++)
	{
                for(int j=0; j<GlobalV::NLOCAL; j++)
		{
			std::cout << wfc_k_grid[ik][i][j] << "\t";
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
			wfc_k_grid[ik][i][j] = ccc[j];
		}	
	}

/*	for(int i=0; i<GlobalV::NBANDS; i++)
	{
                for(int j=0; j<GlobalV::NLOCAL; j++)
		{
			std::cout << wfc_k_grid[ik][i][j] << "\t";
		}
		std::cout <<std::endl;
	}
*/

/*      std::cout << " c: " <<std::endl;
        for(int i=0; i<GlobalV::NLOCAL; i++)
        {
                for(int j=0; j<GlobalV::NLOCAL; j++)
                {
                        std::cout << wfc_k_grid[ik][i][j] <<"\t";
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
                        std::cout << wfc_k_grid[ik][i][j] << "\t";
                }
                std::cout <<std::endl;
        }
        std::cout << std::endl;
*/

//	delete[] work;
//	delete[] ipiv;

	return;
}

#ifdef __MPI
int Evolve_LCAO_Matrix::using_ScaLAPACK_complex(const int &ik, ModuleBase::ComplexMatrix &wfc_2d)const
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

	//cout << "begin1: " <<endl;

	char uplo = 'U';
	const int inc = 1;

	const int  one_int = 1;
	
	//int nprocs, myid;
	//MPI_status status;
	//MPI_Comm_size(comm_2D, &nprocs);
	//MPI_Comm_rank(comm_2D, &myid);

	int loc_pos;
	//complex<double>* Stmp = this->LM->Sdiag2;
	//complex<double>* Htmp1 = this->LM->Hdiag2;
	//complex<double>* Htmp2 = this->LM->Sdiag2;
	//complex<double>* Htmp3 = this->LM->Sdiag2;

	complex<double>* Stmp = new complex<double> [this->LM->ParaV->nloc];
	complex<double>* Htmp1 = new complex<double> [this->LM->ParaV->nloc];
	complex<double>* Htmp2 = new complex<double> [this->LM->ParaV->nloc];
	complex<double>* Htmp3 = new complex<double> [this->LM->ParaV->nloc];
	ModuleBase::GlobalFunc::ZEROS(Stmp,this->LM->ParaV->nloc);
	ModuleBase::GlobalFunc::ZEROS(Htmp1,this->LM->ParaV->nloc);
	ModuleBase::GlobalFunc::ZEROS(Htmp2,this->LM->ParaV->nloc);
	ModuleBase::GlobalFunc::ZEROS(Htmp3,this->LM->ParaV->nloc);
	
	//cout << "this->LM->Hloc2" << *this->LM->Hloc2 << endl;
	//cout << "*Htmp2: " << *Htmp2 << endl;

        double *eigen = new double[GlobalV::NLOCAL];
        ModuleBase::GlobalFunc::ZEROS(eigen, GlobalV::NLOCAL);

	//cout << "GlobalV::NLOCAL : " << GlobalV::NLOCAL << endl;
	//cout << "nloc : " << this->LM->ParaV->nloc << endl;

	//cout << "begin02:" <<endl;

        //assert(loc_size > 0);
        //complex<double>* Z = new complex<double>[this->loc_size * NLOCAL];
        //ModuleBase::GlobalFunc::ZEROS(Z, this->loc_size * NLOCAL);

	zcopy_(&this->LM->ParaV->nloc, this->LM->Sloc2.data(), &inc, Stmp, &inc);
	zcopy_(&this->LM->ParaV->nloc, this->LM->Hloc2.data(), &inc, Htmp1, &inc);
	zcopy_(&this->LM->ParaV->nloc, this->LM->Hloc2.data(), &inc, Htmp2, &inc);
	zcopy_(&this->LM->ParaV->nloc, this->LM->Hloc2.data(), &inc, Htmp3, &inc);

	//cout << "*Htmp2: " << *Htmp2 << endl;

	complex<double> alpha = {1.0, 0.0};
	char transa = 'N';
	int desca = 0; 
	complex<double> beta = {0.0, -0.5*0.02*41.34};  // this need modify
	int descc = 0;

	//cout << "begin03:" << endl;

        pzgeadd_(
		&transa,
		&GlobalV::NLOCAL, &GlobalV::NLOCAL,
		&alpha,
		Stmp, &one_int, &one_int, this->LM->ParaV->desc,
                &beta,
		Htmp1, &one_int, &one_int, this->LM->ParaV->desc);

	//beta = (0.0, 0.5)*INPUT.md_dt;
	beta = {0.0, 0.5*0.02*41.34}; // this need modify

	//cout << "*Htmp1: " << *Htmp1 << endl;
	//cout << "Htmp2: " << *Htmp2 << endl;

	pzgeadd_(
		&transa,
		&GlobalV::NLOCAL, &GlobalV::NLOCAL,
		&alpha,
		Stmp, &one_int, &one_int, this->LM->ParaV->desc, 
		&beta,
		Htmp2, &one_int, &one_int, this->LM->ParaV->desc);

	//Next, invert the denominator
	int *ipiv = new int[ this->LM->ParaV->nloc ];
	int info;

	//cout << "*Htmp2: " << *Htmp2 << endl;
	//cout << "begin04:" << endl;

	pzgetrf_(
		&GlobalV::NLOCAL, &GlobalV::NLOCAL, 
		Htmp2, &one_int, &one_int, this->LM->ParaV->desc,
		ipiv,  &info);

	int LWORK=-1, liWORK=-1;
	std::vector<std::complex<double>> WORK(1,0);
	std::vector<int> iWORK(1,0);


	//cout << "begin05:" << endl;

	pzgetri_(
		&GlobalV::NLOCAL, 
		Htmp2, &one_int, &one_int, this->LM->ParaV->desc,
		ipiv,  WORK.data(),  &LWORK, iWORK.data(), &liWORK, &info);

	LWORK = WORK[0].real();
	WORK.resize(LWORK, 0);
	liWORK = iWORK[0];
	iWORK.resize(liWORK, 0);

	pzgetri_(
		&GlobalV::NLOCAL, 
		Htmp2, &one_int, &one_int, this->LM->ParaV->desc,
		ipiv,  WORK.data(),  &LWORK, iWORK.data(), &liWORK, &info);

	//alpha = (1.0, 0.0);
	//beta = (0.0, 0.0);
	char transb = 'T'; //This place requires subsequent testing of different transb.
	int descb = 0; 

	double alpha_1[2] = {1.0, 0.0};
	double beta_1[2] = {0.0, 0.0};

	//cout << "*Htmp2: " << *Htmp2 << endl;
	//cout << "begin06:" << endl;
	

	pzgemm_(
		&transa, &transb,
		&GlobalV::NLOCAL, &GlobalV::NLOCAL, &GlobalV::NLOCAL,
		&alpha_1[0],
		Htmp2, &one_int, &one_int, this->LM->ParaV->desc,
		Htmp1, &one_int, &one_int, this->LM->ParaV->desc, 
		&beta_1[0],
		Htmp3, &one_int, &one_int, this->LM->ParaV->desc);

	

	//cout << "U_operator Success!!!" <<endl;

	pzgemv_(
		&transa,
		&GlobalV::NLOCAL, &GlobalV::NLOCAL, 
		&alpha_1[0],
		Htmp3, &one_int, &one_int, this->LM->ParaV->desc,
		wfc_2d.c, &one_int, &one_int, this->LM->ParaV->desc, &one_int, 
		&beta_1[0],
		wfc_2d.c, &one_int, &one_int, this->LM->ParaV->desc, &one_int
        );



        // the eigenvalues.
        //dcopy_(&NBANDS, eigen, &inc, ekb, &inc);
        delete[] eigen;

        // Z is delete in gath_eig
        //ModuleBase::timer::tick("Evolve_LCAO_Matrix","gath_eig_complex",'G');
	return 0;

}
#endif
