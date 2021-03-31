#include "LCAO_evolve.h"
#include "../src_pw/global.h"
#include "../src_pdiag/pdiag_double.h"
#include "../src_io/write_HS.h"
#include"../input.h"
#include <complex>
//fuxiang add 2016-10-28

Evolve_LCAO_Matrix::Evolve_LCAO_Matrix(){}
Evolve_LCAO_Matrix::~Evolve_LCAO_Matrix(){}

void Evolve_LCAO_Matrix::evolve_complex_matrix(const int &ik, complex<double>** cc, complex<double>** cc_init)const
{
	TITLE("Evolve_LCAO_Matrix","evolve_complex_matrix");
	time_t time_start = time(NULL);
	ofs_running << " Start Time : " << ctime(&time_start);

	if (tddft==1)
	{
/*
#ifdef __MPI
		this->using_ScaLAPACK_complex(ik, cc, cc_init);
#else
		this->using_LAPACK_complex(ik, cc, cc_init);
#endif
*/
		this->using_LAPACK_complex(ik, cc, cc_init);
	}
	else
	{
		WARNING_QUIT("Evolve_LCAO_Matrix::evolve_complex_matrix","only tddft==1 cando evolve");
	}

	time_t time_end = time(NULL);
	OUT_TIME("evolve(complex)", time_start, time_end);
	
	return;
}

void Evolve_LCAO_Matrix::using_LAPACK_complex(const int &ik, complex<double>** c, complex<double>** c_init)const
{
	                                                                                                                     TITLE("Evolve_LCAO_Matrix","using_LAPACK_complex");

//	Calculate the U operator

	ComplexMatrix Htmp(NLOCAL,NLOCAL);
	ComplexMatrix Stmp(NLOCAL,NLOCAL);
	for(int i=0; i<NLOCAL; i++)
        {
        	for(int j=0; j<NLOCAL; j++)
                {
                	Htmp(i,j) = LM.Hloc2[i*NLOCAL+j];
                        Stmp(i,j) = LM.Sloc2[i*NLOCAL+j];
                }
        }

/*
	cout << " Htmp: " <<endl;
        for(int i=0; i<NLOCAL; i++)
        {
                for(int j=0; j<NLOCAL; j++)
                {
                        cout << Htmp(i,j) <<"\t";
                }
                cout <<endl;
        }
        cout <<endl;

	cout << " Stmp: " <<endl;
        for(int i=0; i<NLOCAL; i++)
        {
                for(int j=0; j<NLOCAL; j++)
                {
                        cout << Stmp(i,j) <<"\t";
                }
                cout <<endl;
        }
        cout <<endl;
*/
				
	int INFO;
/*	LapackConnector::zpotrf('U',NLOCAL,Stmp,NLOCAL,&INFO);
	for (int i=0;i<NLOCAL; i++)
	{
		for (int j=0; j<NLOCAL; j++)
		{
			Stmp(i,j)=Stmp(j,i);
		}
	}
	LapackConnector::zpotri('U',NLOCAL,Stmp,NLOCAL,&INFO);
*/

        int LWORK=3*NLOCAL-1; //tmp
        complex<double> * WORK = new complex<double>[LWORK];
        ZEROS(WORK, LWORK);
        int IPIV[NLOCAL];

        LapackConnector::zgetrf( NLOCAL, NLOCAL, Stmp, NLOCAL, IPIV, &INFO);
        LapackConnector::zgetri( NLOCAL, Stmp, NLOCAL, IPIV, WORK, LWORK, &INFO);

/*
        cout << " S^-1: " <<endl;
        for(int i=0; i<NLOCAL; i++)
        {
                for(int j=0; j<NLOCAL; j++)
                {
                        cout << Stmp(i,j) <<"\t";
                }
                cout <<endl;
        }
        cout <<endl;
*/

	ComplexMatrix S_plus_H(NLOCAL,NLOCAL);
	S_plus_H = Stmp*Htmp;

/*
        cout << " S^-1*H: " <<endl;
        for(int i=0; i<NLOCAL; i++)
        {
                for(int j=0; j<NLOCAL; j++)
                {
                        cout << S_plus_H(i,j) <<"\t";
                }
                cout <<endl;
        }
        cout <<endl;
*/

	ComplexMatrix Denominator(NLOCAL,NLOCAL);
	for (int i=0; i<NLOCAL; i++)
       	{
               	for (int j=0; j<NLOCAL; j++)
                {
                     /*   real(Denominator(i,j)) = -imag(S_plus_H(i,j));
                        imag(Denominator(i,j)) = real(S_plus_H(i,j));*/
                          Denominator(i,j) = std::complex<double>( -S_plus_H(i,j).imag(), S_plus_H(i,j).real() );
                }
        }
        
        ComplexMatrix Idmat(NLOCAL,NLOCAL);
        for(int i=0; i<NLOCAL; i++)
        {
                for(int j=0; j<NLOCAL; j++)
                {
                        if(i==j) Idmat(i,j) = complex<double>(1.0, 0.0);
                       	else Idmat(i,j) = complex<double>(0.0, 0.0);
                }
        }
        double delta_t;
//      delta_t = 0.2;	//identity: fs;
        ComplexMatrix Numerator(NLOCAL,NLOCAL);
        Numerator = Idmat - 0.5*INPUT.md_dt*41.34*Denominator;
        Denominator = Idmat + 0.5*INPUT.md_dt*41.34*Denominator;

	int info;
        int lwork=3*NLOCAL-1; //tmp
        complex<double> * work = new complex<double>[lwork];
        ZEROS(work, lwork);
        int ipiv[NLOCAL];
        
        LapackConnector::zgetrf( NLOCAL, NLOCAL, Denominator, NLOCAL, ipiv, &info);
        LapackConnector::zgetri( NLOCAL, Denominator, NLOCAL, ipiv, work, lwork, &info);

        ComplexMatrix U_operator(NLOCAL,NLOCAL);
/*
        cout << " Numerator: " <<endl;
        for(int i=0; i<NLOCAL; i++)
        {
                for(int j=0; j<NLOCAL; j++)
                {
                        cout << Numerator(i,j) <<"\t";
                }
                cout <<endl;
        }
        cout <<endl;

        cout << " Denominator: " <<endl;
        for(int i=0; i<NLOCAL; i++)
        {
                for(int j=0; j<NLOCAL; j++)
                {
                        cout << Denominator(i,j) <<"\t";
                }
                cout <<endl;
        }
        cout <<endl;
*/

        U_operator = Numerator*Denominator;
/*
	cout << "U_operator Success!!!" <<endl;
        for(int i=0; i<NLOCAL; i++)
        {
                for(int j=0; j<NLOCAL; j++)
                {
                        cout << U_operator(i,j) <<"\t";
                }
                cout <<endl;
        }
	cout <<endl;
*/

// Calculate wave function at t+delta t
				
//	cout << "wave function coe at t+delta t !" << endl;

/*	for(int i=0; i<NBANDS; i++)
	{
		for(int j=0; j<NLOCAL; j++)
		{
			cout << c[i][j] << "\t";
		}
		cout <<endl;
	}
	cout << endl;
*/

/*	for(int i=0; i<NBANDS; i++)
	{
                for(int j=0; j<NLOCAL; j++)
		{
			cout << LOWF.WFC_K[ik][i][j] << "\t";
		}
		cout <<endl;
	}
	cout <<endl;
*/

	for(int i=0; i<NBANDS; i++)
	{
		complex<double> ccc[NLOCAL];
		for(int j=0; j<NLOCAL; j++)
		{	
			ccc[j] = (0.0,0.0);
			for(int k=0; k<NLOCAL; k++)
			{
				 ccc[j] += U_operator(j,k)*c_init[i][k];
			}
		}
		for(int j=0; j<NLOCAL; j++)
		{
			c[i][j] = ccc[j];
			LOWF.WFC_K[ik][i][j] = ccc[j];
		}	
	}

/*	for(int i=0; i<NBANDS; i++)
	{
                for(int j=0; j<NLOCAL; j++)
		{
			cout << LOWF.WFC_K[ik][i][j] << "\t";
		}
		cout <<endl;
	}
*/

/*      cout << " c: " <<endl;
        for(int i=0; i<NLOCAL; i++)
        {
                for(int j=0; j<NLOCAL; j++)
                {
                        cout << c[i][j] <<"\t";
                }
                cout <<endl;
        }
        cout <<endl;
*/
/*
	for(int i=0; i<NBANDS; i++)
        {
                for(int j=0; j<NLOCAL; j++)
                {
                        cout << c[i][j] << "\t";
                }
                cout <<endl;
        }
        cout << endl;
*/

//	delete[] work;
//	delete[] ipiv;

	return;
}

/*
#include <cmath>
#include <mpi.h>
extern "C"
{
    #include "/export/soft/intel2015/composer_xe_2015.1.133/mkl/include/mkl_pblas.h"
    #include "/export/soft/intel2015/composer_xe_2015.1.133/mkl/include/mkl_blacs.h"
    #include "/export/soft/intel2015/composer_xe_2015.1.133/mkl/include/mkl_scalapack.h"
    #include "/export/soft/intel2015/composer_xe_2015.1.133/mkl/include/mkl_blas.h"
}
#include "src_global/blas_connector.h"

int Evolve_LCAO_Matrix::using_ScaLAPACK_complex(const int &ik, complex<double>** c, complex<double>** c_init)const
{
	TITLE("Evolve_LCAO_Matrix", "using_scalapack_complex");
	ComplexMatrix Htmp(NLOCAL,NLOCAL);
	ComplexMatrix Stmp(NLOCAL,NLOCAL);
	for(int i=0; i<NLOCAL; i++)
        {
        		for(int j=0; j<NLOCAL; j++)
                	{
                		Htmp(i,j) = LM.Hloc2[i*NLOCAL+j];
                           	Stmp(i,j) = LM.Sloc2[i*NLOCAL+j];
                	}
        }

        for(int i=0; i<NLOCAL; i++)
        {
                for(int j=0; j<NLOCAL; j++)
                {
                        cout << Stmp(i,j) <<"\t";
                }
                cout <<endl;
        }
	cout <<endl;

	cout << "11111" << endl;

	int nFull, nblk;
	int myid, nprocs, myprow, nprows, mypcol, npcols;
	int my_blacs_ctxt;
	int mpierr, info;
	int narows, nacols;
	int MPIROOT=0;
	char BLACS_LAYOUT='R', BLACS_SCOPE='A';
	int ISRCPROC=0;
	int desc[9];
	
	nFull=NLOCAL;
	nblk=NLOCAL;	

	// set mpi enviroment
	int argc;
	char **argv;
//	MPI_Init(&argc, &argv);
//    	MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
//    	MPI_Comm_rank(MPI_COMM_WORLD, &myid);
	
	mpierr=MPI_Bcast(&nFull, 1, MPI_INT, MPIROOT, MPI_COMM_WORLD);
	mpierr=MPI_Bcast(&nblk, 1, MPI_INT, MPIROOT, MPI_COMM_WORLD);

	// set blacs parameters 
	for(npcols=int(sqrt(double(NPROC))); npcols>=2; --npcols)
    	{
        	if(NPROC%npcols==0) break;
    	}
    	nprows=NPROC/npcols;
	
	cout << "22222" <<endl;

//	Cblacs_get(MPI_COMM_WORLD, 0, &my_blacs_ctxt);
//    	Cblacs_gridinit(&my_blacs_ctxt, &BLACS_LAYOUT, nprows, npcols);
//	Cblacs_gridinfo(my_blacs_ctxt, &nprows, &npcols, &myprow, &mypcol);

//	blacs_get_(MPI_COMM_WORLD, 0, &my_blacs_ctxt);

//	blacs_gridinit_(&my_blacs_ctxt, &BLACS_LAYOUT, &nprows, &npcols);
//	blacs_gridinfo_(&my_blacs_ctxt, &nprows, &npcols, &myprow, &mypcol);


	narows=numroc_(&nFull, &nblk, &myprow, &ISRCPROC, &nprows);
    	nacols=numroc_(&nFull, &nblk, &mypcol, &ISRCPROC, &npcols);
    	descinit_(desc, &nFull, &nFull, &nblk, &nblk, &ISRCPROC, &ISRCPROC, &my_blacs_ctxt, &narows, &info);


	//init main matrix

//	int subMatrixSize=narows*nacols;

	int subMatrixSize=nFull*nFull;
	cout << "333333" << endl;

	// start
	int allinfo;
    	char transa, transb, side, uplo, diag;
    	int isrc=1, jsrc=1;
    	double alpha, beta;
    	int matrixCols=nacols;
    	int lda;
	int ipiv;
	int lwork, iwork, liwork;


	double *S, *work, *a, *b;
	S=new double[subMatrixSize];
	work=new double[subMatrixSize];

	MPI_Barrier(MPI_COMM_WORLD);

	cout << subMatrixSize << endl;

	for(int i=0; i<subMatrixSize; ++i) S[i]= real(LM.Sloc2[i]);
//	for(int i=0; i<subMatrixSize; ++i) cout << S[i] <<endl;

	for(int i=0; i<subMatrixSize; ++i) c[i]= real(LM.Sloc2[i]);

	int inc=1;
	DCOPY(&subMatrixSize, &c, &inc, S, &inc);


	cout << "44444" << endl;

//	pdgetrf(&nFull, &nFull, S, &isrc, &jsrc, desc, &ipiv, &info);

	cout << "55555"	<<endl;

	MPI_Allreduce(&info, &allinfo, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
	
	cout << "66666" <<endl;

	//calculate the inverse
//	pdgetri(&nFull, S, &isrc, &jsrc, desc, &ipiv, work, &lwork, &iwork, &liwork, &info);
	
	cout << "77777" << endl;

//    	pdgemm_(&S, &Htmp, &nFull, &nFull, &nFull,
//            	&alpha, work, &isrc, &jsrc, desc,
//                    	q,    &isrc, &jsrc, desc,
//            	&beta,  b,    &isrc, &jsrc, desc);

	

	MPI_Allreduce(&info, &allinfo, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

	cout << "88888"	<< endl;

	MPI_Barrier(MPI_COMM_WORLD);

	cout << "99999" << endl;

//	blacs_gridexit_(&my_blacs_ctxt);

//    	MPI_Finalize();

	cout <<"end!!!!!"<< endl;
    	return 0;
}

int globalIndex(int localIndex, int nblk, int nprocs, int myproc)
{
    int iblock, gIndex;
    iblock=localIndex/nblk;
    gIndex=(iblock*nprocs+myproc)*nblk+localIndex%nblk;
    return gIndex;
}


int localIndex(int globalIndex, int nblk, int nprocs, int& myproc)
{
    myproc=int((globalIndex%(nblk*nprocs))/nblk);
    return int(globalIndex/(nblk*nprocs))*nblk+globalIndex%nblk;
}

*/
