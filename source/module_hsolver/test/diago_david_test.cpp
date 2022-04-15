#include"module_hsolver/diagodavid.h"
#include"module_hsolver/iterdiagcon.h"
#include"module_hamilt/hamiltpw.h"
#include "src_pw/hamilt_pw.h"
#include"diago_mock.h"
#include "module_psi/psi.h"
#include"gtest/gtest.h"
#include"mpi.h"

#include "../../module_base/inverse_matrix.h"
#include "../../module_base/lapack_connector.h"

#define CONVTHRESHOLD 1e-3
#define DETAILINFO false


/************************************************
*  unit test of class Diago_David
***********************************************/

/**
 * Class Diago_David is used to solve the eigenvalues
 * This unittest test the function Diago_David::diag()
 * with different examples.
 * 	- the hamilt matrix (npw=100,500,1000) produced by random with sparsity of 90%
 *  - the hamilt matrix (npw=100,500,1000) produced by random with sparsity of 50%
 *  - the hamilt matrix (npw=100,500,1000) produced by random with sparsity of 0%
 *  - the hamilt matrix read from "data-H"
 * 
 * The test is passed when the eignvalues are closed to these calculated by LAPACK.
 *  
 */

//use lapack to calcualte eigenvalue of matrix hm
//NOTE: after finish this function, hm stores the eigen vectors.
void lapackEigen(int &npw, ModuleBase::ComplexMatrix &hm, ModuleBase::ComplexMatrix &ev,double * e, bool outtime=false)
{
	int lwork = 2 * npw;
	std::complex<double> *work2= new std::complex<double>[lwork];
	double* rwork = new double[3*npw-2];
	int info = 0;

	ModuleBase::ComplexMatrix tmp = hm;

	clock_t start,end;
	start = clock();
	LapackConnector::zheev('V', 'U', npw, tmp, npw, e, work2, lwork, rwork, &info);
	end = clock();
	if (outtime) std::cout<<"Lapack Run time: "<<(double)(end - start) / CLOCKS_PER_SEC<<" S"<<std::endl;

	ev = tmp;

	delete [] rwork;
	delete [] work2;
}

class DiagoDavPrepare 
{
public:
	DiagoDavPrepare(int nband, int npw, int sparsity, int order,double eps,int maxiter):
		nband(nband),npw(npw),sparsity(sparsity),order(order),eps(eps),maxiter(maxiter) {} 

	int nband, npw, sparsity, order, maxiter, notconv;
	double eps, avg_iter;

	void CompareEigen(ModulePsi::Psi<std::complex<double>> &phi, double *precondition)
	{
		//calculate eigenvalues by LAPACK;
		double* e_lapack = new double[npw];
		ModuleBase::ComplexMatrix ev;
		lapackEigen(npw, DIAGOTEST::hmatrix, ev, e_lapack,DETAILINFO);

		//do Diago_David::diag()
		double* en = new double[npw];
		
		Hamilt_PW hpw;
		ModuleHamilt::Hamilt *phm;
		phm = new ModuleHamilt::HamiltPW;
		ModuleHSolver::DiagoDavid dav(&hpw, precondition);
		ModuleHSolver::DiagoDavid::PW_DIAG_NDIM = order;
		ModuleHSolver::IterDiagControl::PW_DIAG_NMAX = maxiter;
		ModuleHSolver::IterDiagControl::PW_DIAG_THR = eps;
		phi.fix_k(0);

		clock_t start,end;
		start = clock();
		dav.diag(phm,phi,en);
		end = clock();

		if (DETAILINFO) std::cout<<"diag Run time: "<<(double)(end - start) / CLOCKS_PER_SEC<<" S, notconv=" 
						   << notconv << ", avg_iter=" << avg_iter << std::endl;

		for(int i=0;i<nband;i++)
		{
			EXPECT_NEAR(en[i],e_lapack[i],CONVTHRESHOLD);
		}

		delete [] en;
		delete [] e_lapack;	
		delete phm;	
	}
};

class DiagoDavTest : public ::testing::TestWithParam<DiagoDavPrepare> {};

TEST_P(DiagoDavTest,RandomHamilt)
{
	DiagoDavPrepare ddp = GetParam();
	if (DETAILINFO) std::cout << "npw=" << ddp.npw << ", nband=" << ddp.nband << ", sparsity=" 
			  << ddp.sparsity << ", eps=" << ddp.eps << std::endl;

	HPsi hpsi(ddp.nband,ddp.npw,ddp.sparsity);
	DIAGOTEST::hmatrix = hpsi.hamilt();
	DIAGOTEST::npw = ddp.npw;
	ModulePsi::Psi<std::complex<double>> psi = hpsi.psi();	

	ddp.CompareEigen(psi,hpsi.precond());
}


INSTANTIATE_TEST_SUITE_P(VerifyDiag,DiagoDavTest,::testing::Values(
		//DiagoDavPrepare(int nband, int npw, int sparsity, int order,double eps,int maxiter)
        DiagoDavPrepare(10,100,0,4,1e-5,500),
        DiagoDavPrepare(10,100,5,4,1e-5,500),
        DiagoDavPrepare(10,100,9,4,1e-5,500),
		DiagoDavPrepare(20,500,0,4,1e-5,500),
        DiagoDavPrepare(20,500,5,4,1e-5,500),
        DiagoDavPrepare(20,500,9,4,1e-5,500),
        DiagoDavPrepare(10,1000,8,4,1e-5,500)
		//DiagoDavPrepare(20,2000,8,4,1e-5,500)
));

TEST(DiagoDavRealSystemTest,dataH)
{
	ModuleBase::ComplexMatrix hmatrix;
	std::ifstream ifs("data-H");
	DIAGOTEST::readh(ifs,hmatrix);
	ifs.close();
	DIAGOTEST::hmatrix = hmatrix;
	DIAGOTEST::npw = hmatrix.nc;

	DiagoDavPrepare ddp(10,DIAGOTEST::npw,0,2,1e-5,500);
	HPsi hpsi(10,DIAGOTEST::npw);
	ModulePsi::Psi<std::complex<double>> psi = hpsi.psi();	
	
	ddp.CompareEigen(psi,hpsi.precond());
}

int main(int argc, char **argv)
{
	MPI_Init(&argc, &argv);
	testing::InitGoogleTest(&argc, argv);
    int result = RUN_ALL_TESTS();
	MPI_Finalize();
	return result;
}
