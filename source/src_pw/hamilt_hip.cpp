#include "global.h"
#include "hamilt.h"
#include "hip/hip_runtime.h"
// 
#include "diago_cg_hip.h"
#include "diago_david.h"
#include "../module_base/timer.h"
#include "hipfft.h"
using namespace HipCheck;

Hamilt::Hamilt()
{
}
Hamilt::~Hamilt()
{
}


__global__ void hamilt_cast_d2f(float *dst, double *src, int size)
{
	int i = blockIdx.x * blockDim.x + threadIdx.x;
	if (i < size)
	{
		dst[i] = __double2float_rn(src[i]);
	}
}

__global__ void hamilt_cast_f2d(double *dst, float *src, int size)
{
	int i = blockIdx.x * blockDim.x + threadIdx.x;
	if (i < size)
	{
		dst[i] = (double)(src[i]);
	}
}

__global__ void hamilt_cast_d2f(float2 *dst, double2 *src, int size)
{
	int i = blockIdx.x * blockDim.x + threadIdx.x;
	if (i < size)
	{
		dst[i].x = __double2float_rn(src[i].x);
		dst[i].y = __double2float_rn(src[i].y);
	}
}

__global__ void hamilt_cast_f2d(double2 *dst, float2 *src, int size)
{
	int i = blockIdx.x * blockDim.x + threadIdx.x;
	if (i < size)
	{
		dst[i].x = (double)(src[i].x);
		dst[i].y = (double)(src[i].y);
	}
}

bool Hamilt::test_exit_cond(const int &ntry, const int &notconv)
{
	//================================================================
	// If this logical function is true, need to do diagH_subspace
	// and cg again.
	//================================================================

	bool scf = true;
	if (GlobalV::CALCULATION == "nscf")
		scf = false;

	// If ntry <=5, try to do it better, if ntry > 5, exit.
	const bool f1 = (ntry <= 5);

	// In non-self consistent calculation, do until totally converged.
	const bool f2 = ((!scf && (notconv > 0)));

	// if self consistent calculation, if not converged > 5,
	// using diagH_subspace and cg method again. ntry++
	const bool f3 = ((scf && (notconv > 5)));
	return (f1 && (f2 || f3));
}

void Hamilt::diagH_subspace(const int ik,
							const int nstart,
							const int n_band,
							const ModuleBase::ComplexMatrix &psi,
							ModuleBase::ComplexMatrix &evc,
							double *en)
{
	if (nstart < n_band)
	{
		ModuleBase::WARNING_QUIT("diagH_subspace", "nstart < n_band!");
	}

	if (GlobalV::BASIS_TYPE == "pw" || GlobalV::BASIS_TYPE == "lcao_in_pw")
	{
		this->hpw.diagH_subspace(ik, nstart, n_band, psi, evc, en);
	}
	else
	{
		ModuleBase::WARNING_QUIT("diagH_subspace", "Check parameters: GlobalV::BASIS_TYPE. ");
	}
	return;
}

void Hamilt::diagH_subspace_cuda(const int ik,
								 const int nstart,
								 const int n_band,
								 const hipblasDoubleComplex *psi,
								 hipblasDoubleComplex *evc,
								 double *en,
								 hipblasDoubleComplex *d_vkb_c)
{
	if (nstart < n_band)
	{
		ModuleBase::WARNING_QUIT("diagH_subspace_cuda", "nstart < n_band!");
	}

	if (GlobalV::BASIS_TYPE == "pw" || GlobalV::BASIS_TYPE == "lcao_in_pw")
	{
		this->hpw.diagH_subspace_cuda(ik, nstart, n_band, psi, evc, en, d_vkb_c);
	}
	else
	{
		ModuleBase::WARNING_QUIT("diagH_subspace_cuda", "Check parameters: GlobalV::BASIS_TYPE. ");
	}
	return;
}

//====================================================================
// calculates eigenvalues and eigenvectors of the generalized problem
// Hv=eSv, with H hermitean matrix, S overlap matrix .
// On output both matrix are unchanged
// LAPACK version - uses both ZHEGV and ZHEGVX
//=====================================================================
void Hamilt::diagH_LAPACK(const int nstart,
						  const int nbands,
						  const ModuleBase::ComplexMatrix &hc, // nstart * nstart
						  const ModuleBase::ComplexMatrix &sc, // nstart * nstart
						  const int ldh, // nstart
						  double *e,
						  ModuleBase::ComplexMatrix &hvec) // nstart * n_band
{
	ModuleBase::TITLE("Hamilt", "diagH_LAPACK");
	ModuleBase::timer::tick("Hamilt", "diagH_LAPACK");

	// Print info ...
	// cout<<"in diagH_lapack"<<endl;
	// cout<<nstart<<endl;
	// cout<<nbands<<endl;
	// cout<<"hc: "<<hc.nr<<" "<<hc.nc<<endl;
	// cout<<"sc: "<<sc.nr<<" "<<sc.nc<<endl;
	// cout<<ldh<<endl;
	// cout<<"hvec: "<<hvec.nr<<" "<<hvec.nc<<endl;

	int lwork = 0;

	ModuleBase::ComplexMatrix sdum(nstart, ldh);
	ModuleBase::ComplexMatrix hdum;

	sdum = sc;

	const bool all_eigenvalues = (nstart == nbands);

	int nb = LapackConnector::ilaenv(1, "ZHETRD", "U", nstart, -1, -1, -1);

	if (nb < 1)
	{
		nb = std::max(1, nstart);
	}

	if (nb == 1 || nb >= nstart)
	{
		lwork = 2 * nstart; // mohan modify 2009-08-02
	}
	else
	{
		lwork = (nb + 1) * nstart;
	}

	std::complex<double> *work = new std::complex<double>[lwork];
	ModuleBase::GlobalFunc::ZEROS(work, lwork);

	//=====================================================================
	// input s and (see below) h are copied so that they are not destroyed
	//=====================================================================

	int info = 0;
	int rwork_dim;
	if (all_eigenvalues)
	{
		rwork_dim = 3 * nstart - 2;
	}
	else
	{
		rwork_dim = 7 * nstart;
	}

	double *rwork = new double[rwork_dim];
	ModuleBase::GlobalFunc::ZEROS(rwork, rwork_dim);

	if (all_eigenvalues)
	{
		//===========================
		// calculate all eigenvalues
		//===========================
		hvec = hc;
		LapackConnector::zhegv(1, 'V', 'U', nstart, hvec, ldh, sdum, ldh, e, work, lwork, rwork, info);
	}
	else
	{
		//=====================================
		// calculate only m lowest eigenvalues
		//=====================================
		int *iwork = new int[5 * nstart];
		int *ifail = new int[nstart];

		ModuleBase::GlobalFunc::ZEROS(rwork, 7 * nstart);
		ModuleBase::GlobalFunc::ZEROS(iwork, 5 * nstart);
		ModuleBase::GlobalFunc::ZEROS(ifail, nstart);

		hdum.create(nstart, ldh);
		hdum = hc;

		//=============================
		// Number of calculated bands
		//=============================
		int mm = nbands;

		LapackConnector::zhegvx(1, // INTEGER
								'V', // CHARACTER*1
								'I', // CHARACTER*1
								'U', // CHARACTER*1
								nstart, // INTEGER
								hdum, // COMPLEX*16 array
								ldh, // INTEGER
								sdum, // COMPLEX*16 array
								ldh, // INTEGER
								0.0, // DOUBLE PRECISION
								0.0, // DOUBLE PRECISION
								1, // INTEGER
								nbands, // INTEGER
								0.0, // DOUBLE PRECISION
								mm, // INTEGER
								e, // DOUBLE PRECISION array
								hvec, // COMPLEX*16 array
								ldh, // INTEGER
								work, // DOUBLE array, dimension (MAX(1,LWORK))
								lwork, // INTEGER
								rwork, // DOUBLE PRECISION array, dimension (7*N)
								iwork, // INTEGER array, dimension (5*N)
								ifail, // INTEGER array, dimension (N)
								info // INTEGER
		);

		delete[] iwork;
		delete[] ifail;
	}
	delete[] rwork;
	delete[] work;

	ModuleBase::timer::tick("Hamilt", "diagH_LAPACK");
	return;
}
