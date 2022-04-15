#ifndef HAMILT_HIP_H
#define HAMILT_HIP_H

#include "../module_base/global_function.h"
#include "../module_base/global_variable.h"
#include "../module_base/complexmatrix.h"

#ifdef __ROCM
#include "hamilt_pw_hip.h"
#else
#include "hamilt_pw.h"
#endif

class Hamilt
{
  public:
	Hamilt();
	~Hamilt();

	void init_before_ions(void);

	void init_before_electrons(void);

	void clear_after_ions(void);

	// diagH_subspace then full space in pw
	void diagH_pw(const int &istep, const int &iter, const int &ik, const double *precondition, double &avg_iter);

	// generate H and S then call diagH_subspace
	void diagH_subspace(const int ik,
						const int nstart,
						const int nbnd,
						const ModuleBase::ComplexMatrix &psi,
						ModuleBase::ComplexMatrix &evc,
						double *en);

	// be called by diagH_subspace
	void diagH_LAPACK(const int n,
					  const int m,
					  const ModuleBase::ComplexMatrix &hc, // Hamiltonian matrix
					  const ModuleBase::ComplexMatrix &sc, // overlap matrix
					  const int ldh,
					  double *e, // output: eigenvalues
					  ModuleBase::ComplexMatrix &hvec); // output: eigenvectors

#ifdef __ROCM

	// rocsolver_handle cusolver_handle;
	void diagH_subspace_cuda(const int ik,
							 const int nstart,
							 const int n_band,
							 const hipblasDoubleComplex *psi,
							 hipblasDoubleComplex *evc,
							 double *en,
							 hipblasDoubleComplex *d_ekb_c);
/*
	void diagH_CUSOLVER(
		const int nstart,
		const int nbands,
		hipblasDoubleComplex* hc,  // nstart * nstart
		hipblasDoubleComplex* sc,  // nstart * nstart
		const int ldh, // nstart
		double *e,
		hipblasDoubleComplex* hvec);
*/
#endif

	Hamilt_PW hpw;

  private:
	bool test_exit_cond(const int &ntry, const int &notconv);
};
#endif
