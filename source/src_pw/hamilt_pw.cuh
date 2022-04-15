#ifndef HAMILT_PW_CUH
#define HAMILT_PW_CUH

#include "../module_base/global_function.h"
#include "../module_base/global_variable.h"
#include "../module_base/complexmatrix.h"

#ifdef __CUDA
#include "cublas_v2.h"
#include "cufft.h"
#include "cusolverDn.h"
#endif

#ifdef __ROCM
#include "hipblas.h"
#include "hipfft.h"
#include "rocsolver.h"
#endif

class Hamilt_PW
{

  public:
	Hamilt_PW();
	~Hamilt_PW();

	static int moved;

	void allocate(const int &npwx, // number of plane wave (max)
				  const int &npol, // polarization
				  const int &nkb, // number of non-local pseudopotential projectors
				  const int &nrxx); // number of grids on this processor

	void cal_err(const int &npw, ModuleBase::ComplexMatrix &psi, const int &nband, double *em, double *err);

	void init_k(const int ik);

	friend class Diago_David;
	// friend class Diago_CG;
	template <class T, class T2> friend class Diago_CG_CUDA;
	friend class Exx_Lip;
	friend class Hamilt;
	friend class Stochastic_Iter;

	// Use hpsi(cpu) in diagH_subspace.
	void diagH_subspace(const int ik,
						const int nstart,
						const int nbnd,
						const ModuleBase::ComplexMatrix &psi,
						ModuleBase::ComplexMatrix &evc,
						double *en);

	// Use hpsi_cuda(gpu) in diagH_subspace instead of hpsi(cpu).
	void diagH_subspace_cuda(const int ik,
							 const int nstart,
							 const int n_band,
							 const double2 *psi_c, // matrix
							 double2 *evc, // matrix
							 double *en,
							 double2 *vkb_c);

	// Calculate hpsi and spsi (FP32)
	void h_1psi_cuda(const int npw, const float2 *psi1d, float2 *hpsi, float2 *spsi, float2 *vkb_c);

	// Calculate hpsi and spsi (FP64)
	void h_1psi_cuda(const int npw, const double2 *psi1d, double2 *hpsi, double2 *spsi, double2 *vkb_c);

	// hpsi (cpu)
	void h_1psi(const int npw,
				const std::complex<double> *psi1d,
				std::complex<double> *hpsi,
				std::complex<double> *spsi);

	// Hpsi operation (FP32)
	void h_psi_cuda(const float2 *psi, float2 *hpsi, float2 *vkb_c, const int m = 1);

	// Hpsi operation (FP64)
	void h_psi_cuda(const double2 *psi, double2 *hpsi, double2 *vkb_c, const int m = 1);

	// In this version, spsi operation just copy psi into spsi.
	void s_1psi_cuda(const int npw, const float2 *psi, float2 *spsi);

	// In this version, spsi operation just copy psi into spsi.
	void s_1psi_cuda(const int npw, const double2 *psi, double2 *spsi);

	// In this version, spsi operation just copy psi into spsi.
	void s_1psi(const int npw, const std::complex<double> *psi, std::complex<double> *spsi);

	void h_psi(const std::complex<double> *psi, std::complex<double> *hpsi, const int m = 1);

  private:
	int *GR_index;

#ifdef __CUDA
	// vectors on device
	int *GR_index_d;
	double2 *d_becp;
	double2 *d_ps;

	cublasHandle_t hpw_handle; // cublas handle
#endif

	// add contributions of h*psi from
	// non-local pseduopotentials
	// add vnl_pp GPU version (FP32)
	void add_nonlocal_pp_cuda(float2 *hpsi_in, const float2 *becp, const float2 *d_vkb_c, const int m);

	// add vnl_pp GPU version (FP64)
	void add_nonlocal_pp_cuda(double2 *hpsi_in, const double2 *becp, const double2 *d_vkb_c, const int m);

	// CPU version
	void add_nonlocal_pp(std::complex<double> *hpsi, const std::complex<double> *becp, const int m);

  private:
	// Ddot operations.
	double ddot_real(const int &npw, const std::complex<double> *psi_L, const std::complex<double> *psi_R) const;

	std::complex<double> ddot(const int &npw,
							  const std::complex<double> *psi_L,
							  const std::complex<double> *psi_R) const;

	std::complex<double> just_ddot(const int &npw,
								   const std::complex<double> *psi_L,
								   const std::complex<double> *psi_R) const;

	std::complex<double> ddot(const int &npw,
							  const ModuleBase::ComplexMatrix &psi,
							  const int &m,
							  const std::complex<double> *psik) const;

  private:
	// Use Lapack to solve eigenvalue problem.
	void diag_zheev(const int &npw, ModuleBase::ComplexMatrix &psi, const int &nband, double *em, double *err);
};

#endif
