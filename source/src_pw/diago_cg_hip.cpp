#include "diago_cg_hip.h"

#include "global.h"
#include "hip/hip_runtime.h"
using namespace HipCheck;

template <class T, class T2, class T3> int Diago_CG_CUDA<T, T2, T3>::moved = 0;

template <class T, class T2, class T3> __global__ void kernel_normalization(T3 *data, int size, T norm)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;
	if (idx < size)
	{
		data[idx].x /= norm;
		data[idx].y /= norm;
	}
}

template <class T, class T2, class T3>
__global__ void kernel_precondition(T3 *res, const T3 *data, const int size, const T *P)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;
	if (idx < size)
	{
		res[idx].x = data[idx].x / P[idx];
		res[idx].y = data[idx].y / P[idx];
	}
}

template <class T, class T2, class T3>
__global__ void kernel_precondition_inverse(T3 *res, const T3 *data, const int size, const T *P)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;
	if (idx < size)
	{
		res[idx].x = data[idx].x * P[idx];
		res[idx].y = data[idx].y * P[idx];
	}
}

template <class T, class T2, class T3> __global__ void kernel_get_gredient(T3 *g, T3 *ppsi, int size, T lambda)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;
	if (idx < size)
	{
		g[idx].x -= lambda * ppsi[idx].x;
		g[idx].y -= lambda * ppsi[idx].y;
	}
}

template <class T, class T2, class T3> __global__ void kernel_get_gammacg(int size, T3 *dst, const T3 *src, T gamma)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;
	if (idx < size)
	{
		dst[idx].x = dst[idx].x * gamma + src[idx].x;
		dst[idx].y = dst[idx].y * gamma + src[idx].y;
	}
}

template <class T, class T2, class T3> __global__ void kernel_get_normacg(int size, T3 *dst, const T3 *src, T norma)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;
	if (idx < size)
	{
		dst[idx].x = dst[idx].x - norma * src[idx].x;
		dst[idx].y = dst[idx].y - norma * src[idx].y;
	}
}

template <class T, class T2, class T3>
__global__ void kernel_multi_add(T3 *dst, T3 *src1, T a1, const T3 *src2, T a2, int size)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;
	if (idx < size)
	{
		dst[idx].x = src1[idx].x * a1 + src2[idx].x * a2;
		dst[idx].y = src1[idx].y * a1 + src2[idx].y * a2;
	}
}

template <class T, class T2, class T3> Diago_CG_CUDA<T, T2, T3>::Diago_CG_CUDA()
{
	test_cg = 0;
	CHECK_CUBLAS(hipblasCreate(&diag_handle));
	// CHECK_CUBLAS(hipblasCreate(&ddot_handle));
}

template <class T, class T2, class T3> Diago_CG_CUDA<T, T2, T3>::~Diago_CG_CUDA()
{
	CHECK_CUBLAS(hipblasDestroy(diag_handle));
	// CHECK_CUBLAS(hipblasDestroy(ddot_handle));
}

template <class T> void test_print(T *data, int size)
{
	T *h_data = (T *)malloc(size * sizeof(T));
	CHECK_CUDA(hipMemcpy(h_data, data, size * sizeof(T), hipMemcpyDeviceToHost));
	cout << sizeof(h_data[0]) << endl;
	for (int i = 0; i < size; i++)
	{
		cout << h_data[i].x << " " << h_data[i].y << endl;
	}
	delete[] h_data;
}

template <class T, class T2, class T3>
void Diago_CG_CUDA<T, T2, T3>::diag(T2 *phi, // matrix nband*dim
									T *e,
									T2 *vkb_c,
									const int &dim,
									const int &dmx,
									const int &n_band,
									const T *precondition,
									const T &eps,
									const int &maxter,
									const bool &reorder,
									int &notconv,
									double &avg_iter)
{
	if (test_cg == 1)
		ModuleBase::TITLE("Diago_CG_HIP", "diag");
	ModuleBase::timer::tick("Diago_CG_HIP", "diag");

	avg_iter = 0.0;
	notconv = 0;

	//-------------------------------------------------------------------
	// "poor man" iterative diagonalization of a complex hermitian matrix
	// through preconditioned conjugate gradient algorithm
	// Band-by-band algorithm with minimal use of memory
	// Calls h_1phi and s_1phi to calculate H|phi> and S|phi>
	// Works for generalized eigenvalue problem (US pseudopotentials) as well
	//-------------------------------------------------------------------

	T2 *sphi;
	T2 *scg;
	T2 *hphi;
	T2 *g;
	T2 *cg;
	T2 *g0;
	T2 *pphi;
	T2 *lagrange;
	T2 *phi_m;

	CHECK_CUDA(hipMalloc((void **)&sphi, dim * sizeof(T2)));
	CHECK_CUDA(hipMalloc((void **)&scg, dim * sizeof(T2)));
	CHECK_CUDA(hipMalloc((void **)&hphi, dim * sizeof(T2)));
	CHECK_CUDA(hipMalloc((void **)&g, dim * sizeof(T2)));
	CHECK_CUDA(hipMalloc((void **)&cg, dim * sizeof(T2)));
	CHECK_CUDA(hipMalloc((void **)&g0, dim * sizeof(T2)));
	CHECK_CUDA(hipMalloc((void **)&pphi, dim * sizeof(T2)));
	CHECK_CUDA(hipMalloc((void **)&lagrange, n_band * sizeof(T2)));
	CHECK_CUDA(hipMalloc((void **)&phi_m, dim * sizeof(T2)));

	for (int m = 0; m < n_band; m++)
	{
		if (test_cg > 2)
			GlobalV::ofs_running << "Diagonal Band : " << m << endl;
		hipDeviceSynchronize();
		CHECK_CUDA(hipMemcpy(phi_m, &phi[m * dmx], dim * sizeof(T2), hipMemcpyDeviceToDevice));
		GlobalC::hm.hpw.s_1psi_cuda(dim, phi_m, sphi);
		this->schmit_orth(dim, dmx, m, phi, sphi, phi_m);
		GlobalC::hm.hpw.h_1psi_cuda(dim, phi_m, hphi, sphi, vkb_c);
		T em_host = 0;
		em_host = ddot_real(dim, phi_m, hphi);

		CHECK_CUDA(hipMemcpy(&e[m], &em_host, sizeof(T), hipMemcpyHostToDevice));

		int iter = 0;
		T gg_last = 0.0;
		T cg_norm = 0.0;
		T theta = 0.0;
		bool converged = false;
		// cg iteration

		for (iter = 0; iter < maxter; iter++)
		{
			// cout<<"******iter:"<<iter<<"******"<<endl;
			this->calculate_gradient(precondition, dim, hphi, sphi, g, pphi);
			this->orthogonal_gradient(dim, dmx, g, scg, lagrange, phi, m);
			this->calculate_gamma_cg(iter,
									 dim,
									 precondition,
									 g,
									 scg,
									 g0,
									 cg,
									 gg_last,
									 cg_norm,
									 theta,
									 phi_m); // scg used as sg
			converged = this->update_psi(dim,
										 cg_norm,
										 theta,
										 pphi,
										 cg,
										 scg,
										 phi_m,
										 em_host,
										 eps,
										 hphi,
										 sphi,
										 vkb_c); // pphi is used as hcg
			CHECK_CUDA(hipMemcpy(&e[m], &em_host, sizeof(T), hipMemcpyHostToDevice));
			if (converged)
				break;
		} // end iter

		CHECK_CUDA(hipMemcpy(&phi[m * dmx], phi_m, dim * sizeof(T2), hipMemcpyDeviceToDevice));

		if (!converged)
		{
			++notconv;
		}

		avg_iter += static_cast<double>(iter) + 1.00;

		if (m > 0 && reorder)
		{
			ModuleBase::GlobalFunc::NOTE("reorder bands!");
			T *e_host;
			e_host = (T *)malloc(n_band * sizeof(T));
			ModuleBase::GlobalFunc::ZEROS(e_host, n_band);
			CHECK_CUDA(hipMemcpy(e_host, e, n_band * sizeof(T), hipMemcpyDeviceToHost));

			if (e_host[m] - e_host[m - 1] < -2.0 * eps)
			{
				// if the last calculated eigenvalue is not the largest...
				int i = 0;
				for (i = m - 2; i >= 0; i--)
				{
					if (e_host[m] - e_host[i] > 2.0 * eps)
						break;
				}
				i++;
				moved++;

				// last calculated eigenvalue should be in the i-th position: reorder
				T e0 = e_host[m];

				CHECK_CUDA(hipMemcpy(pphi, &phi[m * dmx], dim * sizeof(T2), hipMemcpyDeviceToDevice));

				for (int j = m; j >= i + 1; j--)
				{
					e_host[j] = e_host[j - 1];
					CHECK_CUDA(
						hipMemcpy(&phi[j * dmx], &phi[(j - 1) * dmx], dim * sizeof(T2), hipMemcpyDeviceToDevice));
				}

				e_host[i] = e0;

				CHECK_CUDA(hipMemcpy(&phi[i * dmx], pphi, dim * sizeof(T2), hipMemcpyDeviceToDevice));
				// this procedure should be good if only a few inversions occur,
				// extremely inefficient if eigenvectors are often in bad order
				// (but this should not happen)
			} // endif

			CHECK_CUDA(hipMemcpy(e, e_host, n_band * sizeof(T), hipMemcpyHostToDevice));
			delete[] e_host;
		} // end reorder

	} // end m

	avg_iter /= n_band;

	CHECK_CUDA(hipFree(lagrange));
	CHECK_CUDA(hipFree(pphi));
	CHECK_CUDA(hipFree(g0));
	CHECK_CUDA(hipFree(cg));
	CHECK_CUDA(hipFree(g));
	CHECK_CUDA(hipFree(hphi));
	CHECK_CUDA(hipFree(scg));
	CHECK_CUDA(hipFree(sphi));
	CHECK_CUDA(hipFree(phi_m));

	ModuleBase::timer::tick("Diago_CG_HIP", "diag");
	return;
} // end subroutine ccgdiagg

template <class T, class T2, class T3>
void Diago_CG_CUDA<T, T2, T3>::calculate_gradient(const T *precondition,
												  const int dim,
												  const T2 *hpsi,
												  const T2 *spsi,
												  T2 *g,
												  T2 *ppsi)
{
	if (test_cg == 1)
		ModuleBase::TITLE("Diago_CG_HIP", "calculate_gradient");
	ModuleBase::timer::tick("Diago_CG_HIP", "calculate_grad");

	int thread = 512;
	int block = (dim + thread - 1) / thread;

	// kernel_precondition(data, res, size, precondition)
	// (2) PH|psi> : g[i] = hpsi[i]/precondition[i]
	hipLaunchKernelGGL(HIP_KERNEL_NAME(kernel_precondition<T, T2, T3>),
					   dim3(block),
					   dim3(thread),
					   0,
					   0,
					   reinterpret_cast<T3 *>(g),
					   reinterpret_cast<const T3 *>(hpsi),
					   dim,
					   precondition);
	// (3) PS|psi> : ppsi[i] = spsi[i]/precondition[i]
	hipLaunchKernelGGL(HIP_KERNEL_NAME(kernel_precondition<T, T2, T3>),
					   dim3(block),
					   dim3(thread),
					   0,
					   0,
					   reinterpret_cast<T3 *>(ppsi),
					   reinterpret_cast<const T3 *>(spsi),
					   dim,
					   precondition);

	// Update lambda !
	// (4) <psi|SPH|psi >
	const T eh = this->ddot_real(dim, spsi, g);
	// (5) <psi|SPS|psi >
	const T es = this->ddot_real(dim, spsi, ppsi);
	const T lambda = eh / es;

	// Update g !
	hipLaunchKernelGGL(HIP_KERNEL_NAME(kernel_get_gredient<T, T2, T3>),
					   dim3(block),
					   dim3(thread),
					   0,
					   0,
					   reinterpret_cast<T3 *>(g),
					   reinterpret_cast<T3 *>(ppsi),
					   dim,
					   lambda);
	// hipLaunchKernelGGL(kernel_multi_add, dim3(block), dim3(thread), 0, 0, g, g, 1, ppsi, -lambda, dim);
	ModuleBase::timer::tick("Diago_CG_HIP", "calculate_grad");
	return;
}

template <class T, class T2, class T3>
void Diago_CG_CUDA<T, T2, T3>::orthogonal_gradient(const int &dim,
												   const int &dmx,
												   hipblasComplex *g,
												   hipblasComplex *sg,
												   hipblasComplex *lagrange,
												   const hipblasComplex *eigenfunction,
												   const int m)
{
	if (test_cg == 1)
		ModuleBase::TITLE("Diago_CG_HIP", "orthogonal_gradient");
	ModuleBase::timer::tick("Diago_CG_HIP", "orth_grad");

	GlobalC::hm.hpw.s_1psi_cuda(dim, g, sg);

	int inc = 1;

	hipblasOperation_t trans1 = HIPBLAS_OP_C;
	hipblasComplex ONE(1, 0);
	hipblasComplex ZERO(0, 0);
	hipblasComplex NEG_ONE(-1, 0);

	CHECK_CUBLAS(hipblasCgemv(diag_handle, trans1, dim, m, &ONE, eigenfunction, dmx, sg, inc, &ZERO, lagrange, inc));
	// Parallel_Reduce::reduce_complex_double_pool(lagrange, m); // todo
	// (3) orthogonal |g> and |Sg> to all states (0~m-1)
	hipblasOperation_t trans2 = HIPBLAS_OP_N;

	CHECK_CUBLAS(hipblasCgemv(diag_handle, trans2, dim, m, &NEG_ONE, eigenfunction, dmx, lagrange, inc, &ONE, g, inc));
	CHECK_CUBLAS(hipblasCgemv(diag_handle, trans2, dim, m, &NEG_ONE, eigenfunction, dmx, lagrange, inc, &ONE, sg, inc));

	ModuleBase::timer::tick("Diago_CG_HIP", "orth_grad");
	// CHECK_CUBLAS(hipblasDestroy(handle));
	return;
}

template <class T, class T2, class T3>
void Diago_CG_CUDA<T, T2, T3>::orthogonal_gradient(const int &dim,
												   const int &dmx,
												   hipblasDoubleComplex *g,
												   hipblasDoubleComplex *sg,
												   hipblasDoubleComplex *lagrange,
												   const hipblasDoubleComplex *eigenfunction,
												   const int m)
{
	if (test_cg == 1)
		ModuleBase::TITLE("Diago_CG_HIP", "orthogonal_gradient");
	ModuleBase::timer::tick("Diago_CG_HIP", "orth_grad");

	GlobalC::hm.hpw.s_1psi_cuda(dim, g, sg);

	int inc = 1;

	hipblasOperation_t trans1 = HIPBLAS_OP_C;

	hipblasDoubleComplex ONE(1, 0);
	hipblasDoubleComplex ZERO(0, 0);
	hipblasDoubleComplex NEG_ONE(-1, 0);

	CHECK_CUBLAS(hipblasZgemv(diag_handle, trans1, dim, m, &ONE, eigenfunction, dmx, sg, inc, &ZERO, lagrange, inc));

	// Parallel_Reduce::reduce_complex_double_pool(lagrange, m); // todo
	// (3) orthogonal |g> and |Sg> to all states (0~m-1)
	hipblasOperation_t trans2 = HIPBLAS_OP_N;

	CHECK_CUBLAS(hipblasZgemv(diag_handle, trans2, dim, m, &NEG_ONE, eigenfunction, dmx, lagrange, inc, &ONE, g, inc));
	CHECK_CUBLAS(hipblasZgemv(diag_handle, trans2, dim, m, &NEG_ONE, eigenfunction, dmx, lagrange, inc, &ONE, sg, inc));

	ModuleBase::timer::tick("Diago_CG_HIP", "orth_grad");
	// CHECK_CUBLAS(hipblasDestroy(handle));
	return;
}

template <class T, class T2, class T3>
void Diago_CG_CUDA<T, T2, T3>::calculate_gamma_cg(const int iter,
												  const int dim,
												  const T *precondition,
												  const T2 *g,
												  const T2 *sg,
												  T2 *psg,
												  T2 *cg,
												  T &gg_last,
												  const T &cg_norm,
												  const T &theta,
												  const T2 *psi_m)
{
	if (test_cg == 1)
		ModuleBase::TITLE("Diago_CG_HIP", "calculate_gamma_cg");
	ModuleBase::timer::tick("Diago_CG_HIP", "gamma_cg");
	T gg_inter;
	if (iter > 0)
	{
		// (1) Update gg_inter!
		// gg_inter = <g|psg>
		// Attention : the 'g' in psg is getted last time
		gg_inter = this->ddot_real(dim, g, psg); // b means before
	}

	int thread = 512;
	int block = (dim + thread - 1) / thread;
	hipLaunchKernelGGL(HIP_KERNEL_NAME(kernel_precondition_inverse<T, T2, T3>),
					   dim3(block),
					   dim3(thread),
					   0,
					   0,
					   reinterpret_cast<T3 *>(psg),
					   reinterpret_cast<const T3 *>(sg),
					   dim,
					   precondition);

	// (3) Update gg_now!
	// gg_now = < g|P|sg > = < g|psg >
	const T gg_now = this->ddot_real(dim, g, psg);

	if (iter == 0)
	{
		// (40) gg_last first value : equal gg_now
		gg_last = gg_now;
		// (50) cg direction first value : |g>
		// |cg> = |g>

		CHECK_CUDA(hipMemcpy(cg, g, dim * sizeof(T2), hipMemcpyDeviceToDevice));
	}
	else
	{
		// (4) Update gamma !
		assert(gg_last != 0.0);
		const T gamma = (gg_now - gg_inter) / gg_last;

		// (5) Update gg_last !
		gg_last = gg_now;

		// (6) Update cg direction !(need gamma and |go> ):
		hipLaunchKernelGGL(HIP_KERNEL_NAME(kernel_get_gammacg<T, T2, T3>),
						   dim3(block),
						   dim3(thread),
						   0,
						   0,
						   dim,
						   reinterpret_cast<T3 *>(cg),
						   reinterpret_cast<const T3 *>(g),
						   gamma);

		const T norma = gamma * cg_norm * sin(theta);
		hipLaunchKernelGGL(HIP_KERNEL_NAME(kernel_get_normacg<T, T2, T3>),
						   dim3(block),
						   dim3(thread),
						   0,
						   0,
						   dim,
						   reinterpret_cast<T3 *>(cg),
						   reinterpret_cast<const T3 *>(psi_m),
						   norma);
	}
	ModuleBase::timer::tick("Diago_CG_HIP", "gamma_cg");
	return;
}

template <class T, class T2, class T3>
bool Diago_CG_CUDA<T, T2, T3>::update_psi(const int dim,
										  T &cg_norm,
										  T &theta,
										  T2 *hcg,
										  const T2 *cg,
										  T2 *scg,
										  T2 *psi_m,
										  T &eigenvalue,
										  const T &threshold,
										  T2 *hpsi,
										  T2 *sphi,
										  T2 *vkb_c)
{
	if (test_cg == 1)
		ModuleBase::TITLE("Diago_CG_HIP", "update_psi");
	ModuleBase::timer::tick("Diago_CG_HIP", "update_psi");
	int thread = 512;
	int block = (dim + thread - 1) / thread;
	// pw.h_1psi(dim, cg, hcg, scg); // TODO
	// to cpu
	GlobalC::hm.hpw.h_1psi_cuda(dim, cg, hcg, scg, vkb_c);
	// hpsi end

	cg_norm = sqrt(this->ddot_real(dim, cg, scg));

	if (cg_norm < 1.0e-10)
		return 1;

	const T a0 = this->ddot_real(dim, psi_m, hcg) * 2.0 / cg_norm;
	const T b0 = this->ddot_real(dim, cg, hcg) / (cg_norm * cg_norm);

	const T e0 = eigenvalue;

	theta = atan(a0 / (e0 - b0)) / 2.0;

	const T new_e = (e0 - b0) * cos(2.0 * theta) + a0 * sin(2.0 * theta);

	const T e1 = (e0 + b0 + new_e) / 2.0;
	const T e2 = (e0 + b0 - new_e) / 2.0;
	if (e1 > e2)
	{
		theta += (T)ModuleBase::PI_HALF;
	}

	eigenvalue = min(e1, e2);

	const T cost = cos(theta);
	const T sint_norm = sin(theta) / cg_norm;
	hipLaunchKernelGGL(HIP_KERNEL_NAME(kernel_multi_add<T, T2, T3>),
					   dim3(block),
					   dim3(thread),
					   0,
					   0,
					   reinterpret_cast<T3 *>(psi_m),
					   reinterpret_cast<T3 *>(psi_m),
					   cost,
					   reinterpret_cast<const T3 *>(cg),
					   sint_norm,
					   dim);

	if (abs(eigenvalue - e0) < threshold)
	{
		ModuleBase::timer::tick("Diago_CG_HIP", "update_psi");
		return 1;
	}
	else
	{
		hipLaunchKernelGGL(HIP_KERNEL_NAME(kernel_multi_add<T, T2, T3>),
						   dim3(block),
						   dim3(thread),
						   0,
						   0,
						   reinterpret_cast<T3 *>(sphi),
						   reinterpret_cast<T3 *>(sphi),
						   cost,
						   reinterpret_cast<const T3 *>(scg),
						   sint_norm,
						   dim);
		hipLaunchKernelGGL(HIP_KERNEL_NAME(kernel_multi_add<T, T2, T3>),
						   dim3(block),
						   dim3(thread),
						   0,
						   0,
						   reinterpret_cast<T3 *>(hpsi),
						   reinterpret_cast<T3 *>(hpsi),
						   cost,
						   reinterpret_cast<const T3 *>(hcg),
						   sint_norm,
						   dim);
		ModuleBase::timer::tick("Diago_CG_HIP", "update_psi");
		return 0;
	}
}

template <class T, class T2, class T3>
void Diago_CG_CUDA<T, T2, T3>::schmit_orth(const int &dim,
										   const int &dmx,
										   const int &m, // end
										   const hipblasComplex *psi, // matrix
										   hipblasComplex *sphi,
										   hipblasComplex *psi_m)
{
	ModuleBase::timer::tick("Diago_CG_HIP", "schmit_orth");
	assert(m >= 0);
	// cout<<"orth, dim="<<dim<<endl;

	hipblasComplex *lagrange;
	CHECK_CUDA(hipMalloc((void **)&lagrange, (m + 1) * sizeof(hipblasComplex)));
	int inc = 1;
	int mp1 = m + 1;

	hipblasOperation_t trans1 = HIPBLAS_OP_C;

	hipblasComplex ONE(1, 0);
	hipblasComplex ZERO(0, 0);
	hipblasComplex NEG_ONE(-1, 0);
	CHECK_CUBLAS(hipblasCgemv(diag_handle, trans1, dim, mp1, &ONE, psi, dmx, sphi, inc, &ZERO, lagrange, inc));

	float psi_norm;
	CHECK_CUDA(hipMemcpy(&psi_norm, &lagrange[m], sizeof(float), hipMemcpyDeviceToHost));
	hipblasOperation_t trans2 = HIPBLAS_OP_N;
	CHECK_CUBLAS(hipblasCgemv(diag_handle, trans2, dim, m, &NEG_ONE, psi, dmx, lagrange, inc, &ONE, psi_m, inc));

	psi_norm -= ddot_real(m, lagrange, lagrange); // next
	psi_norm = sqrt(psi_norm);

	int thread = 512;
	int block = (dim + thread - 1) / thread;
	// float test_real = psi_m[0].real();
	hipLaunchKernelGGL(HIP_KERNEL_NAME(kernel_normalization<float, float2>),
					   dim3(block),
					   dim3(thread),
					   0,
					   0,
					   reinterpret_cast<float2 *>(psi_m),
					   dim,
					   psi_norm);

	GlobalC::hm.hpw.s_1psi_cuda(dim, psi_m, sphi);

	ModuleBase::timer::tick("Diago_CG_HIP", "schmit_orth");
	CHECK_CUDA(hipFree(lagrange));
	return;
}

template <class T, class T2, class T3>
void Diago_CG_CUDA<T, T2, T3>::schmit_orth(const int &dim,
										   const int &dmx,
										   const int &m, // end
										   const hipblasDoubleComplex *psi, // matrix
										   hipblasDoubleComplex *sphi,
										   hipblasDoubleComplex *psi_m)
{
	ModuleBase::timer::tick("Diago_CG_HIP", "schmit_orth");
	assert(m >= 0);

	hipblasDoubleComplex *lagrange;
	CHECK_CUDA(hipMalloc((void **)&lagrange, (m + 1) * sizeof(hipblasDoubleComplex)));
	int inc = 1;
	int mp1 = m + 1;

	hipblasOperation_t trans1 = HIPBLAS_OP_C;

	hipblasDoubleComplex ONE(1, 0);
	hipblasDoubleComplex ZERO(0, 0);
	hipblasDoubleComplex NEG_ONE(-1, 0);
	CHECK_CUBLAS(hipblasZgemv(diag_handle, trans1, dim, mp1, &ONE, psi, dmx, sphi, inc, &ZERO, lagrange, inc));

	double psi_norm;
	CHECK_CUDA(hipMemcpy(&psi_norm, &lagrange[m], sizeof(double), hipMemcpyDeviceToHost));
	hipblasOperation_t trans2 = HIPBLAS_OP_N;
	CHECK_CUBLAS(hipblasZgemv(diag_handle, trans2, dim, m, &NEG_ONE, psi, dmx, lagrange, inc, &ONE, psi_m, inc));

	psi_norm -= ddot_real(m, lagrange, lagrange); // next
	psi_norm = sqrt(psi_norm);

	int thread = 512;
	int block = (dim + thread - 1) / thread;
	hipLaunchKernelGGL(HIP_KERNEL_NAME(kernel_normalization<double, double2>),
					   dim3(block),
					   dim3(thread),
					   0,
					   0,
					   reinterpret_cast<double2 *>(psi_m),
					   dim,
					   psi_norm);

	GlobalC::hm.hpw.s_1psi_cuda(dim, psi_m, sphi);

	ModuleBase::timer::tick("Diago_CG_HIP", "schmit_orth");
	CHECK_CUDA(hipFree(lagrange));
	return;
}

template <class T, class T2, class T3>
float Diago_CG_CUDA<T, T2, T3>::ddot_real(const int &dim,
										  const hipblasComplex *psi_L,
										  const hipblasComplex *psi_R,
										  const bool reduce)
{
	int dim2 = 2 * dim;
	float result;
	CHECK_CUBLAS(hipblasSdot(diag_handle, dim2, (float *)psi_L, 1, (float *)psi_R, 1, &result));
	return result;
}

template <class T, class T2, class T3>
double Diago_CG_CUDA<T, T2, T3>::ddot_real(const int &dim,
										   const hipblasDoubleComplex *psi_L,
										   const hipblasDoubleComplex *psi_R,
										   const bool reduce)
{
	int dim2 = 2 * dim;
	double result;
	CHECK_CUBLAS(hipblasDdot(diag_handle, dim2, (double *)psi_L, 1, (double *)psi_R, 1, &result));
	return result;
}

template <class T, class T2, class T3>
hipblasComplex Diago_CG_CUDA<T, T2, T3>::ddot(const int &dim, const hipblasComplex *psi_L, const hipblasComplex *psi_R)
{
	hipblasComplex result;
	CHECK_CUBLAS(hipblasCdotc(diag_handle, dim, psi_L, 1, psi_R, 1, &result));
	return result;
} // end of ddot

template <class T, class T2, class T3>
hipblasDoubleComplex Diago_CG_CUDA<T, T2, T3>::ddot(const int &dim,
													const hipblasDoubleComplex *psi_L,
													const hipblasDoubleComplex *psi_R)
{
	hipblasDoubleComplex result;
	CHECK_CUBLAS(hipblasZdotc(diag_handle, dim, psi_L, 1, psi_R, 1, &result));
	return result;
} // end of ddot

template <class T, class T2, class T3>
hipblasComplex Diago_CG_CUDA<T, T2, T3>::ddot(const int &dim,
											  const hipblasComplex *psi, // complex
											  const int &m,
											  hipblasComplex *psik)
{
	hipblasComplex result;
	CHECK_CUBLAS(hipblasCdotc(diag_handle, dim, &psi[m * dim], 1, psik, 1, &result));
	return result;
} // end of ddot

template <class T, class T2, class T3>
hipblasDoubleComplex Diago_CG_CUDA<T, T2, T3>::ddot(const int &dim,
													const hipblasDoubleComplex *psi, // complex
													const int &m,
													hipblasDoubleComplex *psik)
{
	hipblasDoubleComplex result;
	CHECK_CUBLAS(hipblasZdotc(diag_handle, dim, &psi[m * dim], 1, psik, 1, &result));
	return result;
} // end of ddot

// this return <psi_L(m) | psi_R(n)>
template <class T, class T2, class T3>
hipblasComplex Diago_CG_CUDA<T, T2, T3>::ddot(const int &dim,
											  const hipblasComplex *psi_L,
											  const int &m,
											  const hipblasComplex *psi_R,
											  const int &n)
{
	hipblasComplex result;
	CHECK_CUBLAS(hipblasCdotc(diag_handle, dim, &psi_L[m * dim], 1, &psi_R[n * dim], 1, &result));
	return result;
} // end of ddot

template <class T, class T2, class T3>
hipblasDoubleComplex Diago_CG_CUDA<T, T2, T3>::ddot(const int &dim,
													const hipblasDoubleComplex *psi_L,
													const int &m,
													const hipblasDoubleComplex *psi_R,
													const int &n)
{
	hipblasDoubleComplex result;
	CHECK_CUBLAS(hipblasZdotc(diag_handle, dim, &psi_L[m * dim], 1, &psi_R[n * dim], 1, &result));
	return result;
} // end of ddot

template class Diago_CG_CUDA<double, hipblasDoubleComplex, double2>;
template class Diago_CG_CUDA<float, hipblasComplex, float2>;
