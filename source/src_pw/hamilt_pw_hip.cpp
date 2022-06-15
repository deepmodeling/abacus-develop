#include "global.h"
#include "hip/hip_runtime.h"
#include "../module_base/global_function.h"
#include "../module_base/global_variable.h"
// #include "hamilt_pw.cuh"
#include "../module_base/blas_connector.h"
#include "../src_io/optical.h" // only get judgement to calculate optical matrix or not.
#include "hamilt_pw_hip.h"
#include "myfunc.h"
#include "../module_base/timer.h"
using namespace HipCheck;

__global__ void cast_d2f(float *dst, double *src, int size)
{
	int i = blockIdx.x * blockDim.x + threadIdx.x;
	if (i < size)
	{
		dst[i] = __double2float_rn(src[i]);
	}
}

__global__ void cast_f2d(double *dst, float *src, int size)
{
	int i = blockIdx.x * blockDim.x + threadIdx.x;
	if (i < size)
	{
		dst[i] = (double)(src[i]);
	}
}

__global__ void cast_d2f(float2 *dst, double2 *src, int size)
{
	int i = blockIdx.x * blockDim.x + threadIdx.x;
	if (i < size)
	{
		dst[i].x = __double2float_rn(src[i].x);
		dst[i].y = __double2float_rn(src[i].y);
	}
}

__global__ void cast_f2d(double2 *dst, float2 *src, int size)
{
	int i = blockIdx.x * blockDim.x + threadIdx.x;
	if (i < size)
	{
		dst[i].x = (double)(src[i].x);
		dst[i].y = (double)(src[i].y);
	}
}

template <class T2> __global__ void kernel_copy(int size, T2 *dst, const T2 *src)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;
	if (idx < size)
	{
		dst[idx].x = src[idx].x;
		dst[idx].y = src[idx].y;
	}
}

template <class T, class T2> __global__ void kernel_get_tmhpsi(int size, T2 *dst, const T2 *src, T *g2kin)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;
	if (idx < size)
	{
		dst[idx].x = src[idx].x * g2kin[idx];
		dst[idx].y = src[idx].y * g2kin[idx];
	}
}

template <class T2> __global__ void kernel_add_tmhpsi(int size, T2 *dst, T2 *src, int *index)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;
	int p = index[idx];
	if (idx < size)
	{
		dst[idx].x += src[p].x;
		dst[idx].y += src[p].y;
	}
}

template <class T, class T2>
__global__ void kernel_addpp(T2 *ps, T *deeq, const T2 *becp, int nproj, int nprojx, int sum, int m, int nkb)
{
	int ip2 = blockDim.x * blockIdx.x + threadIdx.x;
	int ib = blockDim.y * blockIdx.y + threadIdx.y;
	if (ip2 < nproj && ib < m)
	{
		ps[(sum + ip2) * m + ib].x = ps[(sum + ip2) * m + ib].y = 0;

		for (int ip = 0; ip < nproj; ip++)
		{
			ps[(sum + ip2) * m + ib].x += deeq[ip * nprojx + ip2] * becp[ib * nkb + sum + ip].x;
			ps[(sum + ip2) * m + ib].y += deeq[ip * nprojx + ip2] * becp[ib * nkb + sum + ip].y;
		}
		// __syncthreads();
	}
}

template <class T> void print_test(T *data, int size)
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

int Hamilt_PW::moved = 0;

Hamilt_PW::Hamilt_PW()
{
	// hpsi = nullptr;
	// spsi = nullptr;
	// GR_index = nullptr;
#ifdef __ROCM
	hipMalloc((void **)&GR_index_d, sizeof(int));
	CHECK_CUBLAS(hipblasCreate(&hpw_handle));
#endif
}

Hamilt_PW::~Hamilt_PW()
{
	// delete[] hpsi;
	// delete[] spsi;
	delete[] GR_index;
#ifdef __ROCM
	CHECK_CUDA(hipFree(GR_index_d));
	if (hpw_handle)
		CHECK_CUBLAS(hipblasDestroy(hpw_handle));
		// CHECK_CUDA(hipFree(GR_index_d));
#endif
}

void Hamilt_PW::allocate(const int &npwx, const int &npol, const int &nkb, const int &nrxx)
{
	ModuleBase::TITLE("Hamilt_PW_HIP", "allocate");

	assert(npwx > 0);
	assert(npol > 0);
	assert(nkb >= 0);

	// delete[] hpsi;
	// delete[] spsi;
	// delete[] GR_index;

	delete[] GR_index;
	GR_index = new int[npwx];

#ifdef __ROCM
	if (GR_index_d)
	{
		CHECK_CUDA(hipFree(GR_index_d));
	}

	// CHECK_CUDA(hipMalloc((void**)&GR_index_d, npwx*sizeof(int)));
#endif
	// ZEROS(this->hpsi, npwx * npol);
	// ZEROS(this->spsi, npwx * npol);
	// ZEROS(this->GR_index, nrxx);

	return;
}

void Hamilt_PW::init_k(const int ik)
{
	ModuleBase::TITLE("Hamilt_PW_HIP", "init_k");
	// mohan add 2010-09-30
	// (1) Which spin to use.
	if (GlobalV::NSPIN == 2)
	{
		GlobalV::CURRENT_SPIN = GlobalC::kv.isk[ik];
	}

	// (3) Take the local potential.
	// cout<<"nrxx="<<GlobalC::wfcpw->nrxx<<endl;

	for (int ir = 0; ir < GlobalC::wfcpw->nrxx; ir++)
	{
		GlobalC::pot.vr_eff1[ir] = GlobalC::pot.vr_eff(GlobalV::CURRENT_SPIN, ir); // mohan add 2007-11-12
	}

	// (4) Calculate nonlocal pseudopotential vkb
	// if (GlobalC::ppcell.nkb > 0 && !LINEAR_SCALING) xiaohui modify 2013-09-02
	if (GlobalC::ppcell.nkb > 0
		&& (GlobalV::BASIS_TYPE == "pw" || GlobalV::BASIS_TYPE == "lcao_in_pw")) // xiaohui add 2013-09-02. Attention...
	{
		GlobalC::ppcell.getvnl(ik);
	}

	// (5) The number of wave functions.
	GlobalC::wf.npw = GlobalC::kv.ngk[ik];

	// (6) The index of plane waves.
	// int *GR_index_tmp = new int[GlobalC::wfcpw->nrxx];
	// cout<<"npw:"<<GlobalC::wf.npw<<endl;
	// cout<<"npwx:"<<GlobalC::wf.npwx<<endl;
	for (int ig = 0; ig < GlobalC::wf.npw; ig++)
	{
		// GR_index[ig] = GlobalC::wfcpw->ig2fftw[GlobalC::wf.igk(ik, ig)];
		// if(ig<20){
		// 	cout<<GR_index[ig]<<" "<< GlobalC::wf.igk(ik, ig) <<" "<<GlobalC::wfcpw->ig2fftw[ GlobalC::wf.igk(ik, ig)
		// ]<<endl;
		// }
	}
	// cout<<"init_K"<<endl;
#ifdef __ROCM
	CHECK_CUDA(hipMalloc((void **)&GR_index_d, GlobalC::wf.npwx * sizeof(int)));
	CHECK_CUDA(hipMemcpy(GR_index_d, GR_index, GlobalC::wf.npwx * sizeof(int), hipMemcpyHostToDevice));
	// delete [] GR_index_tmp;
#endif

	// (7) ik
	GlobalV::CURRENT_K = ik;

	return;
}

//----------------------------------------------------------------------
// Hamiltonian diagonalization in the subspace spanned
// by nstart states psi (atomic or random wavefunctions).
// Produces on output n_band eigenvectors (n_band <= nstart) in evc.
//----------------------------------------------------------------------
void Hamilt_PW::diagH_subspace(const int ik,
							   const int nstart,
							   const int n_band,
							   const ModuleBase::ComplexMatrix &psi,
							   ModuleBase::ComplexMatrix &evc,
							   double *en)
{
	ModuleBase::TITLE("Hamilt_PW_HIP", "diagH_subspace");
	ModuleBase::timer::tick("Hamilt_PW_HIP", "diagH_subspace");

	assert(nstart != 0);
	assert(n_band != 0);

	ModuleBase::ComplexMatrix hc(nstart, nstart);
	ModuleBase::ComplexMatrix sc(nstart, nstart);
	ModuleBase::ComplexMatrix hvec(nstart, n_band);

	int dmin = 0;
	int dmax = 0;
	const int npw = GlobalC::kv.ngk[ik];

	if (GlobalV::NSPIN != 4)
	{
		dmin = npw;
		dmax = GlobalC::wf.npwx;
	}
	else
	{
		dmin = GlobalC::wf.npwx * GlobalV::NPOL;
		dmax = GlobalC::wf.npwx * GlobalV::NPOL;
	}

	// qianrui improve this part 2021-3-14
	std::complex<double> *aux = new std::complex<double>[dmax * nstart];
	// std::complex<double> *paux = aux;
	// std::complex<double> *ppsi = psi.c;

	// qianrui replace it
	this->h_psi(psi.c, aux, nstart);

	char trans1 = 'C';
	char trans2 = 'N';
	zgemm_(&trans1,
		   &trans2,
		   &nstart,
		   &nstart,
		   &dmin,
		   &ModuleBase::ONE,
		   psi.c,
		   &dmax,
		   aux,
		   &dmax,
		   &ModuleBase::ZERO,
		   hc.c,
		   &nstart);
	hc = transpose(hc, false);

	zgemm_(&trans1,
		   &trans2,
		   &nstart,
		   &nstart,
		   &dmin,
		   &ModuleBase::ONE,
		   psi.c,
		   &dmax,
		   psi.c,
		   &dmax,
		   &ModuleBase::ZERO,
		   sc.c,
		   &nstart);
	sc = transpose(sc, false);

	delete[] aux;

	// Peize Lin add 2019-03-09
#ifdef __LCAO
	if (GlobalV::BASIS_TYPE == "lcao_in_pw")
	{
		auto add_Hexx = [&](const double alpha) {
			for (int m = 0; m < nstart; ++m)
			{
				for (int n = 0; n < nstart; ++n)
				{
					hc(m, n) += alpha * GlobalC::exx_lip.get_exx_matrix()[ik][m][n];
				}
			}
		};
		if(XC_Functional::get_func_type()==4)
		{
			if ( Exx_Global::Hybrid_Type::HF   == GlobalC::exx_lcao.info.hybrid_type ) // HF
			{
				add_Hexx(1);
			}
			else if (Exx_Global::Hybrid_Type::PBE0 == GlobalC::exx_lcao.info.hybrid_type || 
					Exx_Global::Hybrid_Type::HSE  == GlobalC::exx_lcao.info.hybrid_type) // PBE0 or HSE
			{
				add_Hexx(GlobalC::exx_global.info.hybrid_alpha);
			}
		}
	}
#endif

	if (GlobalV::NPROC_IN_POOL > 1)
	{
		Parallel_Reduce::reduce_complex_double_pool(hc.c, nstart * nstart);
		Parallel_Reduce::reduce_complex_double_pool(sc.c, nstart * nstart);
	}

	// after generation of H and S matrix, diag them
	GlobalC::hm.diagH_LAPACK(nstart, n_band, hc, sc, nstart, en, hvec);

	// Peize Lin add 2019-03-09
#ifdef __LCAO
	if ("lcao_in_pw" == GlobalV::BASIS_TYPE)
	{
		switch (GlobalC::exx_global.info.hybrid_type)
		{
		case Exx_Global::Hybrid_Type::HF:
		case Exx_Global::Hybrid_Type::PBE0:
		case Exx_Global::Hybrid_Type::HSE:
			GlobalC::exx_lip.k_pack->hvec_array[ik] = hvec;
			break;
		}
	}
#endif

	//=======================
	// diagonize the H-matrix
	//=======================

	// for tests
	/*
			std::cout << std::setprecision(3);
			out.printV3(GlobalV::ofs_running,GlobalC::kv.kvec_c[ik]);
			out.printcm_norm("sc",sc,1.0e-4);
			out.printcm_norm("hvec",hvec,1.0e-4);
			out.printcm_norm("hc",hc,1.0e-4);
			std::cout << std::endl;
	*/

	std::cout << std::setprecision(5);

	//--------------------------
	// KEEP THIS BLOCK FOR TESTS
	//--------------------------
	/*
		std::cout << "  hc matrix" << std::endl;
		for(int i=0; i<GlobalV::NLOCAL; i++)
		{
			for(int j=0; j<GlobalV::NLOCAL; j++)
			{
				double a = hc(i,j).real();
				if(abs(a) < 1.0e-5) a = 0;
				std::cout << std::setw(6) << a;
			}
			std::cout << std::endl;
		}

		std::cout << "  sc matrix" << std::endl;
		for(int i=0; i<GlobalV::NLOCAL; i++)
		{
			for(int j=0; j<GlobalV::NLOCAL; j++)
			{
				double a = sc(i,j).real();
				if(abs(a) < 1.0e-5) a = 0;
				std::cout << std::setw(6) << a;
			}
			std::cout << std::endl;
		}

		std::cout << "\n Band Energy" << std::endl;
		for(int i=0; i<GlobalV::NBANDS; i++)
		{
			std::cout << " e[" << i+1 << "]=" << en[i] * Ry_to_eV << std::endl;
		}
	*/
	//--------------------------
	// KEEP THIS BLOCK FOR TESTS
	//--------------------------

	if ((GlobalV::BASIS_TYPE == "lcao" || GlobalV::BASIS_TYPE == "lcao_in_pw") && GlobalV::CALCULATION == "nscf"
		&& !Optical::opt_epsilon2)
	{
		GlobalV::ofs_running << " Not do zgemm to get evc." << std::endl;
	}
	else if ((GlobalV::BASIS_TYPE == "lcao" || GlobalV::BASIS_TYPE == "lcao_in_pw")
			 && (GlobalV::CALCULATION == "scf" || GlobalV::CALCULATION == "md"
				 || GlobalV::CALCULATION == "relax")) // pengfei 2014-10-13
	{
		// because psi and evc are different here,
		// I think if psi and evc are the same,
		// there may be problems, mohan 2011-01-01
		char transa = 'N';
		char transb = 'T';
		zgemm_(&transa,
			   &transb,
			   &dmax, // m: row of A,C
			   &n_band, // n: col of B,C
			   &nstart, // k: col of A, row of B
			   &ModuleBase::ONE, // alpha
			   psi.c, // A
			   &dmax, // LDA: if(N) max(1,m) if(T) max(1,k)
			   hvec.c, // B
			   &n_band, // LDB: if(N) max(1,k) if(T) max(1,n)
			   &ModuleBase::ZERO, // belta
			   evc.c, // C
			   &dmax); // LDC: if(N) max(1, m)
	}
	else
	{
		// As the evc and psi may refer to the same matrix, we first
		// create a temporary matrix to story the result. (by wangjp)
		// qianrui improve this part 2021-3-13
		char transa = 'N';
		char transb = 'T';
		ModuleBase::ComplexMatrix evctmp(n_band, dmin, false);
		zgemm_(&transa,
			   &transb,
			   &dmin,
			   &n_band,
			   &nstart,
			   &ModuleBase::ONE,
			   psi.c,
			   &dmax,
			   hvec.c,
			   &n_band,
			   &ModuleBase::ZERO,
			   evctmp.c,
			   &dmin);
		for (int ib = 0; ib < n_band; ib++)
		{
			for (int ig = 0; ig < dmin; ig++)
			{
				evc(ib, ig) = evctmp(ib, ig);
			}
		}
	}
	// out.printr1_d("en",en,n_band);

	//	std::cout << "\n bands" << std::endl;
	//	for(int ib=0; ib<n_band; ib++)
	//	{
	//		std::cout << " ib=" << ib << " " << en[ib] * Ry_to_eV << std::endl;
	//	}

	// out.printcm_norm("hvec",hvec,1.0e-8);

	ModuleBase::timer::tick("Hamilt_PW_HIP", "diagH_subspace");
	return;
}

void Hamilt_PW::diagH_subspace_cuda(const int ik,
									const int nstart,
									const int n_band,
									const hipblasDoubleComplex *psi_c, // matrix
									hipblasDoubleComplex *evc, // matrix
									double *en,
									hipblasDoubleComplex *vkb_c)
{
	ModuleBase::TITLE("Hamilt_PW_HIP", "diagH_subspace_hip");
	ModuleBase::timer::tick("Hamilt_PW_HIP", "diagH_subspace_hip");

	assert(nstart != 0);
	assert(n_band != 0);

	hipblasDoubleComplex *hc;
	hipblasDoubleComplex *tmp_hc;
	hipblasDoubleComplex *sc;
	hipblasDoubleComplex *hvec;
	CHECK_CUDA(hipMalloc((void **)&hc, nstart * nstart * sizeof(hipblasDoubleComplex)));
	CHECK_CUDA(hipMalloc((void **)&tmp_hc, nstart * nstart * sizeof(hipblasDoubleComplex)));
	CHECK_CUDA(hipMalloc((void **)&sc, nstart * nstart * sizeof(hipblasDoubleComplex)));
	CHECK_CUDA(hipMalloc((void **)&hvec, nstart * n_band * sizeof(hipblasDoubleComplex)));

	int dmin = 0;
	int dmax = 0;
	const int npw = GlobalC::kv.ngk[ik];

	// cout<<"npw:"<<npw<<endl;
	// cout<<"npwx"<<GlobalC::wf.npwx<<endl;

	if (GlobalV::NSPIN != 4)
	{
		dmin = npw;
		dmax = GlobalC::wf.npwx;
	}
	else
	{
		dmin = GlobalC::wf.npwx * GlobalV::NPOL;
		dmax = GlobalC::wf.npwx * GlobalV::NPOL;
	}

	// //qianrui improve this part 2021-3-14
	// std::complex<double> *aux=new std::complex<double> [dmax*nstart];
	// std::complex<double> *paux = aux;
	// std::complex<double> *ppsi = psi.c;

	hipblasDoubleComplex *aux;
	CHECK_CUDA(hipMalloc((void **)&aux, dmax * nstart * sizeof(hipblasDoubleComplex)));
	// hipblasDoubleComplex* paux = aux;
	// const hipblasDoubleComplex* ppsi = psi_c; // ?

	// qianrui replace it
	// cout<<"psi before hpsi"<<endl;
	// print_test<hipblasDoubleComplex>((hipblasDoubleComplex*)psi_c, 10);
	this->h_psi_cuda(psi_c, aux, vkb_c, nstart); // TODO: vkb_c  nstart!=1 ?
	// cout<<"aux after hpsi"<<endl;
	// print_test<hipblasDoubleComplex>(aux, 10);

	// double2 ONE, ZERO;
	// ONE.y = ZERO.x = ZERO.y = 0.0;
	// ONE.x = 1.0;

	hipblasDoubleComplex ONE(1, 0);
	hipblasDoubleComplex ZERO(0, 0);
	hipblasDoubleComplex NEG_ONE(-1, 0);

	hipblasOperation_t trans1 = HIPBLAS_OP_C;
	hipblasOperation_t trans2 = HIPBLAS_OP_N;
	CHECK_CUBLAS(hipblasZgemm(hpw_handle,
							  trans1,
							  trans2,
							  nstart,
							  nstart,
							  dmin,
							  &ONE,
							  psi_c,
							  dmax,
							  aux,
							  dmax,
							  &ZERO,
							  tmp_hc,
							  nstart));
	// hc=transpose(hc,false); // TODO: transpose

	// use 'geam' API todo transpose.
	// double2 t_alpha, t_beta;
	// t_alpha.y = t_beta.x = t_beta.y = 0.0;
	// t_alpha.x = 1.0;

	hipblasDoubleComplex t_alpha(1, 0);
	hipblasDoubleComplex t_beta(0, 0);

	CHECK_CUBLAS(hipblasZgeam(hpw_handle,
							  HIPBLAS_OP_T,
							  HIPBLAS_OP_T,
							  nstart,
							  nstart,
							  &t_alpha,
							  tmp_hc,
							  nstart,
							  &t_beta,
							  tmp_hc,
							  nstart,
							  hc,
							  nstart));

	CHECK_CUBLAS(hipblasZgemm(hpw_handle,
							  trans1,
							  trans2,
							  nstart,
							  nstart,
							  dmin,
							  &ONE,
							  psi_c,
							  dmax,
							  psi_c,
							  dmax,
							  &ZERO,
							  tmp_hc,
							  nstart));
	// sc=transpose(sc,false); // TODO: transpose
	CHECK_CUBLAS(hipblasZgeam(hpw_handle,
							  HIPBLAS_OP_T,
							  HIPBLAS_OP_T,
							  nstart,
							  nstart,
							  &t_alpha,
							  tmp_hc,
							  nstart,
							  &t_beta,
							  tmp_hc,
							  nstart,
							  sc,
							  nstart));

	CHECK_CUDA(hipFree(aux));

	// if(GlobalV::NPROC_IN_POOL>1)
	// {
	// 	Parallel_Reduce::reduce_complex_double_pool( hc.c, nstart*nstart );
	// 	Parallel_Reduce::reduce_complex_double_pool( sc.c, nstart*nstart );
	// }

	// after generation of H and S matrix, diag them

	// Method1 : Do with diagH_LAPACK

	ModuleBase::ComplexMatrix h_hc(nstart, nstart);
	ModuleBase::ComplexMatrix h_sc(nstart, nstart);
	ModuleBase::ComplexMatrix h_hvec(nstart, n_band);

	double *h_en = new double[n_band];

	CHECK_CUDA(hipMemcpy(h_hc.c, hc, nstart * nstart * sizeof(hipblasDoubleComplex), hipMemcpyDeviceToHost));
	CHECK_CUDA(hipMemcpy(h_sc.c, sc, nstart * nstart * sizeof(hipblasDoubleComplex), hipMemcpyDeviceToHost));

	GlobalC::hm.diagH_LAPACK(nstart, n_band, h_hc, h_sc, nstart, h_en, h_hvec);
	CHECK_CUDA(hipMemcpy(hvec, h_hvec.c, nstart * n_band * sizeof(hipblasDoubleComplex), hipMemcpyHostToDevice));
	CHECK_CUDA(hipMemcpy(en, h_en, n_band * sizeof(double), hipMemcpyHostToDevice));
	delete[] h_en;

	// Method2 : Do with diagH_CUSOLVER
	// GlobalC::hm.diagH_CUSOLVER(nstart, n_band, hc, sc, nstart, en, hvec);

	// Peize Lin add 2019-03-09
	/*
	#ifdef __LCAO
		if("lcao_in_pw"==GlobalV::BASIS_TYPE)
		{
			switch(GlobalC::exx_global.info.hybrid_type)
			{
				case Exx_Global::Hybrid_Type::HF:
				case Exx_Global::Hybrid_Type::PBE0:
				case Exx_Global::Hybrid_Type::HSE:
					GlobalC::exx_lip.k_pack->hvec_array[ik] = hvec;
					break;
			}
		}
	#endif
	*/

	std::cout << std::setprecision(5);

	if ((GlobalV::BASIS_TYPE == "lcao" || GlobalV::BASIS_TYPE == "lcao_in_pw") && GlobalV::CALCULATION == "nscf"
		&& !Optical::opt_epsilon2)
	{
		GlobalV::ofs_running << " Not do zgemm to get evc." << std::endl;
	}
	else if ((GlobalV::BASIS_TYPE == "lcao" || GlobalV::BASIS_TYPE == "lcao_in_pw")
			 && (GlobalV::CALCULATION == "scf" || GlobalV::CALCULATION == "md"
				 || GlobalV::CALCULATION == "relax")) // pengfei 2014-10-13
	{
		// because psi and evc are different here,
		// I think if psi and evc are the same,
		// there may be problems, mohan 2011-01-01
		// char transa = 'N';
		// char transb = 'T';
		// zgemm_( &transa,
		// 		&transb,
		// 		&dmax, // m: row of A,C
		// 		&n_band, // n: col of B,C
		// 		&nstart, // k: col of A, row of B
		// 		&ModuleBase::ONE, // alpha
		// 		psi.c, // A
		// 		&dmax, // LDA: if(N) max(1,m) if(T) max(1,k)
		// 		hvec.c, // B
		// 		&n_band, // LDB: if(N) max(1,k) if(T) max(1,n)
		// 		&ModuleBase::ZERO,  // belta
		// 		evc.c, // C
		// 		&dmax ); // LDC: if(N) max(1, m)
		hipblasOperation_t transa = HIPBLAS_OP_N;
		hipblasOperation_t transb = HIPBLAS_OP_T;
		// double2 ONE, ZERO;
		// ONE.y = ZERO.x = ZERO.y = 0.0;
		// ONE.x = 1.0;

		hipblasDoubleComplex ONE(1, 0);
		hipblasDoubleComplex ZERO(0, 0);
		hipblasDoubleComplex NEG_ONE(-1, 0);

		CHECK_CUBLAS(hipblasZgemm(hpw_handle,
								  transa,
								  transb,
								  dmax,
								  n_band,
								  nstart,
								  &ONE,
								  psi_c,
								  dmax,
								  hvec,
								  n_band,
								  &ZERO,
								  evc,
								  dmax));
	}
	else
	{
		// As the evc and psi may refer to the same matrix, we first
		// create a temporary matrix to story the result. (by wangjp)
		// qianrui improve this part 2021-3-13
		// char transa = 'N';
		// char transb = 'T';
		hipblasOperation_t transa = HIPBLAS_OP_N;
		hipblasOperation_t transb = HIPBLAS_OP_T;

		// double2 ONE, ZERO;
		// ONE.y = ZERO.x = ZERO.y = 0.0;
		// ONE.x = 1.0;
		hipblasDoubleComplex ONE(1, 0);
		hipblasDoubleComplex ZERO(0, 0);
		hipblasDoubleComplex NEG_ONE(-1, 0);

		hipblasDoubleComplex *evctmp;
		CHECK_CUDA(hipMalloc((void **)&evctmp, n_band * dmin * sizeof(hipblasDoubleComplex)));
		CHECK_CUBLAS(hipblasZgemm(hpw_handle,
								  transa,
								  transb,
								  dmin,
								  n_band,
								  nstart,
								  &ONE,
								  psi_c,
								  dmax,
								  hvec,
								  n_band,
								  &ZERO,
								  evctmp,
								  dmin));
		// cout<<"evctmp before cpy back"<<endl;
		// print_test<hipblasDoubleComplex>(evctmp, 15);
		for (int ib = 0; ib < n_band; ib++)
		{
			// for(int ig=0; ig<dmin; ig++)
			// {
			// 	evc(ib,ig) = evctmp(ib,ig);
			// }
			CHECK_CUDA(hipMemcpy(&evc[ib * dmax],
								 &evctmp[ib * dmin],
								 dmin * sizeof(hipblasDoubleComplex),
								 hipMemcpyDeviceToDevice));
		}
		CHECK_CUDA(hipFree(evctmp));
		hipDeviceSynchronize();
	}

	CHECK_CUDA(hipFree(hc));
	CHECK_CUDA(hipFree(tmp_hc));
	CHECK_CUDA(hipFree(sc));
	CHECK_CUDA(hipFree(hvec));

	ModuleBase::timer::tick("Hamilt_PW_HIP", "diagH_subspace_hip");
	return;
}

void Hamilt_PW::h_1psi_cuda(const int npw_in,
							const hipblasComplex *psi,
							hipblasComplex *hpsi,
							hipblasComplex *spsi,
							hipblasComplex *vkb_c)
{
	this->h_psi_cuda(psi, hpsi, vkb_c);

	int thread = 512;
	int block = (npw_in + thread - 1) / thread;
	hipLaunchKernelGGL(kernel_copy,
					   dim3(thread),
					   dim3(block),
					   0,
					   0,
					   npw_in,
					   reinterpret_cast<float2 *>(spsi),
					   reinterpret_cast<const float2 *>(psi));
	return;
}

void Hamilt_PW::h_1psi_cuda(const int npw_in,
							const hipblasDoubleComplex *psi,
							hipblasDoubleComplex *hpsi,
							hipblasDoubleComplex *spsi,
							hipblasDoubleComplex *vkb_c)
{
	this->h_psi_cuda(psi, hpsi, vkb_c);

	int thread = 512;
	int block = (npw_in + thread - 1) / thread;
	hipLaunchKernelGGL(kernel_copy,
					   dim3(thread),
					   dim3(block),
					   0,
					   0,
					   npw_in,
					   reinterpret_cast<double2 *>(spsi),
					   reinterpret_cast<const double2 *>(psi));
	return;
}

void Hamilt_PW::s_1psi_cuda(const int dim, const hipblasComplex *psi, hipblasComplex *spsi)
{
	CHECK_CUDA(hipMemcpy(spsi, psi, dim * sizeof(hipblasComplex), hipMemcpyDeviceToDevice));
	return;
}

void Hamilt_PW::s_1psi_cuda(const int dim, const hipblasDoubleComplex *psi, hipblasDoubleComplex *spsi)
{
	CHECK_CUDA(hipMemcpy(spsi, psi, dim * sizeof(hipblasDoubleComplex), hipMemcpyDeviceToDevice));
	return;
}

void Hamilt_PW::h_1psi(const int npw_in,
					   const std::complex<double> *psi,
					   std::complex<double> *hpsi,
					   std::complex<double> *spsi)
{
	this->h_psi(psi, hpsi);

	for (int i = 0; i < npw_in; i++)
	{
		spsi[i] = psi[i];
	}
	return;
}

void Hamilt_PW::s_1psi(const int dim, const std::complex<double> *psi, std::complex<double> *spsi)
{
	for (int i = 0; i < dim; i++)
	{
		spsi[i] = psi[i];
	}
	return;
}

void Hamilt_PW::h_psi_cuda(const hipblasComplex *psi_in, hipblasComplex *hpsi, hipblasComplex *vkb_c, const int m)
{
	ModuleBase::timer::tick("Hamilt_PW_HIP", "h_psi");
	// int i = 0;
	// int j = 0;
	// int ig= 0;

	// if(NSPIN!=4) ZEROS(hpsi, wf.npw);
	// else ZEROS(hpsi, wf.npwx * NPOL);//added by zhengdy-soc
	int dmax = GlobalC::wf.npwx * GlobalV::NPOL;

	// cout<<"dim inside="<<GlobalC::wf.npwx * GlobalV::NPOL<<endl;

	//------------------------------------
	//(1) the kinetical energy.
	//------------------------------------
	hipblasComplex *tmhpsi;
	const hipblasComplex *tmpsi_in;
	ModuleBase::timer::tick("Hamilt_PW_HIP", "kinetic");
	if (GlobalV::T_IN_H)
	{
		tmhpsi = hpsi;
		tmpsi_in = psi_in;

		double *d_g2kin;
		float *f_g2kin;

		CHECK_CUDA(hipMalloc((void **)&d_g2kin, GlobalC::wf.npwx * sizeof(double)));
		CHECK_CUDA(hipMalloc((void **)&f_g2kin, GlobalC::wf.npwx * sizeof(float)));
		// CHECK_CUDA(hipMemcpy(d_g2kin, GlobalC::wf.g2kin, GlobalC::wf.npw * sizeof(double), hipMemcpyHostToDevice));

		int thread = 512;
		int block = (GlobalC::wf.npw + thread - 1) / thread;
		hipLaunchKernelGGL(cast_d2f, dim3(block), dim3(thread), 0, 0, f_g2kin, d_g2kin, GlobalC::wf.npw);

		for (int ib = 0; ib < m; ++ib)
		{
			// cout<<"in hpsi-Kinetic, iband = "<<ib<<endl;

			// todo: template
			hipLaunchKernelGGL(HIP_KERNEL_NAME(kernel_get_tmhpsi<float, float2>),
							   dim3(block),
							   dim3(thread),
							   0,
							   0,
							   GlobalC::wf.npw,
							   reinterpret_cast<float2 *>(tmhpsi),
							   reinterpret_cast<const float2 *>(tmpsi_in),
							   f_g2kin);

			// if(GlobalC::NSPIN==4){
			//     for(ig=GlobalC::wf.npw; ig < GlobalC::wf.npwx; ++ig)
			//     {
			//         tmhpsi[ig] = 0;
			//     }
			//     tmhpsi += GlobalC::wf.npwx;
			//     tmpsi_in += GlobalC::wf.npwx;
			//     for (ig = 0;ig < GlobalC::wf.npw ;++ig)
			//     {
			//         tmhpsi[ig] = GlobalC::wf.g2kin[ig] * tmpsi_in[ig];
			//     }
			//     // TODO: setup with 0
			//     for(ig=GlobalC::wf.npw; ig < GlobalC::wf.npwx; ++ig)
			//     {
			//         tmhpsi[ig] =0;
			//     }
			// }

			tmhpsi += GlobalC::wf.npwx;
			tmpsi_in += GlobalC::wf.npwx;
		}
		CHECK_CUDA(hipFree(d_g2kin));
		CHECK_CUDA(hipFree(f_g2kin));
	}

	ModuleBase::timer::tick("Hamilt_PW_HIP", "kinetic");

	// cout<<"======after hpsi part I======="<<endl;
	// print_test<hipblasComplex>(hpsi, 15);

	//------------------------------------
	//(2) the local potential.
	//-----------------------------------
	ModuleBase::timer::tick("Hamilt_PW_HIP", "vloc");
	//  ...
	if (GlobalV::VL_IN_H)
	{
		tmhpsi = hpsi;
		tmpsi_in = psi_in;
		// int *d_GR_index;
		double *d_vr_eff1;
		float *f_vr_eff1;
		hipblasComplex *f_porter;

		// CHECK_CUDA(hipMalloc((void**)&d_GR_index, GlobalC::wf.npwx * sizeof(int)));
		CHECK_CUDA(hipMalloc((void **)&d_vr_eff1, GlobalC::wfcpw->nrxx * sizeof(double)));
		CHECK_CUDA(hipMalloc((void **)&f_vr_eff1, GlobalC::wfcpw->nrxx * sizeof(float)));
		CHECK_CUDA(hipMalloc((void **)&f_porter, GlobalC::wfcpw->nrxx * sizeof(hipblasComplex)));

		CHECK_CUDA(
			hipMemcpy(d_vr_eff1, GlobalC::pot.vr_eff1, GlobalC::wfcpw->nrxx * sizeof(double), hipMemcpyHostToDevice));

		int thread2 = 512;
		int block2 = (GlobalC::wfcpw->nrxx + thread2 - 1) / thread2;
		hipLaunchKernelGGL(cast_d2f, dim3(block2), dim3(thread2), 0, 0, f_vr_eff1, d_vr_eff1, GlobalC::wfcpw->nrxx);

		// cout<<"NSPIN = "<<GlobalV::NSPIN<<endl;
		for (int ib = 0; ib < m; ++ib)
		{
			// cout<<"in hpsi:loacl_pot, iband = "<<ib<<endl;
			// if(NSPIN!=4){
			// ZEROS( UFFT.porter, pw.nrxx);
			CHECK_CUDA(hipMemset(f_porter, 0, GlobalC::wfcpw->nrxx * sizeof(hipblasComplex)));

			// todo
			// GlobalC::UFFT.RoundTrip(tmpsi_in, f_vr_eff1, GR_index_d, f_porter);
			GlobalC::wfcpw->recip2real(tmpsi_in, f_porter, ik);
			for (int ir=0; ir< nrxx; ir++)
			{
				f_porter[ir] *=  GlobalC::pot.vr_eff1[ir];
			}
			GlobalC::wfcpw->real2recip(f_porter, tmhpsi, ik, true);

			// for (j = 0;j < wf.npw;j++)
			// {
			//     tmhpsi[j] += UFFT.porter[ GR_index[j] ];
			// }
			int thread = 512;
			int block = (GlobalC::wf.npw + thread - 1) / thread;
			hipLaunchKernelGGL(HIP_KERNEL_NAME(kernel_add_tmhpsi<float2>),
							   dim3(block),
							   dim3(thread),
							   0,
							   0,
							   GlobalC::wf.npw,
							   reinterpret_cast<float2 *>(tmhpsi),
							   reinterpret_cast<float2 *>(f_porter),
							   GR_index_d);

			tmhpsi += dmax;
			tmpsi_in += dmax;
		}
		// CHECK_CUDA(hipFree(d_GR_index));
		CHECK_CUDA(hipFree(d_vr_eff1));
		CHECK_CUDA(hipFree(f_vr_eff1));
		CHECK_CUDA(hipFree(f_porter));
	}
	ModuleBase::timer::tick("Hamilt_PW_HIP", "vloc");
	//------------------------------------
	// (3) the nonlocal pseudopotential.
	//------------------------------------
	ModuleBase::timer::tick("Hamilt_PW_HIP", "vnl");

	// cout<<"======after hpsi part II======="<<endl;
	// print_test<hipblasComplex>(hpsi, 15);

	if (GlobalV::VNL_IN_H)
	{
		if (GlobalC::ppcell.nkb > 0)
		{
			int nkb = GlobalC::ppcell.nkb;
			hipblasComplex *becp;
			// hipblasDoubleComplex *d_vkb_c;
			// hipblasComplex *f_vkb_c;
			CHECK_CUDA(hipMalloc((void **)&becp, GlobalV::NPOL * m * nkb * sizeof(hipblasComplex)));
			// CHECK_CUDA(hipMalloc((void**)&d_vkb_c, GlobalC::wf.npwx*nkb*sizeof(hipblasDoubleComplex)));
			// CHECK_CUDA(hipMalloc((void**)&f_vkb_c, GlobalC::wf.npwx*nkb*sizeof(hipblasComplex)));

			// CHECK_CUDA(hipMemcpy(d_vkb_c, GlobalC::ppcell.vkb.c, GlobalC::wf.npwx*nkb*sizeof(hipblasDoubleComplex),
			// hipMemcpyHostToDevice)); int thread = 512; int block = (GlobalC::wf.npwx*nkb + thread - 1) / thread;
			// hipLaunchKernelGGL(cast_d2f, dim3(block), dim3(thread), 0, 0, f_vkb_c, d_vkb_c, GlobalC::wf.npwx*nkb);

			hipblasOperation_t transa = HIPBLAS_OP_C;
			hipblasOperation_t transb = HIPBLAS_OP_N;
			// hipblasHandle_t handle;
			// CHECK_CUBLAS(hipblasCreate(&handle));

			// float2 ONE, ZERO;
			// ONE.y = ZERO.x = ZERO.y = 0.0;
			// ONE.x = 1.0;

			hipblasComplex ZERO(0, 0);
			hipblasComplex ONE(1, 0);
			// NEG_ONE.x = -1.0;

			// cout<<"===== vkbc ===="<<endl;
			// print_test<hipblasComplex>(f_vkb_c, 15);

			if (m == 1 && GlobalV::NPOL == 1)
			{
				int inc = 1;
				CHECK_CUBLAS(hipblasCgemv(hpw_handle,
										  transa,
										  GlobalC::wf.npw,
										  nkb,
										  &ONE,
										  vkb_c,
										  GlobalC::wf.npwx,
										  psi_in,
										  inc,
										  &ZERO,
										  becp,
										  inc));
			}
			else
			{
				int npm = GlobalV::NPOL * m;
				CHECK_CUBLAS(hipblasCgemm(hpw_handle,
										  transa,
										  transb,
										  nkb,
										  npm,
										  GlobalC::wf.npw,
										  &ONE,
										  vkb_c,
										  GlobalC::wf.npwx,
										  psi_in,
										  GlobalC::wf.npwx,
										  &ZERO,
										  becp,
										  nkb));
			}

			// complex<double> *hpsi_cpu = new complex<double>[GlobalC::wf.npw*GlobalV::NPOL];
			// complex<double> *becp_cpu = new complex<double>[GlobalV::NPOL*m*nkb];

			// CHECK_CUDA(hipMemcpy(becp_cpu, becp, GlobalV::NPOL*m*nkb*sizeof(hipblasComplex), hipMemcpyDeviceToHost));

			// CHECK_CUDA(hipMemcpy(hpsi_cpu, hpsi, GlobalC::wf.npw*GlobalV::NPOL*sizeof(hipblasComplex),
			// hipMemcpyDeviceToHost));

			// this->add_nonlocal_pp(hpsi_cpu, becp_cpu, m);

			// CHECK_CUDA(hipMemcpy(hpsi, hpsi_cpu, GlobalC::wf.npw*GlobalV::NPOL*sizeof(hipblasComplex),
			// hipMemcpyHostToDevice));

			// delete [] hpsi_cpu;
			// delete [] becp_cpu;

			// cout<<"===== becp before add nonloaclpp ===="<<endl;
			// print_test<hipblasComplex>(becp, 15);

			this->add_nonlocal_pp_cuda(hpsi, becp, vkb_c, m);

			// CHECK_CUBLAS(hipblasDestroy(handle));
			CHECK_CUDA(hipFree(becp));
			// CHECK_CUDA(hipFree(d_vkb_c));
			// CHECK_CUDA(hipFree(f_vkb_c));
			// cout<<"nonlocal end"<<endl;
		}
	}

	ModuleBase::timer::tick("Hamilt_PW_HIP", "vnl");

	// cout<<"======after hpsi part III======="<<endl;
	// print_test<hipblasComplex>(hpsi, 15);

	//------------------------------------
	// (4) the metaGGA part
	//------------------------------------
	// TODO: add metaGGA part

	ModuleBase::timer::tick("Hamilt_PW_HIP", "h_psi");
	return;
}

void Hamilt_PW::h_psi_cuda(const hipblasDoubleComplex *psi_in,
						   hipblasDoubleComplex *hpsi,
						   hipblasDoubleComplex *vkb_c,
						   const int m)
{
	ModuleBase::timer::tick("Hamilt_PW_HIP", "h_psi");
	// int i = 0;
	// int j = 0;
	// int ig= 0;

	// if(NSPIN!=4) ZEROS(hpsi, wf.npw);
	// else ZEROS(hpsi, wf.npwx * NPOL);//added by zhengdy-soc
	int dmax = GlobalC::wf.npwx * GlobalV::NPOL;

	// cout<<"dim inside="<<GlobalC::wf.npwx * GlobalV::NPOL<<endl;

	//------------------------------------
	//(1) the kinetical energy.
	//------------------------------------
	hipblasDoubleComplex *tmhpsi;
	const hipblasDoubleComplex *tmpsi_in;
	ModuleBase::timer::tick("Hamilt_PW_HIP", "kinetic");
	if (GlobalV::T_IN_H)
	{
		tmhpsi = hpsi;
		tmpsi_in = psi_in;

		double *d_g2kin;
		CHECK_CUDA(hipMalloc((void **)&d_g2kin, GlobalC::wf.npwx * sizeof(double)));
		// CHECK_CUDA(hipMemcpy(d_g2kin, GlobalC::wf.g2kin, GlobalC::wf.npw * sizeof(double), hipMemcpyHostToDevice));
		for (int ib = 0; ib < m; ++ib)
		{
			// cout<<"in hpsi-Kinetic, iband = "<<ib<<endl;

			int thread = 512;
			int block = (GlobalC::wf.npw + thread - 1) / thread;
			hipLaunchKernelGGL(HIP_KERNEL_NAME(kernel_get_tmhpsi<double, double2>),
							   dim3(block),
							   dim3(thread),
							   0,
							   0,
							   GlobalC::wf.npw,
							   reinterpret_cast<double2 *>(tmhpsi),
							   reinterpret_cast<const double2 *>(tmpsi_in),
							   d_g2kin);

			// if(GlobalC::NSPIN==4){
			//     for(ig=GlobalC::wf.npw; ig < GlobalC::wf.npwx; ++ig)
			//     {
			//         tmhpsi[ig] = 0;
			//     }
			//     tmhpsi += GlobalC::wf.npwx;
			//     tmpsi_in += GlobalC::wf.npwx;
			//     for (ig = 0;ig < GlobalC::wf.npw ;++ig)
			//     {
			//         tmhpsi[ig] = GlobalC::wf.g2kin[ig] * tmpsi_in[ig];
			//     }
			//     // TODO: setup with 0
			//     for(ig=GlobalC::wf.npw; ig < GlobalC::wf.npwx; ++ig)
			//     {
			//         tmhpsi[ig] =0;
			//     }
			// }

			tmhpsi += GlobalC::wf.npwx;
			tmpsi_in += GlobalC::wf.npwx;
		}
		CHECK_CUDA(hipFree(d_g2kin));
	}

	ModuleBase::timer::tick("Hamilt_PW_HIP", "kinetic");

	// cout<<"======after hpsi part I======="<<endl;
	// print_test<hipblasDoubleComplex>(hpsi, 10);

	//------------------------------------
	//(2) the local potential.
	//-----------------------------------
	ModuleBase::timer::tick("Hamilt_PW_HIP", "vloc");
	//  ...
	if (GlobalV::VL_IN_H)
	{
		tmhpsi = hpsi;
		tmpsi_in = psi_in;
		// int *d_GR_index;
		double *d_vr_eff1;
		hipblasDoubleComplex *d_porter;

		// CHECK_CUDA(hipMalloc((void**)&d_GR_index, GlobalC::wf.npwx * sizeof(int)));
		CHECK_CUDA(hipMalloc((void **)&d_vr_eff1, GlobalC::wfcpw->nrxx * sizeof(double)));
		CHECK_CUDA(hipMalloc((void **)&d_porter, GlobalC::wfcpw->nrxx * sizeof(hipblasDoubleComplex)));

		CHECK_CUDA(
			hipMemcpy(d_vr_eff1, GlobalC::pot.vr_eff1, GlobalC::wfcpw->nrxx * sizeof(double), hipMemcpyHostToDevice));
		// cout<<"NSPIN = "<<GlobalV::NSPIN<<endl;
		for (int ib = 0; ib < m; ++ib)
		{
			// cout<<"in hpsi:loacl_pot, iband = "<<ib<<endl;
			// if(NSPIN!=4){
			// ZEROS( UFFT.porter, pw.nrxx);
			CHECK_CUDA(hipMemset(d_porter, 0, GlobalC::wfcpw->nrxx * sizeof(hipblasDoubleComplex)));
			// cout<<"m:"<<m<<endl;
			// cout<<"grindex:"<<endl;
			// int *GR_index_h = new int[10];
			// CHECK_CUDA(hipMemcpy(GR_index_h, GR_index_d, 10*sizeof(int), hipMemcpyDeviceToHost));
			// for(int i=0;i<10;i++)
			// {
			// 	cout<<GR_index_h[i]<<endl;
			// }
			// GlobalC::UFFT.RoundTrip(tmpsi_in, d_vr_eff1, GR_index_d, d_porter);
			GlobalC::wfcpw->recip2real(tmpsi_in, f_porter, ik);
			for (int ir=0; ir< nrxx; ir++)
			{
				f_porter[ir] *=  GlobalC::pot.vr_eff1[ir];
			}
			GlobalC::wfcpw->real2recip(f_porter, tmhpsi, ik, true);

			// cout<<"dporter"<<endl;
			// print_test<hipblasDoubleComplex>(d_porter, 10);

			// for (j = 0;j < wf.npw;j++)
			// {
			//     tmhpsi[j] += UFFT.porter[ GR_index[j] ];
			// }
			int thread = 512;
			int block = (GlobalC::wf.npw + thread - 1) / thread;
			hipLaunchKernelGGL(HIP_KERNEL_NAME(kernel_add_tmhpsi<double2>),
							   dim3(block),
							   dim3(thread),
							   0,
							   0,
							   GlobalC::wf.npw,
							   reinterpret_cast<double2 *>(tmhpsi),
							   reinterpret_cast<double2 *>(d_porter),
							   GR_index_d);

			tmhpsi += dmax;
			tmpsi_in += dmax;
		}
		// CHECK_CUDA(hipFree(d_GR_index));
		CHECK_CUDA(hipFree(d_vr_eff1));
		CHECK_CUDA(hipFree(d_porter));
	}
	ModuleBase::timer::tick("Hamilt_PW_HIP", "vloc");

	// cout<<"======after hpsi part II======="<<endl;
	// print_test<hipblasDoubleComplex>(hpsi, 10);
	//------------------------------------
	// (3) the nonlocal pseudopotential.
	//------------------------------------
	ModuleBase::timer::tick("Hamilt_PW_HIP", "vnl");

	if (GlobalV::VNL_IN_H)
	{
		if (GlobalC::ppcell.nkb > 0)
		{
			ModuleBase::timer::tick("PART VNL", "malloc vectors");
			int nkb = GlobalC::ppcell.nkb;
			hipblasDoubleComplex *becp;
			// hipblasDoubleComplex *d_vkb_c;
			CHECK_CUDA(hipMalloc((void **)&becp, GlobalV::NPOL * m * nkb * sizeof(hipblasDoubleComplex)));
			// CHECK_CUDA(hipMalloc((void**)&d_vkb_c, GlobalC::wf.npwx*nkb*sizeof(hipblasDoubleComplex)));

			// CHECK_CUDA(hipMemcpy(d_vkb_c, GlobalC::ppcell.vkb.c, GlobalC::wf.npwx*nkb*sizeof(hipblasDoubleComplex),
			// hipMemcpyHostToDevice));

			ModuleBase::timer::tick("PART VNL", "malloc vectors");
			hipblasOperation_t transa = HIPBLAS_OP_C;
			hipblasOperation_t transb = HIPBLAS_OP_N;
			// hipblasHandle_t handle;
			// CHECK_CUBLAS(hipblasCreate(&handle));

			// cout<<"===== vkbc ===="<<endl;
			// print_test<hipblasDoubleComplex>(d_vkb_c, 15);

			// double2 ONE, ZERO;
			// ONE.y = ZERO.x = ZERO.y = 0.0;
			// ONE.x = 1.0;

			hipblasDoubleComplex ONE(1, 0);
			hipblasDoubleComplex ZERO(0, 0);
			hipblasDoubleComplex NEG_ONE(-1, 0);
			// NEG_ONE.x = -1.0;

			ModuleBase::timer::tick("PART VNL", "hipblas gemv");
			if (m == 1 && GlobalV::NPOL == 1)
			{
				int inc = 1;
				CHECK_CUBLAS(hipblasZgemv(hpw_handle,
										  transa,
										  GlobalC::wf.npw,
										  nkb,
										  &ONE,
										  vkb_c,
										  GlobalC::wf.npwx,
										  psi_in,
										  inc,
										  &ZERO,
										  becp,
										  inc));
			}
			else
			{
				int npm = GlobalV::NPOL * m;
				CHECK_CUBLAS(hipblasZgemm(hpw_handle,
										  transa,
										  transb,
										  nkb,
										  npm,
										  GlobalC::wf.npw,
										  &ONE,
										  vkb_c,
										  GlobalC::wf.npwx,
										  psi_in,
										  GlobalC::wf.npwx,
										  &ZERO,
										  becp,
										  nkb));
			}
			ModuleBase::timer::tick("PART VNL", "hipblas gemv");

			// complex<double> *hpsi_cpu = new complex<double>[GlobalC::wf.npw*GlobalV::NPOL];
			// complex<double> *becp_cpu = new complex<double>[GlobalV::NPOL*m*nkb];

			// CHECK_CUDA(hipMemcpy(becp_cpu, becp, GlobalV::NPOL*m*nkb*sizeof(hipblasDoubleComplex),
			// hipMemcpyDeviceToHost));

			// CHECK_CUDA(hipMemcpy(hpsi_cpu, hpsi, GlobalC::wf.npw*GlobalV::NPOL*sizeof(hipblasDoubleComplex),
			// hipMemcpyDeviceToHost));

			// this->add_nonlocal_pp(hpsi_cpu, becp_cpu, m);

			// CHECK_CUDA(hipMemcpy(hpsi, hpsi_cpu, GlobalC::wf.npw*GlobalV::NPOL*sizeof(hipblasDoubleComplex),
			// hipMemcpyHostToDevice));

			// delete [] hpsi_cpu;
			// delete [] becp_cpu;

			// cout<<"===== becp before add nonloaclpp ===="<<endl;
			// print_test<hipblasDoubleComplex>(becp, 15);

			this->add_nonlocal_pp_cuda(hpsi, becp, vkb_c, m);

			// CHECK_CUBLAS(hipblasDestroy(handle));
			CHECK_CUDA(hipFree(becp));
			// CHECK_CUDA(hipFree(d_vkb_c));
			// cout<<"nonlocal end"<<endl;
		}
	}

	ModuleBase::timer::tick("Hamilt_PW_HIP", "vnl");
	// cout<<"======after hpsi part III======="<<endl;
	// print_test<hipblasDoubleComplex>(hpsi, 10);

	//------------------------------------
	// (4) the metaGGA part
	//------------------------------------
	// TODO: add metaGGA part

	ModuleBase::timer::tick("Hamilt_PW_HIP", "h_psi");
	return;
}

void Hamilt_PW::h_psi(const std::complex<double> *psi_in, std::complex<double> *hpsi, const int m)
{
	ModuleBase::timer::tick("Hamilt_PW", "h_psi_cpu");
	// int i = 0;
	int j = 0;
	int ig = 0;

	// if(GlobalV::NSPIN!=4) ZEROS(hpsi, GlobalC::wf.npw);
	// else ZEROS(hpsi, GlobalC::wf.npwx * GlobalV::NPOL);//added by zhengdy-soc
	const int ik = this->current_ik;
	const double tpiba2 = GlobalC::ucell.tpiba2;
	const int npwx = GlobalC::wf.npwx;
	const int npw = npw;
	const int nrxx = GlobalC::rhopw->nrxx;
	int dmax = GlobalC::wf.npwx * GlobalV::NPOL;

	//------------------------------------
	//(1) the kinetical energy.
	//------------------------------------
	std::complex<double> *tmhpsi;
	const std::complex<double> *tmpsi_in;
	if (GlobalV::T_IN_H)
	{
		tmhpsi = hpsi;
		tmpsi_in = psi_in;
		for (int ib = 0; ib < m; ++ib)
		{
			for (ig = 0; ig < GlobalC::wf.npw; ++ig)
			{
				tmhpsi[ig] = GlobalC::wfcpw->getgk2(ik,ig) * GlobalC::ucell.tpiba2 * tmpsi_in[ig];
			}
			if (GlobalV::NSPIN == 4)
			{
				for (ig = GlobalC::wf.npw; ig < GlobalC::wf.npwx; ++ig)
				{
					tmhpsi[ig] = 0;
				}
				tmhpsi += GlobalC::wf.npwx;
				tmpsi_in += GlobalC::wf.npwx;
				for (ig = 0; ig < GlobalC::wf.npw; ++ig)
				{
					tmhpsi[ig] = GlobalC::wfcpw->getgk2(ik,ig) * GlobalC::ucell.tpiba2 * tmpsi_in[ig];
				}
				for (ig = GlobalC::wf.npw; ig < GlobalC::wf.npwx; ++ig)
				{
					tmhpsi[ig] = 0;
				}
			}
			tmhpsi += GlobalC::wf.npwx;
			tmpsi_in += GlobalC::wf.npwx;
		}
	}

	//------------------------------------
	//(2) the local potential.
	//-----------------------------------
	ModuleBase::timer::tick("Hamilt_PW", "vloc");
	if (GlobalV::VL_IN_H)
	{
		tmhpsi = hpsi;
		tmpsi_in = psi_in;
		for(int ib = 0 ; ib < m; ++ib)
		{
			if(GlobalV::NSPIN!=4){
				GlobalC::wfcpw->recip2real(tmpsi_in, porter, ik);
				for (int ir=0; ir< nrxx; ir++)
				{
					porter[ir] *=  GlobalC::pot.vr_eff1[ir];
				}
				GlobalC::wfcpw->real2recip(porter, tmhpsi, ik, true);
			}
			else
			{
				std::complex<double>* porter1 = new std::complex<double>[nrxx];
				// fft to real space and doing things.
				GlobalC::wfcpw->recip2real(tmpsi_in, porter, ik);
				GlobalC::wfcpw->recip2real(tmpsi_in+npwx, porter1, ik);
				std::complex<double> sup,sdown;
				for (int ir=0; ir< nrxx; ir++)
				{
					sup = porter[ir] * (GlobalC::pot.vr_eff(0,ir) + GlobalC::pot.vr_eff(3,ir)) +
						porter1[ir] * (GlobalC::pot.vr_eff(1,ir) - std::complex<double>(0.0,1.0) * GlobalC::pot.vr_eff(2,ir));
					sdown = porter1[ir] * (GlobalC::pot.vr_eff(0,ir) - GlobalC::pot.vr_eff(3,ir)) +
					porter[ir] * (GlobalC::pot.vr_eff(1,ir) + std::complex<double>(0.0,1.0) * GlobalC::pot.vr_eff(2,ir));
					porter[ir] = sup;
					porter1[ir] = sdown;
				}
				// (3) fft back to G space.
				GlobalC::wfcpw->real2recip(porter, tmhpsi, ik, true);
				GlobalC::wfcpw->real2recip(porter1, tmhpsi+npwx, ik, true);

				delete[] porter1;
			}
			tmhpsi += dmax;
			tmpsi_in += dmax;
		}
	}
	ModuleBase::timer::tick("Hamilt_PW", "vloc");

	//------------------------------------
	// (3) the nonlocal pseudopotential.
	//------------------------------------
	ModuleBase::timer::tick("Hamilt_PW", "vnl");
	if (GlobalV::VNL_IN_H)
	{
		if (GlobalC::ppcell.nkb > 0)
		{
			//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
			// qianrui optimize 2021-3-31
			int nkb = GlobalC::ppcell.nkb;
			ModuleBase::ComplexMatrix becp(GlobalV::NPOL * m, nkb, false);
			char transa = 'C';
			char transb = 'N';
			if (m == 1 && GlobalV::NPOL == 1)
			{
				int inc = 1;
				zgemv_(&transa,
					   &GlobalC::wf.npw,
					   &nkb,
					   &ModuleBase::ONE,
					   GlobalC::ppcell.vkb.c,
					   &GlobalC::wf.npwx,
					   psi_in,
					   &inc,
					   &ModuleBase::ZERO,
					   becp.c,
					   &inc);
			}
			else
			{
				int npm = GlobalV::NPOL * m;
				zgemm_(&transa,
					   &transb,
					   &nkb,
					   &npm,
					   &GlobalC::wf.npw,
					   &ModuleBase::ONE,
					   GlobalC::ppcell.vkb.c,
					   &GlobalC::wf.npwx,
					   psi_in,
					   &GlobalC::wf.npwx,
					   &ModuleBase::ZERO,
					   becp.c,
					   &nkb);
				// add_nonlocal_pp is moddified, thus tranpose not needed here.
				// if(GlobalV::NONCOLIN)
				//{
				//	ComplexMatrix partbecp(GlobalV::NPOL, nkb ,false);
				//	for(int ib = 0; ib < m; ++ib)
				//	{
				//
				//		for ( i = 0;i < GlobalV::NPOL;i++)
				//			for (j = 0;j < nkb;j++)
				//				partbecp(i, j) = tmbecp[i*nkb+j];
				//		for (j = 0; j < nkb; j++)
				//			for (i = 0;i < GlobalV::NPOL;i++)
				//				tmbecp[j*GlobalV::NPOL+i] = partbecp(i, j);
				//		tmbecp += GlobalV::NPOL * nkb;
				//	}
				//}
			}

			Parallel_Reduce::reduce_complex_double_pool(becp.c, nkb * GlobalV::NPOL * m);

			this->add_nonlocal_pp(hpsi, becp.c, m);
			//======================================================================
			/*std::complex<double> *becp = new std::complex<double>[ GlobalC::ppcell.nkb * GlobalV::NPOL ];
			ZEROS(becp,GlobalC::ppcell.nkb * GlobalV::NPOL);
			for (i=0;i< GlobalC::ppcell.nkb;i++)
			{
				const std::complex<double>* p = &GlobalC::ppcell.vkb(i,0);
				const std::complex<double>* const p_end = p + GlobalC::wf.npw;
				const std::complex<double>* psip = psi_in;
				for (;p<p_end;++p,++psip)
				{
					if(!GlobalV::NONCOLIN) becp[i] += psip[0]* conj( p[0] );
					else{
						becp[i*2] += psip[0]* conj( p[0] );
						becp[i*2+1] += psip[GlobalC::wf.npwx]* conj( p[0] );
					}
				}
			}
			Parallel_Reduce::reduce_complex_double_pool( becp, GlobalC::ppcell.nkb * GlobalV::NPOL);
			this->add_nonlocal_pp(hpsi, becp);
			delete[] becp;*/
			//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
		}
	}
	ModuleBase::timer::tick("Hamilt_PW", "vnl");
	//------------------------------------
	// (4) the metaGGA part
	//------------------------------------
	// timer::tick("Hamilt_PW","meta");
	if (XC_Functional::get_func_type() == 3)
	{
		tmhpsi = hpsi;
		tmpsi_in = psi_in;
		for(int ib = 0; ib < m; ++ib)
		{
			for(int j=0; j<3; j++)
			{
				for (int ig = 0;ig < GlobalC::kv.ngk[GlobalV::CURRENT_K] ; ig++)
				{
					double fact = GlobalC::wfcpw->getgpluskcar(ik,ig)[j] * tpiba2;
					porter[ig] = tmpsi_in[ig] * complex<double>(0.0,fact);
				}
				
				GlobalC::wfcpw->recip2real(porter,porter,ik);

				for (int ir = 0; ir < nrxx; ir++)
				{
					porter[ir] *= GlobalC::pot.vofk(GlobalV::CURRENT_SPIN,ir);
				}
				GlobalC::wfcpw->real2recip(porter,porter,ik);

				for (int ig = 0;ig < npw ; ig++)
				{
					double fact = GlobalC::wfcpw->getgpluskcar(ik,ig)[j] * tpiba2;
					tmhpsi[ig] -= complex<double>(0.0,fact) * porter[ig];
				}
			}//x,y,z directions
		}
	}
	// timer::tick("Hamilt_PW","meta");
	ModuleBase::timer::tick("Hamilt_PW", "h_psi_cpu");
	return;
}

void Hamilt_PW::add_nonlocal_pp_cuda(hipblasComplex *hpsi_in,
									 const hipblasComplex *becp,
									 const hipblasComplex *f_vkb_c,
									 const int m)
{
	ModuleBase::timer::tick("Hamilt_PW_HIP", "add_nonlocal_pp");

	// number of projectors
	int nkb = GlobalC::ppcell.nkb;

	// complex<double> *ps  = new complex<double> [nkb * GlobalV::NPOL * m];
	// ZEROS(ps, GlobalV::NPOL * m * nkb);
	hipblasComplex *ps;
	CHECK_CUDA(hipMalloc((void **)&ps, nkb * GlobalV::NPOL * m * sizeof(hipblasComplex)));
	CHECK_CUDA(hipMemset(ps, 0, GlobalV::NPOL * m * sizeof(hipblasComplex)));

	int sum = 0;
	int iat = 0;
	// if(GlobalV::NSPIN!=4)
	// {
	for (int it = 0; it < GlobalC::ucell.ntype; it++)
	{

		const int nproj = GlobalC::ucell.atoms[it].nh;
		const int nprojx = GlobalC::ppcell.nhm;
		double *cur_deeq;
		float *f_cur_deeq;
		CHECK_CUDA(hipMalloc((void **)&cur_deeq, nprojx * nprojx * sizeof(double)));
		CHECK_CUDA(hipMalloc((void **)&f_cur_deeq, nprojx * nprojx * sizeof(float)));

		for (int ia = 0; ia < GlobalC::ucell.atoms[it].na; ia++)
		{
			CHECK_CUDA(hipMemcpy(cur_deeq,
								 &(GlobalC::ppcell.deeq(GlobalV::CURRENT_SPIN, iat, 0, 0)),
								 nprojx * nprojx * sizeof(double),
								 hipMemcpyHostToDevice));

			int thread0 = 512;
			int block0 = (nprojx * nprojx + thread0 - 1) / thread0;
			hipLaunchKernelGGL(cast_d2f, dim3(block0), dim3(thread0), 0, 0, f_cur_deeq, cur_deeq, nprojx * nprojx);

			int thread_x = 16;
			dim3 thread(thread_x, thread_x);
			dim3 block((nproj + thread_x - 1) / thread_x, (m + thread_x - 1) / thread_x);
			// dim3 block(1, 1, 1);

			// cout<<"===== ps before add pp kernel ===="<<endl;
			// print_test<hipblasComplex>(ps, 15);

			// cout<<"===== becp before add pp kernel ===="<<endl;
			// print_test<hipblasComplex>((hipblasComplex*)becp, 15);

			// cout<<"===== f_deeq before add pp kernel ===="<<endl;
			// float* test_deeq = (float*)malloc(15*sizeof(float));
			// CHECK_CUDA(hipMemcpy(test_deeq, f_cur_deeq, 15*sizeof(float), hipMemcpyDeviceToHost));
			// for(int i=0;i<15;i++){
			// 	cout<<test_deeq[i]<<endl;
			// }
			// cout<<endl;
			// delete [] test_deeq;

			hipLaunchKernelGGL(HIP_KERNEL_NAME(kernel_addpp<float, float2>),
							   dim3(block),
							   dim3(thread),
							   0,
							   0,
							   reinterpret_cast<float2 *>(ps),
							   f_cur_deeq,
							   reinterpret_cast<const float2 *>(becp),
							   nproj,
							   nprojx,
							   sum,
							   m,
							   nkb);

			sum += nproj;
			++iat;
		} // end na
	} // end nt
	// }

	hipblasOperation_t transa = HIPBLAS_OP_N;
	hipblasOperation_t transb = HIPBLAS_OP_T;
	// hipblasHandle_t handle;
	// CHECK_CUBLAS(hipblasCreate(&handle);
	// float2 ONE;
	// ONE.y = 0.0;
	// ONE.x = 1.0;
	hipblasComplex ONE(1, 0);
	if (GlobalV::NPOL == 1 && m == 1)
	{
		int inc = 1;
		CHECK_CUBLAS(hipblasCgemv(hpw_handle,
								  transa,
								  GlobalC::wf.npw,
								  GlobalC::ppcell.nkb,
								  &ONE,
								  f_vkb_c,
								  GlobalC::wf.npwx,
								  ps,
								  inc,
								  &ONE,
								  hpsi_in,
								  inc));
	}
	else
	{
		int npm = GlobalV::NPOL * m;
		CHECK_CUBLAS(hipblasCgemm(hpw_handle,
								  transa,
								  transb,
								  GlobalC::wf.npw,
								  npm,
								  GlobalC::ppcell.nkb,
								  &ONE,
								  f_vkb_c,
								  GlobalC::wf.npwx,
								  ps,
								  npm,
								  &ONE,
								  hpsi_in,
								  GlobalC::wf.npwx));
	}

	// delete[] ps;
	CHECK_CUDA(hipFree(ps));
	// CHECK_CUBLAS(hipblasDestroy(handle));
	ModuleBase::timer::tick("Hamilt_PW_HIP", "add_nonlocal_pp");
	return;
}

void Hamilt_PW::add_nonlocal_pp_cuda(hipblasDoubleComplex *hpsi_in,
									 const hipblasDoubleComplex *becp,
									 const hipblasDoubleComplex *d_vkb_c,
									 const int m)
{
	ModuleBase::timer::tick("Hamilt_PW_HIP", "add_nonlocal_pp");

	// number of projectors
	int nkb = GlobalC::ppcell.nkb;

	// complex<double> *ps  = new complex<double> [nkb * GlobalV::NPOL * m];
	// ZEROS(ps, GlobalV::NPOL * m * nkb);
	hipblasDoubleComplex *ps;
	CHECK_CUDA(hipMalloc((void **)&ps, nkb * GlobalV::NPOL * m * sizeof(hipblasDoubleComplex)));
	CHECK_CUDA(hipMemset(ps, 0, GlobalV::NPOL * m * sizeof(hipblasDoubleComplex)));

	int sum = 0;
	int iat = 0;
	// if(GlobalV::NSPIN!=4)
	// {
	for (int it = 0; it < GlobalC::ucell.ntype; it++)
	{

		const int nproj = GlobalC::ucell.atoms[it].nh;
		const int nprojx = GlobalC::ppcell.nhm;
		double *cur_deeq;
		CHECK_CUDA(hipMalloc((void **)&cur_deeq, nprojx * nprojx * sizeof(double)));
		for (int ia = 0; ia < GlobalC::ucell.atoms[it].na; ia++)
		{
			CHECK_CUDA(hipMemcpy(cur_deeq,
								 &(GlobalC::ppcell.deeq(GlobalV::CURRENT_SPIN, iat, 0, 0)),
								 nprojx * nprojx * sizeof(double),
								 hipMemcpyHostToDevice));

			int thread_x = 16;
			dim3 thread(thread_x, thread_x);
			dim3 block((nproj + thread_x - 1) / thread_x, (m + thread_x - 1) / thread_x);
			// dim3 block(1, 1, 1);

			// cout<<"===== ps before add pp kernel ===="<<endl;
			// print_test<hipblasDoubleComplex>(ps, 15);

			// cout<<"===== becp before add pp kernel ===="<<endl;
			// print_test<hipblasDoubleComplex>((hipblasDoubleComplex*)becp, 15);

			// cout<<"===== f_deeq before add pp kernel ===="<<endl;
			// double* test_deeq = (double*)malloc(15*sizeof(double));
			// CHECK_CUDA(hipMemcpy(test_deeq, cur_deeq, 15*sizeof(double), hipMemcpyDeviceToHost));
			// for(int i=0;i<15;i++){
			// 	cout<<test_deeq[i]<<endl;
			// }
			// cout<<endl;
			// delete [] test_deeq;

			hipLaunchKernelGGL(HIP_KERNEL_NAME(kernel_addpp<double, double2>),
							   dim3(block),
							   dim3(thread),
							   0,
							   0,
							   reinterpret_cast<double2 *>(ps),
							   cur_deeq,
							   reinterpret_cast<const double2 *>(becp),
							   nproj,
							   nprojx,
							   sum,
							   m,
							   nkb);

			sum += nproj;
			++iat;
		} // end na
	} // end nt
	// }

	hipblasOperation_t transa = HIPBLAS_OP_N;
	hipblasOperation_t transb = HIPBLAS_OP_T;
	// hipblasHandle_t handle;
	// CHECK_CUBLAS(hipblasCreate(&handle));
	// double2 ONE;
	// ONE.y = 0.0;
	// ONE.x = 1.0;
	hipblasDoubleComplex ONE(1, 0);

	if (GlobalV::NPOL == 1 && m == 1)
	{
		int inc = 1;
		CHECK_CUBLAS(hipblasZgemv(hpw_handle,
								  transa,
								  GlobalC::wf.npw,
								  GlobalC::ppcell.nkb,
								  &ONE,
								  d_vkb_c,
								  GlobalC::wf.npwx,
								  ps,
								  inc,
								  &ONE,
								  hpsi_in,
								  inc));
	}
	else
	{
		int npm = GlobalV::NPOL * m;
		CHECK_CUBLAS(hipblasZgemm(hpw_handle,
								  transa,
								  transb,
								  GlobalC::wf.npw,
								  npm,
								  GlobalC::ppcell.nkb,
								  &ONE,
								  d_vkb_c,
								  GlobalC::wf.npwx,
								  ps,
								  npm,
								  &ONE,
								  hpsi_in,
								  GlobalC::wf.npwx));
	}

	// delete[] ps;
	CHECK_CUDA(hipFree(ps));
	// CHECK_CUBLAS(hipblasDestroy(handle));
	ModuleBase::timer::tick("Hamilt_PW_HIP", "add_nonlocal_pp");
	return;
}

void Hamilt_PW::add_nonlocal_pp(std::complex<double> *hpsi_in, const std::complex<double> *becp, const int m)
{
	ModuleBase::timer::tick("Hamilt_PW", "add_nonlocal_pp");

	// number of projectors
	int nkb = GlobalC::ppcell.nkb;

	std::complex<double> *ps = new std::complex<double>[nkb * GlobalV::NPOL * m];
	ModuleBase::GlobalFunc::ZEROS(ps, GlobalV::NPOL * m * nkb);

	int sum = 0;
	int iat = 0;
	if (GlobalV::NSPIN != 4)
	{
		for (int it = 0; it < GlobalC::ucell.ntype; it++)
		{
			const int nproj = GlobalC::ucell.atoms[it].nh;
			for (int ia = 0; ia < GlobalC::ucell.atoms[it].na; ia++)
			{
				// each atom has nproj, means this is with structure factor;
				// each projector (each atom) must multiply coefficient
				// with all the other projectors.
				for (int ip = 0; ip < nproj; ip++)
				{
					for (int ip2 = 0; ip2 < nproj; ip2++)
					{
						for (int ib = 0; ib < m; ++ib)
						{
							ps[(sum + ip2) * m + ib] += GlobalC::ppcell.deeq(GlobalV::CURRENT_SPIN, iat, ip, ip2)
														* becp[ib * nkb + sum + ip];
						} // end ib
					} // end ih
				} // end jh
				sum += nproj;
				++iat;
			} // end na
		} // end nt
	}
	else
	{
		for (int it = 0; it < GlobalC::ucell.ntype; it++)
		{
			int psind = 0;
			int becpind = 0;
			std::complex<double> becp1 = std::complex<double>(0.0, 0.0);
			std::complex<double> becp2 = std::complex<double>(0.0, 0.0);

			const int nproj = GlobalC::ucell.atoms[it].nh;
			for (int ia = 0; ia < GlobalC::ucell.atoms[it].na; ia++)
			{
				// each atom has nproj, means this is with structure factor;
				// each projector (each atom) must multiply coefficient
				// with all the other projectors.
				for (int ip = 0; ip < nproj; ip++)
				{
					for (int ip2 = 0; ip2 < nproj; ip2++)
					{
						for (int ib = 0; ib < m; ++ib)
						{
							psind = (sum + ip2) * 2 * m + ib * 2;
							becpind = ib * nkb * 2 + sum + ip;
							becp1 = becp[becpind];
							becp2 = becp[becpind + nkb];
							ps[psind] += GlobalC::ppcell.deeq_nc(0, iat, ip2, ip) * becp1
										 + GlobalC::ppcell.deeq_nc(1, iat, ip2, ip) * becp2;
							ps[psind + 1] += GlobalC::ppcell.deeq_nc(2, iat, ip2, ip) * becp1
											 + GlobalC::ppcell.deeq_nc(3, iat, ip2, ip) * becp2;
						} // end ib
					} // end ih
				} // end jh
				sum += nproj;
				++iat;
			} // end na
		} // end nt
	}

	/*
	for (int ig=0;ig<GlobalC::wf.npw;ig++)
	{
		for (int i=0;i< GlobalC::ppcell.nkb;i++)
		{
			hpsi_in[ig]+=ps[i]*GlobalC::ppcell.vkb(i,ig);
		}
	}
	*/

	// use simple method.
	//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
	// qianrui optimize 2021-3-31
	char transa = 'N';
	char transb = 'T';
	if (GlobalV::NPOL == 1 && m == 1)
	{
		int inc = 1;
		zgemv_(&transa,
			   &GlobalC::wf.npw,
			   &GlobalC::ppcell.nkb,
			   &ModuleBase::ONE,
			   GlobalC::ppcell.vkb.c,
			   &GlobalC::wf.npwx,
			   ps,
			   &inc,
			   &ModuleBase::ONE,
			   hpsi_in,
			   &inc);
	}
	else
	{
		int npm = GlobalV::NPOL * m;
		zgemm_(&transa,
			   &transb,
			   &GlobalC::wf.npw,
			   &npm,
			   &GlobalC::ppcell.nkb,
			   &ModuleBase::ONE,
			   GlobalC::ppcell.vkb.c,
			   &GlobalC::wf.npwx,
			   ps,
			   &npm,
			   &ModuleBase::ONE,
			   hpsi_in,
			   &GlobalC::wf.npwx);
	}

	//======================================================================
	/*if(!GlobalV::NONCOLIN)
	for(int i=0; i<GlobalC::ppcell.nkb; i++)
	{
		std::complex<double>* p = &GlobalC::ppcell.vkb(i,0);
		std::complex<double>* p_end = p + GlobalC::wf.npw;
		std::complex<double>* hp = hpsi_in;
		std::complex<double>* psp = &ps[i];
		for (;p<p_end;++p,++hp)
		{
			hp[0] += psp[0] * p[0];
		}
	}
	else
	for(int i=0; i<GlobalC::ppcell.nkb; i++)
	{
		std::complex<double>* p = &GlobalC::ppcell.vkb(i,0);
		std::complex<double>* p_end = p + GlobalC::wf.npw;
		std::complex<double>* hp = hpsi_in;
		std::complex<double>* hp1 = hpsi_in + GlobalC::wf.npwx;
		std::complex<double>* psp = &ps[i*2];
		for (;p<p_end;p++,++hp,++hp1)
		{
			hp[0] += psp[0] * (p[0]);
			hp1[0] += psp[1] * (p[0]);
		}
	}*/
	//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

	delete[] ps;
	ModuleBase::timer::tick("Hamilt_PW", "add_nonlocal_pp");
	return;
}
