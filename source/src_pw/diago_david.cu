#include "diago_david.cuh"
#include "diago_cg.cuh"
#include "global.h"
using namespace CudaCheck;

template<class T, class T2>
__global__ void kernel_precondition_david(T2 *res, const T2 *data, const int size, const T *P)
{
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if(idx < size)
    {
        res[idx].x = data[idx].x / P[idx];
        res[idx].y = data[idx].y / P[idx];
    }
}

template<class T, class T2>
__global__ void kernel_add_respsi(T2* respsi, int nbase, int npw, T2* vc, T2* hp, T2* sp, int unconv_m, T e_unconv_m);
{
    int ig = blockDim.x * blockIdx.x + threadIdx.x;
    for ( int i = 0; i < nbase; i++ ) 
    {
        respsi[ig] += vc[i*npw + unconv_m] * (hp[i*npw + ig] - e_unconv_m * sp[i*npw + ig]);
    }
}

template<class T, class T2>
__global__ void kernel_normalization_david(T2 *data, int size, T norm)
{
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if(idx < size)
    {
        data[idx].x /= norm;
        data[idx].y /= norm;
    }
}

Diago_David_CUDA::Diago_David_CUDA()
{
    test_david = 2;
    // 1: check which function is called and which step is executed
    // 2: check the eigenvalues of the result of each iteration
    // 3: check the eigenvalues and errors of the last result
    // default: no check
}

Diago_David_CUDA::~Diago_David_CUDA() 
{
}


void Diago_David_CUDA::diag
(
    double2* psi, // matrix
    double* en,
    const int& npw,
    const int& nband,
    const double* precondition,
    const int order,
    const double& eps,
    const int& maxiter,
    int& notconv,
    double& avg_iter
)
{
    if (test_david==1) ModuleBase::TITLE("Diago_David_CUDA","diag");
    ModuleBase::timer::tick("Diago_David_CUDA", "diag");

    assert( order > 1 );
    assert( order*nband < npw * GlobalV::NPROC_IN_POOL ); 
    //In strictly speaking, it shoule be order*nband < npw sum of all pools. We roughly estimate it here.
    //However, in most cases, total number of plane waves should be much larger than nband*order

    int nbase_x = order * nband ;				// maximum dimension of the reduced basis set

    // ModuleBase::ComplexMatrix basis( nbase_x, npw );		// the reduced basis set
    // ModuleBase::ComplexMatrix hp( nbase_x, npw );			// the product of H and psi in the reduced basis set
    // ModuleBase::ComplexMatrix sp( nbase_x, npw );			// the Product of S and psi in the reduced basis set

    double2* basis;
    double2* hp;
    double2* sp;
    CHECK_CUDA(cudaMalloc((void**)&basis, nbase_x*npw*sizeof(double2)));
    CHECK_CUDA(cudaMalloc((void**)&hp, nbase_x*npw*sizeof(double2)));
    CHECK_CUDA(cudaMalloc((void**)&sp, nbase_x*npw*sizeof(double2)));

    double2* hc;
    double2* sc;
    double2* vc;
    double* e;
    CHECK_CUDA(cudaMalloc((void**)&hc, nbase_x*nbase_x*sizeof(double2)));
    CHECK_CUDA(cudaMalloc((void**)&sc, nbase_x*nbase_x*sizeof(double2)));
    CHECK_CUDA(cudaMalloc((void**)&vc, nbase_x*nbase_x*sizeof(double2)));
    // CHECK_CUDA(cudaMalloc((void**)&e, nbase_x*sizeof(double)));
    assert(e != 0);

    double2* psi_m;
    double2* hpsi;
    double2* spsi;
    double2* ppsi;
    double2* respsi;
    CHECK_CUDA(cudaMalloc((void**)&psi_m, npw*sizeof(double2)));
    CHECK_CUDA(cudaMalloc((void**)&hpsi, npw*sizeof(double2)));
    CHECK_CUDA(cudaMalloc((void**)&spsi, npw*sizeof(double2)));
    CHECK_CUDA(cudaMalloc((void**)&ppsi, npw*sizeof(double2)));
    CHECK_CUDA(cudaMalloc((void**)&respsi, npw*sizeof(double2)));

    // bool* convflag = new bool[nband] ;	// convflag[m] = true if the m th band is convergent
    // assert(convflag != 0) ;
    // int* unconv = new int[nband] ;		// unconv[m] store the number of the m th unconvergent band
    // assert(unconv != 0) ;

    bool* convflag = new bool[nband];
    int* unconv = new int[nband];
    int nbase = 0;						// the dimension of the reduced basis set
    notconv = nband;					// the number of the unconvergent bands
    // ModuleBase::GlobalFunc::ZEROS( convflag, nband );
    for ( int m = 0 ; m < nband; m++ ) unconv[m] = m;

    ModuleBase::timer::tick("Diago_David_CUDA","first");
    // orthogonalise the initial trial psi(0~nband-1)
    for (int m = 0; m < nband; m++)
    {
        // psi_m = psi(m)
        CHECK_CUDA(cudaMemcpy(psi_m, psi[npw*m], npw*sizeof(double2), cudaMemcpyDeviceToDevice));

        this->SchmitOrth(npw, nband, m, basis, psi_m, spsi);

        GlobalC::hm.hpw.h_1psi_cuda(npw, psi_m, hpsi, spsi);

        // basis(m) = psi_m, hp(m) = H |psi_m>, sp(m) = S |psi_m>
        CHECK_CUDA(cudaMemcpy(basis[npw*m], psi_m, npw*sizeof(double2), cudaMemcpyDeviceToDevice));
        CHECK_CUDA(cudaMemcpy(hp[npw*m], hpsi, npw*sizeof(double2),cudaMemcpyDeviceToDevice));
        CHECK_CUDA(cudaMemcpy(sp[npw*m], spsi, npw*sizeof(double2), cudaMemcpyDeviceToDevice));
    }

    // TODO: zero_out ...
    hc.zero_out();
    sc.zero_out();

    this->cal_elem( npw, nbase, nbase_x, notconv, basis, hp, sp, hc, sc );

    this->diag_zhegvx( nbase, nband, hc, sc, nbase_x, e, vc );

    for ( int m = 0; m < nband; m++ ) 
	{
		en[m] = e[m];
	}
    // (en & e is on host.)
    // CHECK_CUDA(cudaMemcpy(en, e, nband*sizeof(double), cudaMemcpyHostToDevice));

    ModuleBase::timer::tick("Diago_David_CUDA","first");

    int dav_iter = 0;
    do
    {
        dav_iter++;

        this->cal_grad( npw, nbase, notconv, basis, hp, sp, vc, 
			unconv, precondition, e, hpsi, spsi, ppsi, respsi );

        this->cal_elem( npw, nbase, notconv, basis, hp, sp, hc, sc );

        this->diag_zhegvx( nbase, nband, hc, sc, nbase_x, e, vc );

        // check convergence and update eigenvalues
        ModuleBase::timer::tick("Diago_David_CUDA","check_update");

        notconv = 0;
        for ( int m = 0 ; m < nband; m++ )
        {
            convflag[m] = ( abs( e[m] - en[m] ) < eps );

            if ( !convflag[m] ) 
			{
                unconv[notconv] = m;
                notconv++;
            }

            en[m] = e[m];
        }

        ModuleBase::timer::tick("Diago_David_CUDA","check_update");

        if ( !notconv || ( nbase + notconv > nbase_x) || (dav_iter == maxiter) )
        {
            ModuleBase::timer::tick("Diago_David_CUDA","last");

            // updata eigenvectors of Hamiltonian
            psi.zero_out();
            for ( int m = 0; m < nband; m++ )
            {
                for ( int j = 0; j < nbase; j++ )
                {
                    for ( int ig = 0; ig < npw; ig++ ) 
					{
						psi(m,ig) += vc(j,m) * basis(j,ig);
					}
                }
            }
            // add kernel

            if ( !notconv || (dav_iter == maxiter) )
            {
                // overall convergence or last iteration: exit the iteration

                ModuleBase::timer::tick("Diago_David_CUDA","last");
                break;
            }
            else
            {
                // if the dimension of the reduced basis set is becoming too large,
                // then replace the first N (=nband) basis vectors with the current
                // estimate of the eigenvectors and set the basis dimension to N;

                this->refresh( npw, nband, nbase, en, psi, basis, hp, sp, hc, sc, vc );
                ModuleBase::timer::tick("Diago_David_CUDA","last");
            }

        }// end of if

    } while (1);

    avg_iter = static_cast<double>(dav_iter);
    /*
    	// if( test_david==3 ) std::cout info of davidson diag
    	if( test_david==3 )
    	{
    		std::cout << "hamilt davidson diag output " << std::endl ;
    		std::cout << "dav_iter = " << dav_iter <<" notconv = " << notconv << std::endl;
    		this->cal_err( npw, nband, nbase, vc, hp, basis, en, respsi );
        	std::cout << "hamilt davidson diag output  end " << std::endl ;
    	}
    */
    delete[] e;
    delete[] psi_m;
    delete[] hpsi;
    delete[] spsi;
    delete[] ppsi;
    delete[] respsi;
    delete[] convflag;
    delete[] unconv;

    ModuleBase::timer::tick("Diago_David_CUDA", "diag");
    return;
}

void Diago_David_CUDA::cal_grad
(
    const int& npw,
    const int& nbase,	// current dimension of the reduced basis
    const int& notconv,
    double2* basis, // matrix
    double2* hp, // matrix
    double2* sp, // matrix
    const double2* vc,
    const int* unconv,
    const double* precondition,
    const double* e,
    double2* hpsi,
    double2* spsi,
    double2* ppsi,
    double2* respsi
)
{
    if ( test_david ==1 ) ModuleBase::TITLE("DIAGO_DAVID","cal_grad");
    ModuleBase::timer::tick("Diago_David_CUDA", "cal_grad"
    );

    // expand the reduced basis set with the new basis vectors P|R(psi)>...
    // in which psi are the last eigenvectors
    // we define |R(psi)> as (H-ES)*|Psi>, E = <psi|H|psi>/<psi|S|psi>
    for ( int m = 0; m < notconv; m++ )
    {
        // ModuleBase::GlobalFunc::ZEROS( respsi, npw );
        // for ( int i = 0; i < nbase; i++ )
        // {
        //     for ( int ig = 0; ig < npw; ig++ )
        //     {
        //         respsi[ig] += vc( i, unconv[m] ) * ( hp(i,ig) - e[ unconv[m] ] * sp(i,ig) ) ;
        //     }
        // }
        int thread = 512;
        int block = (npw + thread - 1) / thread;
        kernel_add_respsi<double, double2><<<block, thread>>>(respsi, nbase, npw, vc, hp, sp, unconv[m], e[unconv[m]]);

        // for ( int ig = 0; ig < npw; ig++ ) 
		// {
		// 	ppsi[ig] = respsi[ig] / precondition[ig] ;
		// } 
        kernel_precondition_david<double, double2><<<block, thread>>>(ppsi, respsi, npw, precondition);

        this->SchmitOrth(npw, nbase+notconv, nbase+m, basis, ppsi, spsi);

        GlobalC::hm.hpw.h_1psi(npw, ppsi, hpsi, spsi);

        // for ( int ig = 0; ig < npw; ig++ )
        // {
        //     basis(nbase+m,ig) = ppsi[ig];
        //     hp(nbase+m,ig) = hpsi[ig];
        //     sp(nbase+m,ig) = spsi[ig];
        // }
        CHECK_CUDA(cudaMemcpy(basis[(nbase+m)*npw], ppsi, npw*sizeof(double2), cudaMemcpyDeviceToDevice));
        CHECK_CUDA(cudaMemcpy(hp[(nbase+m)*npw], hpsi, npw*sizeof(double2), cudaMemcpyDeviceToDevice));
        CHECK_CUDA(cudaMemcpy(sp[(nbase+m)*npw], spsi, npw*sizeof(double2), cudaMemcpyDeviceToDevice));
    }

    ModuleBase::timer::tick("Diago_David_CUDA","cal_grad");
    return;
}

void Diago_David_CUDA::cal_elem
(
    const int& npw,
    int& nbase,			// current dimension of the reduced basis
    int& nbase_x,
    const int& notconv,	// number of newly added basis vectors
    const double2* basis, // matrix
    const double2* hp, // matrix
    const double2* sp, // matrix
    double2* hc, // nbase_x * nbase_x
    double2* sc // nbase_x * nbase_x
)
{
    if ( test_david ==1 ) ModuleBase::TITLE("DIAGO_DAVID_CUDA","cal_elem");
    ModuleBase::timer::tick("Diago_David_CUDA","cal_elem");

    // updat the reduced Hamiltonian
    int offset_h = nbase * nbase_x ;
    int offset_s = nbase * nbase_x ;
//	ModuleBase::GlobalFunc::ZEROS( hc.c+offset_h, notconv*hc.nr );
//	ModuleBase::GlobalFunc::ZEROS( sc.c+offset_s, notconv*sc.nr );

    for ( int i = nbase; i < nbase+notconv; i++ )
    {
        for ( int j = 0; j <= i; j++ )
        {
            hc[j*nbase_x+i] = ddot(npw, basis, j, hp, i);
            sc[j*nbase_x+i] = ddot(npw, basis, j, sp, i);
        }
    }

    // Parallel_Reduce::reduce_complex_double_pool( hc.c+offset_h, notconv*hc.nr );
    // Parallel_Reduce::reduce_complex_double_pool( sc.c+offset_s, notconv*sc.nr );
    /*
    	for( int i = nbase; i < nbase+notconv; i++ )
    	{
    		for( int j = 0; j <i; j++ )
    		{
    			hc(j,i) = conj( hc(i,j) );
    			sc(j,i) = conj( sc(i,j) );
    		}
    	}
    */
    nbase += notconv;
    ModuleBase::timer::tick("Diago_David_CUDA","cal_elem");
    return;
}

void Diago_David_CUDA::diag_zhegvx
(
    const int& n,
    const int& m,
    const double2* hc, // matrix nbase_x * nbase_x
    const double2* sc, // matrix nbase_x * nbase_x
    const int& ldh,
    double* e,
    double2* vc // matrix nbase_x * nbase_x
)
{
//	ModuleBase::TITLE("DIAGO_DAVID","diag_zhegvx");
    ModuleBase::timer::tick("Diago_David_CUDA","diag_zhegvx");
    assert( ldh >= max(1,n) );
    int lwork ;
    int info = 0;
    int mm = m;
    std::string name1 = "ZHETRD";
    std::string name2 = "L";

    int nb = LapackConnector::ilaenv(1, name1.c_str(), name2.c_str(), n, -1, -1, -1);
    if (nb < 1)
    {
        nb = max(1, n);
    }
    if (nb == 1 || nb >= n)
    {
        lwork = 2 * n; //qianrui fix a bug 2021-7-25 : lwork should be at least max(1,2*n)
    }
    else
    {
        lwork = (nb + 1) * n;
    }
    double2 *work = new double2[2*lwork];
    assert(work != 0);
    double *rwork = new double[7*n];
    assert(rwork != 0);
    int *iwork = new int [5*n];
    assert(iwork != 0);
    int *ifail = new int[n];
    assert(ifail != 0);
    ModuleBase::GlobalFunc::ZEROS(work,lwork); //qianrui change it, only first lwork numbers are used in zhegvx
    ModuleBase::GlobalFunc::ZEROS(rwork,7*n);
    ModuleBase::GlobalFunc::ZEROS(iwork,5*n);
    ModuleBase::GlobalFunc::ZEROS(ifail,n);

	//ModuleBase::WARNING_QUIT("divid","open zhegvx!");
	
	LapackConnector::zhegvx(1, 'V', 'I', 'L', n, hc, n, sc, n,
           0.0, 0.0, 1, m, 0.0, mm, e, vc, n,
           work, lwork, rwork, iwork, ifail, info);
/*
		double2 vc_norm = 0.0;
		for(int ib =0; ib < n; ib++)
		{
			vc_norm += conj(vc(ib, 0)) * vc(ib, 0);
		}
		assert( vc_norm.real() > 0.0);
*/
    delete[] work;
    delete[] rwork;
    delete[] iwork;
    delete[] ifail;
    ModuleBase::timer::tick("Diago_David_CUDA","diag_zhegvx");
    return;
}

void Diago_David_CUDA::refresh
(
    const int& npw,
    const int& nband,
    int& nbase,
    const double* en,
    const ModuleBase::ComplexMatrix &psi,
    ModuleBase::ComplexMatrix &basis,
    ModuleBase::ComplexMatrix &hp,
    ModuleBase::ComplexMatrix &sp,
    ModuleBase::ComplexMatrix &hc,
    ModuleBase::ComplexMatrix &sc,
    ModuleBase::ComplexMatrix &vc
)
{
    if ( test_david==1 ) ModuleBase::TITLE("Diago_David_CUDA","refresh");
    ModuleBase::timer::tick("Diago_David_CUDA","refresh");

    // update hp,sp
    basis.zero_out();
    for ( int m = 0; m < nband; m++ )
    {
        for ( int j = 0; j < nbase; j++ )
        {
            for ( int ig = 0; ig < npw; ig++ )
            {
                basis(m, ig) += vc(j,m) * hp(j,ig);
                basis(m+nband, ig) +=vc(j,m) * sp(j,ig);
            }
        }
    }

    for ( int m = 0; m < nband; m++ )
    {
        for (int ig = 0; ig < npw; ig++)
        {
            hp(m,ig) = basis(m, ig);
            sp(m,ig) = basis(m+nband, ig);
        }
    }

    // update basis
    basis.zero_out();
    for ( int m = 0; m < nband; m++ )
    {
        for ( int ig = 0; ig < npw; ig++ ) basis(m,ig) = psi(m,ig);
    }

    // updata the reduced Hamiltonian
    nbase = nband;
    hc.zero_out();
    sc.zero_out();
    for ( int i = 0; i < nbase; i++ )
    {
        hc(i,i) = en[i];
        sc(i,i) = ModuleBase::ONE;
        vc(i,i) = ModuleBase::ONE;
    }

    ModuleBase::timer::tick("Diago_David_CUDA","refresh");
    return;
}

void Diago_David_CUDA::cal_err
(
    const int& npw,
    const int& nband,
    const int& nbase,
    const ModuleBase::ComplexMatrix &vc,
    const ModuleBase::ComplexMatrix &hp,
    const ModuleBase::ComplexMatrix &basis,
    const double* en,
    double2* respsi
)
{
    ModuleBase::timer::tick("Diago_David_CUDA","cal_err");
    double *err = new double[nband];
    assert(err != 0);

    for ( int m = 0; m < nband; m++ )
    {
        ModuleBase::GlobalFunc::ZEROS( respsi, npw );
        for ( int j = 0; j < nbase; j++ )
        {
            for ( int ig=0; ig<npw; ig++ )
            {
                respsi[ig] +=  vc(j,m)*( hp(j,ig) - en[m] * basis(j,ig) );
            }
        }

        err[m] = ModuleBase::GlobalFunc::ddot_real( npw, respsi, respsi );
        err[m] = sqrt( err[m] );
    }

    std::cout << "i       eigenvalues       err (||*||)   " << std::endl ;

    for (int i = 0; i < nband ; i++)
    {
        std::cout << i << "\t\t"  << en[i]  << "\t\t" << err[i] << std::endl ;
    }

    delete[] err;
    ModuleBase::timer::tick("Diago_David_CUDA","cal_err");
    return;
}

void Diago_David_CUDA::SchmitOrth
(
    const int& npw,
    const int n_band,
    const int m,
    const double2 *psi, // matrix
    double2* psi_m,
    double2* spsi
)
{
//	if(test_david == 1) ModuleBase::TITLE("Diago_David_CUDA","SchmitOrth");
    ModuleBase::timer::tick("Diago_David_CUDA","SchmitOrth");

    // orthogonalize starting eigenfunction to those already calculated
    // psi_m orthogonalize to psi(0) ~ psi(m-1)
    // Attention, the orthogonalize here read as
    // psi(m) -> psi(m) - \sum_{i < m} \langle psi(i)|S|psi(m) \rangle psi(i)
    // so the orthogonalize is performed about S.

    // assert(psi.nr >= n_band); // TODO: Add 'nr' parameter?...
    assert(m >= 0);
    assert(m < n_band);

    GlobalC::hm.hpw.s_1psi_cuda(npw, psi_m, spsi);

    // double2* lagrange = new double2[m+1];
    // ModuleBase::GlobalFunc::ZEROS( lagrange, m+1 );

    double2* lagrange;
    CHECK_CUDA(cudaMalloc((void**)&lagrange, (m+1)*sizeof(double2)));
    // TODO: Memset

    // replace with 'gemv'.
    // TODO: Global handle ...
    cublasHandle_t diag_handle;
    CHECK_CUBLAS(cublasCreate(&diag_handle));
    cublasOperation_t trans1 = CUBLAS_OP_C;

    double2 ONE, ZERO, NEG_ONE;
    ONE.y = ZERO.x = ZERO.y = 0.0;
    ONE.x = 1.0;
    NEG_ONE.x = -1.0;
    CHECK_CUBLAS(cublasZgemv(diag_handle, trans1, npw, m+1, &ONE, psi, npw, spsi, 1, &ZERO, lagrange, 1));
    
//	out.printr1_d("lagrange", lagrange, m+1 );

    double psi_norm;
    CHECK_CUDA(cudaMemcpy(&psi_norm, &lagrange[m], sizeof(double), cudaMemcpyDeviceToHost));
    assert(psi_norm > 0.0);
//	std::cout << "m = " << m << std::endl;

    cublasOperation_t trans2 = CUBLAS_OP_N;
    CHECK_CUBLAS(cublasZgemv(diag_handle, trans2, npw, m, &NEG_ONE, psi, npw, lagrange, 1, &ONE, psi_m, 1));
    // TODO: ddot_real
    psi_norm -= ddot_real(m, lagrange, lagrange);
    assert(psi_norm > 0.0);
    psi_norm = sqrt(psi_norm);

    if (psi_norm < 1.0e-12 ) 
	{
        std::cout << "Diago_David_CUDA::SchmitOrth:aborted for psi_norm <1.0e-12" << std::endl;
        std::cout << "n_band = " << n_band << std::endl;
        std::cout << "m = " << m << std::endl;
        exit(0);
    }
    else 
	{
        int thread = 512;
        int block = (npw + thread - 1) / thread;
        kernel_normalization_david<double, double2><<<block, thread>>>(psi_m, dim, psi_norm);
    }

    GlobalC::hm.hpw.s_1psi_cuda(npw, psi_m, spsi);

    CHECK_CUDA(cudaFree(lagrange));
    ModuleBase::timer::tick("Diago_David_CUDA","SchmitOrth");
    return;
}
