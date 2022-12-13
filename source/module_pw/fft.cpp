#include "fft.h"

#include "module_base/global_variable.h"

#ifdef _OPENMP
#include <omp.h>
#endif

namespace ModulePW
{

FFT::FFT()
{
// need to call multi-threads init
// ref: https://www.fftw.org/fftw3_doc/Usage-of-Multi_002dthreaded-FFTW.html
#ifdef _OPENMP
	fftw_init_threads();
	fftw_plan_with_nthreads(omp_get_max_threads());
#endif
}

FFT::~FFT()
{
	this->clear();
// need to call multi-threads cleanup
// ref: https://www.fftw.org/fftw3_doc/Usage-of-Multi_002dthreaded-FFTW.html
#ifdef _OPENMP
	fftw_cleanup_threads();
#endif
}
void FFT::clear()
{
	this->cleanFFT();
	if(auxg!=nullptr) {fftw_free(auxg); auxg = nullptr;}
	if(auxr!=nullptr) {fftw_free(auxr); auxr = nullptr;}
	r_rspace = nullptr;
#if defined(__CUDA) || defined(__UT_USE_CUDA)
    if (GlobalV::device_flag == "gpu") {
        if (auxg_3d != nullptr) {
            cudaFree(auxg_3d);
            auxg_3d = nullptr;
        }
        if (auxr_3d != nullptr) {
            cudaFree(auxr_3d);
            auxr_3d = nullptr;
        }
    }
#endif
#ifdef __MIX_PRECISION
	this->cleanfFFT();
	if(auxfg!=nullptr) {fftw_free(auxfg); auxfg = nullptr;}
	if(auxfr!=nullptr) {fftw_free(auxfr); auxfr = nullptr;}
	rf_rspace = nullptr;
#endif
}

void FFT:: initfft(int nx_in, int ny_in, int nz_in, int lixy_in, int rixy_in, int ns_in, int nplane_in, int nproc_in, bool gamma_only_in, bool xprime_in, bool mpifft_in)
{
	this->gamma_only = gamma_only_in;
	this->xprime = xprime_in;
	this->fftnx = this->nx = nx_in;
	this->fftny = this->ny = ny_in;
	if(this->gamma_only)
	{
		if(xprime) 	this->fftnx = int(nx/2) +1;
		else		this->fftny = int(ny/2) +1;
	} 
	this->nz = nz_in;
	this->ns = ns_in;
	this->lixy = lixy_in;
	this->rixy = rixy_in;
	this->nplane = nplane_in;
	this->nproc = nproc_in;
	this->mpifft = mpifft_in;
	this->nxy = this->nx * this-> ny;
	this->fftnxy = this->fftnx * this->fftny;
	// this->maxgrids = (this->nz * this->ns > this->nxy * nplane) ? this->nz * this->ns : this->nxy * nplane;
	const int nrxx = this->nxy * this->nplane;
	const int nsz = this->nz * this->ns;
	int maxgrids = (nsz > nrxx) ? nsz : nrxx;
	if(!this->mpifft)
	{
		auxg  = (std::complex<double> *) fftw_malloc(sizeof(fftw_complex) * maxgrids);
		auxr  = (std::complex<double> *) fftw_malloc(sizeof(fftw_complex) * maxgrids);
		r_rspace = (double *) auxg;
        // auxg_3d = static_cast<std::complex<double> *>(
        //     fftw_malloc(sizeof(fftw_complex) * (this->nx * this->ny * this->nz)));
        // auxr_3d = static_cast<std::complex<double> *>(
        //     fftw_malloc(sizeof(fftw_complex) * (this->nx * this->ny * this->nz)));
        #if defined(__CUDA) || defined(__UT_USE_CUDA)
        if (GlobalV::device_flag == "gpu") {
            cudaMalloc(reinterpret_cast<void **>(&auxg_3d),
                       this->nx * this->ny * this->nz * sizeof(std::complex<double>));
            cudaMalloc(reinterpret_cast<void **>(&auxr_3d),
                       this->nx * this->ny * this->nz * sizeof(std::complex<double>));
        }
        #endif
#ifdef __MIX_PRECISION
		auxfg  = (std::complex<float> *) fftw_malloc(sizeof(fftwf_complex) * maxgrids);
		auxfr  = (std::complex<float> *) fftw_malloc(sizeof(fftwf_complex) * maxgrids);
		rf_rspace = (float *) auxfg;
#endif
	}
	else
	{
		
	}
	
}

void FFT:: setupFFT()
{
	if(!this->mpifft)
	{
		this->initplan();
#ifdef __MIX_PRECISION
		this->initplanf();
#endif
	}
#if defined(__FFTW3_MPI) && defined(__MPI)
	else
	{
		this->initplan_mpi();
#ifdef __MIX_PRECISION
		this->initplanf_mpi();
#endif
	}
#endif
	return;
}
	
void FFT :: initplan()
{
	//---------------------------------------------------------
	//                              1 D - Z
	//---------------------------------------------------------

	//               fftw_plan_many_dft(int rank,          const int *n,       int howmany,
	//					                fftw_complex *in,  const int *inembed, int istride, int idist, 
	//					                fftw_complex *out, const int *onembed, int ostride, int odist, int sign, unsigned flags);
	
	this->planzfor = fftw_plan_many_dft(     1,    &this->nz,  this->ns,  
					    (fftw_complex*) auxg,  &this->nz,  1,  this->nz,
					    (fftw_complex*) auxg,  &this->nz,  1,  this->nz,  FFTW_FORWARD,  FFTW_MEASURE);
	
	this->planzbac = fftw_plan_many_dft(     1,    &this->nz,  this->ns,  
						(fftw_complex*) auxg,  &this->nz,  1,  this->nz,
						(fftw_complex*) auxg,  &this->nz,  1,  this->nz,  FFTW_BACKWARD,  FFTW_MEASURE);

	//---------------------------------------------------------
	//                              2 D - XY
	//---------------------------------------------------------
	//1D+1D is much faster than 2D FFT!
	//in-place fft is better for c2c and out-of-place fft is better for c2r
	int *embed = nullptr;
	int npy = this->nplane * this->ny;
	if(this->xprime)
	{
		this->planyfor  = fftw_plan_many_dft(  1, &this->ny,	this->nplane,	 (fftw_complex *)auxr, 		  embed, nplane,     1,
				(fftw_complex *)auxr, 	 embed, nplane,		1,		 FFTW_FORWARD,	FFTW_MEASURE   );
		this->planybac  = fftw_plan_many_dft(  1, &this->ny,	this->nplane,	 (fftw_complex *)auxr, 		  embed, nplane,     1,
				(fftw_complex *)auxr, 	 embed, nplane,		1,		 FFTW_BACKWARD,	FFTW_MEASURE   );
		if(this->gamma_only)
		{
			this->planxr2c  = fftw_plan_many_dft_r2c(  1, &this->nx, npy,	r_rspace , embed, npy,      1,
				(fftw_complex*)auxr, embed, npy,		1,	FFTW_MEASURE   );
			this->planxc2r  = fftw_plan_many_dft_c2r(  1, &this->nx, npy,	(fftw_complex*)auxr , embed, npy,      1,
				r_rspace, embed, npy,		1,			FFTW_MEASURE   );
		}
		else
		{
			this->planxfor1  = fftw_plan_many_dft(  1, &this->nx,	npy,	 (fftw_complex *)auxr, 		  embed, npy,     1,
				(fftw_complex *)auxr, 	 embed, npy,		1,		 FFTW_FORWARD,	FFTW_MEASURE   );
			this->planxbac1  = fftw_plan_many_dft(  1, &this->nx,	npy,	 (fftw_complex *)auxr, 		  embed, npy,     1,
				(fftw_complex *)auxr, 	 embed, npy,		1,		 FFTW_BACKWARD,	FFTW_MEASURE   );
		}
		
	}
	else
	{
		this->planxfor1  = fftw_plan_many_dft(  1, &this->nx,	this->nplane * (lixy + 1),	 (fftw_complex *)auxr, 		  embed, npy,     1,
				(fftw_complex *)auxr, 	 embed, npy,		1,		 FFTW_FORWARD,	FFTW_MEASURE   );
		this->planxbac1  = fftw_plan_many_dft(  1, &this->nx,	this->nplane * (lixy + 1),	 (fftw_complex *)auxr, 		  embed, npy,     1,
				(fftw_complex *)auxr, 	 embed, npy,		1,		 FFTW_BACKWARD,	FFTW_MEASURE   );
		if(this->gamma_only)
		{
			this->planyr2c  = fftw_plan_many_dft_r2c(  1, &this->ny, this->nplane,	r_rspace , embed, this->nplane,      1,
				(fftw_complex*)auxr, embed, this->nplane,		1,	FFTW_MEASURE   );
			this->planyc2r  = fftw_plan_many_dft_c2r(  1, &this->ny, this->nplane,	(fftw_complex*)auxr , embed, this->nplane,      1,
				r_rspace, embed, this->nplane,		1,			FFTW_MEASURE   );
		}
		else
		{

			this->planxfor2  = fftw_plan_many_dft(  1, &this->nx,	this->nplane * (ny - rixy),	 (fftw_complex *)auxr, 		  embed, npy,     1,
					(fftw_complex *)auxr, 	 embed, npy,		1,		 FFTW_FORWARD,	FFTW_MEASURE   );
			this->planxbac2  = fftw_plan_many_dft(  1, &this->nx,	this->nplane * (ny - rixy),	 (fftw_complex *)auxr, 		  embed, npy,     1,
					(fftw_complex *)auxr, 	 embed, npy,		1,		 FFTW_BACKWARD,	FFTW_MEASURE   );
			this->planyfor  = fftw_plan_many_dft(  1, &this->ny, this->nplane,	(fftw_complex*)auxr , embed, this->nplane,      1,
				(fftw_complex*)auxr, embed, this->nplane,		1,		 FFTW_FORWARD,	FFTW_MEASURE   );
			this->planybac  = fftw_plan_many_dft(  1, &this->ny, this->nplane,	(fftw_complex*)auxr , embed, this->nplane,      1,
				(fftw_complex*)auxr, embed, this->nplane,		1,		 FFTW_BACKWARD,	FFTW_MEASURE   );
		}
	}

    //---------------------------------------------------------
    //                              3 D - XYZ
    //---------------------------------------------------------
    //in-place fft test
    //this->plan3dforward = fftw_plan_dft_3d(
    //    this->nx, this->ny, this->nz,
    //    reinterpret_cast<fftw_complex *>(auxr_3d),
    //    reinterpret_cast<fftw_complex *>(auxr_3d),
    //    FFTW_FORWARD, FFTW_MEASURE);
    //this->plan3dbackward = fftw_plan_dft_3d(
    //    this->nx, this->ny, this->nz,
    //    reinterpret_cast<fftw_complex *>(auxr_3d),
    //    reinterpret_cast<fftw_complex *>(auxr_3d),
    //    FFTW_BACKWARD, FFTW_MEASURE);

    #if defined(__CUDA) || defined(__UT_USE_CUDA)
    if (GlobalV::device_flag == "gpu") {
        cufftPlan3d(&fft_handle, this->nx, this->ny, this->nz, CUFFT_Z2Z);
    }
    #endif

	destroyp = false;
}

#ifdef __MIX_PRECISION
void FFT :: initplanf()
{
	//---------------------------------------------------------
	//                              1 D
	//---------------------------------------------------------

	//               fftw_plan_many_dft(int rank,          const int *n,       int howmany,
	//					                fftw_complex *in,  const int *inembed, int istride, int idist, 
	//					                fftw_complex *out, const int *onembed, int ostride, int odist, int sign, unsigned flags);
	
	this->planfzfor = fftwf_plan_many_dft(     1,    &this->nz,  this->ns,  
					    (fftwf_complex*) auxfg,  &this->nz,  1,  this->nz,
					    (fftwf_complex*) auxfg,  &this->nz,  1,  this->nz,  FFTW_FORWARD,  FFTW_MEASURE);
	
	this->planfzbac = fftwf_plan_many_dft(     1,    &this->nz,  this->ns,  
						(fftwf_complex*) auxfg,  &this->nz,  1,  this->nz,
						(fftwf_complex*) auxfg,  &this->nz,  1,  this->nz,  FFTW_BACKWARD,  FFTW_MEASURE);
	//---------------------------------------------------------
	//                              2 D
	//---------------------------------------------------------
	
	int *embed = nullptr;
	int npy = this->nplane * this->ny;
	if(this->xprime)
	{
		this->planfyfor  = fftwf_plan_many_dft(  1, &this->ny,	this->nplane,	 (fftwf_complex *)auxfr, 		  embed, nplane,     1,
				(fftwf_complex *)auxfr, 	 embed, nplane,		1,		 FFTW_FORWARD,	FFTW_MEASURE   );
		this->planfybac  = fftwf_plan_many_dft(  1, &this->ny,	this->nplane,	 (fftwf_complex *)auxfr, 		  embed, nplane,     1,
				(fftwf_complex *)auxfr, 	 embed, nplane,		1,		 FFTW_BACKWARD,	FFTW_MEASURE   );
		if(this->gamma_only)
		{
			this->planfxr2c  = fftwf_plan_many_dft_r2c(  1, &this->nx, npy,	rf_rspace , embed, npy,      1,
				(fftwf_complex*)auxfr, embed, npy,		1,	FFTW_MEASURE   );
			this->planfxc2r  = fftwf_plan_many_dft_c2r(  1, &this->nx, npy,	(fftwf_complex*)auxfr , embed, npy,      1,
				rf_rspace, embed, npy,		1,			FFTW_MEASURE   );
		}
		else
		{
			this->planfxfor1  = fftwf_plan_many_dft(  1, &this->nx,	npy,	 (fftwf_complex *)auxfr, 		  embed, npy,     1,
				(fftwf_complex *)auxfr, 	 embed, npy,		1,		 FFTW_FORWARD,	FFTW_MEASURE   );
			this->planfxbac1  = fftwf_plan_many_dft(  1, &this->nx,	npy,	 (fftwf_complex *)auxfr, 		  embed, npy,     1,
				(fftwf_complex *)auxfr, 	 embed, npy,		1,		 FFTW_BACKWARD,	FFTW_MEASURE   );
		}
		
	}
	else
	{
		this->planfxfor1  = fftwf_plan_many_dft(  1, &this->nx,	this->nplane * (lixy + 1),	 (fftwf_complex *)auxfr, 		  embed, npy,     1,
				(fftwf_complex *)auxfr, 	 embed, npy,		1,		 FFTW_FORWARD,	FFTW_MEASURE   );
		this->planfxbac1  = fftwf_plan_many_dft(  1, &this->nx,	this->nplane * (lixy + 1),	 (fftwf_complex *)auxfr, 		  embed, npy,     1,
				(fftwf_complex *)auxfr, 	 embed, npy,		1,		 FFTW_BACKWARD,	FFTW_MEASURE   );
		if(this->gamma_only)
		{
			this->planfyr2c  = fftwf_plan_many_dft_r2c(  1, &this->ny, this->nplane,	rf_rspace , embed, this->nplane,      1,
				(fftwf_complex*)auxfr, embed, this->nplane,		1,	FFTW_MEASURE   );
			this->planfyc2r  = fftwf_plan_many_dft_c2r(  1, &this->ny, this->nplane,	(fftwf_complex*)auxfr , embed, this->nplane,      1,
				rf_rspace, embed, this->nplane,		1,			FFTW_MEASURE   );
		}
		else
		{
			this->planfxfor2  = fftwf_plan_many_dft(  1, &this->nx,	this->nplane * (this->ny - rixy),	 (fftwf_complex *)auxfr, 		  embed, npy,     1,
				(fftwf_complex *)auxfr, 	 embed, npy,		1,		 FFTW_FORWARD,	FFTW_MEASURE   );
			this->planfxbac2  = fftwf_plan_many_dft(  1, &this->nx,	this->nplane * (this->ny - rixy),	 (fftwf_complex *)auxfr, 		  embed, npy,     1,
				(fftwf_complex *)auxfr, 	 embed, npy,		1,		 FFTW_BACKWARD,	FFTW_MEASURE   );
			this->planfyfor  = fftwf_plan_many_dft(  1, &this->ny, this->nplane,	(fftwf_complex*)auxfr , embed, this->nplane,      1,
				(fftwf_complex*)auxfr, embed, this->nplane,		1,		 FFTW_FORWARD,	FFTW_MEASURE   );
			this->planfybac  = fftwf_plan_many_dft(  1, &this->ny, this->nplane,	(fftwf_complex*)auxfr , embed, this->nplane,      1,
				(fftwf_complex*)auxfr, embed, this->nplane,		1,		 FFTW_BACKWARD,	FFTW_MEASURE   );
		}
	}
	destroypf = false;
}
#endif

void FFT :: initplan_mpi()
{

}

#ifdef __MIX_PRECISION
void FFT :: initplanf_mpi()
{
	
}
#endif

void FFT:: cleanFFT()
{
	if(destroyp==true) return;
	fftw_destroy_plan(planzfor);
	fftw_destroy_plan(planzbac);
	if(this->xprime)
	{
		fftw_destroy_plan(planyfor);
		fftw_destroy_plan(planybac);
		if(this->gamma_only)
		{
			fftw_destroy_plan(planxr2c);
			fftw_destroy_plan(planxc2r);
		}
		else
		{
			fftw_destroy_plan(planxfor1);
			fftw_destroy_plan(planxbac1);
		}
	}
	else
	{
		fftw_destroy_plan(planxfor1);
		fftw_destroy_plan(planxbac1);
		if(this->gamma_only)
		{
			fftw_destroy_plan(planyr2c);
			fftw_destroy_plan(planyc2r);
		}
		else
		{
			fftw_destroy_plan(planxfor2);
			fftw_destroy_plan(planxbac2);
			fftw_destroy_plan(planyfor);
			fftw_destroy_plan(planybac);
		}
	}
    // fftw_destroy_plan(this->plan3dforward);
    // fftw_destroy_plan(this->plan3dbackward);
    #if defined(__CUDA) || defined(__UT_USE_CUDA)
    if (GlobalV::device_flag == "gpu") {
        cufftDestroy(fft_handle);
    }
    #endif
	destroyp = true;
	return;
}
#ifdef __MIX_PRECISION
void FFT:: cleanfFFT()
{
	if(destroypf==true) return;
	fftwf_destroy_plan(planfzfor);
	fftwf_destroy_plan(planfzbac);
	if(this->xprime)
	{
		fftwf_destroy_plan(planfyfor);
		fftwf_destroy_plan(planfybac);
		if(this->gamma_only)
		{
			fftwf_destroy_plan(planfxr2c);
			fftwf_destroy_plan(planfxc2r);
		}
		else
		{
			fftwf_destroy_plan(planfxfor1);
			fftwf_destroy_plan(planfxbac1);
		}
	}
	else
	{
		fftwf_destroy_plan(planfxfor1);
		fftwf_destroy_plan(planfxbac1);
		if(this->gamma_only)
		{
			fftwf_destroy_plan(planfyr2c);
			fftwf_destroy_plan(planfyc2r);
		}
		else
		{
			fftwf_destroy_plan(planfxfor2);
			fftwf_destroy_plan(planfxbac2);
			fftwf_destroy_plan(planfyfor);
			fftwf_destroy_plan(planfybac);
		}
	}
	destroypf = true;
	return;
}
#endif


void FFT::fftzfor(std::complex<double>* & in, std::complex<double>* & out)
{
	fftw_execute_dft(this->planzfor,(fftw_complex *)in,(fftw_complex *)out);
	return;
}

void FFT::fftzbac(std::complex<double>* & in, std::complex<double>* & out)
{
	fftw_execute_dft(this->planzbac,(fftw_complex *)in, (fftw_complex *)out);
	return;
}

void FFT::fftxyfor(std::complex<double>* & in, std::complex<double>* & out)
{
	int npy = this->nplane * this-> ny;
	if(this->xprime)
	{
		fftw_execute_dft( this->planxfor1, (fftw_complex *)in, (fftw_complex *)out);

		for(int i = 0 ; i < this->lixy + 1; ++i)
		{
			fftw_execute_dft( this->planyfor, (fftw_complex *)&in[i*npy], (fftw_complex *)&out[i*npy]);
		}
		for(int i = rixy ; i < this->nx; ++i)
		{
			fftw_execute_dft( this->planyfor, (fftw_complex *)&in[i*npy], (fftw_complex *)&out[i*npy]);
		}
	}
	else
	{
		for (int i=0; i<this->nx;++i)
		{
			fftw_execute_dft( this->planyfor, (fftw_complex *)&in[i*npy], (fftw_complex *)&out[i*npy]);
		}

		fftw_execute_dft( this->planxfor1, (fftw_complex *)in, (fftw_complex *)out);
		fftw_execute_dft( this->planxfor2, (fftw_complex *)&in[rixy*nplane], (fftw_complex *)&out[rixy*nplane]);
	}
	return;
}

void FFT::fftxybac(std::complex<double>* & in, std::complex<double>* & out)
{
	int npy = this->nplane * this-> ny;
	if(this->xprime)
	{
		for(int i = 0 ; i < this->lixy + 1; ++i)
		{
			fftw_execute_dft( this->planybac, (fftw_complex*)&in[i*npy], (fftw_complex*)&out[i*npy] );
		}
		for(int i = rixy ; i < this->nx; ++i)
		{
			fftw_execute_dft( this->planybac, (fftw_complex*)&in[i*npy], (fftw_complex*)&out[i*npy] );
		}

		fftw_execute_dft( this->planxbac1, (fftw_complex *)in, (fftw_complex *)out);
	}
	else
	{
		fftw_execute_dft( this->planxbac1, (fftw_complex *)in, (fftw_complex *)out);
		fftw_execute_dft( this->planxbac2, (fftw_complex *)&in[rixy*nplane], (fftw_complex *)&out[rixy*nplane]);

		for (int i=0; i<this->nx;++i)
		{
			fftw_execute_dft( this->planybac, (fftw_complex*)&in[i*npy], (fftw_complex*)&out[i*npy] );
		}
	}
	return;
}

void FFT::fftxyr2c(double* &in, std::complex<double>* & out)
{
	int npy = this->nplane * this-> ny;
	if(this->xprime)
	{
		fftw_execute_dft_r2c( this->planxr2c, in, (fftw_complex *)out);

		for(int i = 0 ; i < this->lixy + 1; ++i)
		{
			fftw_execute_dft( this->planyfor, (fftw_complex *)&out[i*npy], (fftw_complex *)&out[i*npy]);
		}
	}
	else
	{
		for (int i=0; i<this->nx;++i)
		{
			fftw_execute_dft_r2c( this->planyr2c, &in[i*npy], (fftw_complex*)&out[i*npy] );
		}

		fftw_execute_dft( this->planxfor1, (fftw_complex *)out, (fftw_complex *)out);
	}
	return;
}


void FFT::fftxyc2r(std::complex<double>* & in, double* & out)
{
	int npy = this->nplane * this-> ny;
	if(this->xprime)
	{
		for(int i = 0 ; i < this->lixy + 1; ++i)
		{
			fftw_execute_dft( this->planybac, (fftw_complex *)&in[i*npy], (fftw_complex *)&in[i*npy]);
		}

		fftw_execute_dft_c2r( this->planxc2r, (fftw_complex *)in, out);
	}
	else
	{
		fftw_execute_dft( this->planxbac1, (fftw_complex *)in, (fftw_complex *)in);

		for (int i=0; i<this->nx;++i)
		{
			fftw_execute_dft_c2r( this->planyc2r, (fftw_complex*)&in[i*npy], &out[i*npy] );
		}
	}
	return;
}

#if defined(__CUDA) || defined(__UT_USE_CUDA)
void FFT::fft3D_forward(std::complex<double>* & in, std::complex<double>* & out)
{
//    fftw_execute_dft(
//        this->plan3dforward,
//        reinterpret_cast<fftw_complex *>(in),
//        reinterpret_cast<fftw_complex *>(out));
    cufftExecZ2Z(this->fft_handle,
          reinterpret_cast<cufftDoubleComplex*>(in),
          reinterpret_cast<cufftDoubleComplex*>(out),
          CUFFT_FORWARD);
    cudaDeviceSynchronize();
}

void FFT::fft3D_backward(std::complex<double>* & in, std::complex<double>* & out)
{
//    fftw_execute_dft(
//        this->plan3dbackward,
//        reinterpret_cast<fftw_complex *>(in),
//        reinterpret_cast<fftw_complex *>(out));
    cufftExecZ2Z(this->fft_handle,
             reinterpret_cast<cufftDoubleComplex*>(in),
             reinterpret_cast<cufftDoubleComplex*>(out),
             CUFFT_INVERSE);
    cudaDeviceSynchronize();
}
#endif


#ifdef __MIX_PRECISION
void FFT::fftfzfor(std::complex<float>* & in, std::complex<float>* & out)
{
	fftwf_execute_dft(this->planfzfor,(fftwf_complex *)in,(fftwf_complex *)out);
	return;
}

void FFT::fftfzbac(std::complex<float>* & in, std::complex<float>* & out)
{
	fftwf_execute_dft(this->planfzbac,(fftwf_complex *)in, (fftwf_complex *)out);
	return;
}

void FFT::fftfxyfor(std::complex<float>* & in, std::complex<float>* & out)
{
	int npy = this->nplane * this-> ny;
	if(this->xprime)
	{
		fftwf_execute_dft( this->planfxfor1, (fftwf_complex *)in, (fftwf_complex *)out);

		for(int i = 0 ; i < this->lixy + 1; ++i)
		{
			fftwf_execute_dft( this->planfyfor, (fftwf_complex *)&in[i*npy], (fftwf_complex *)&out[i*npy]);
		}
		for(int i = rixy ; i < this->nx; ++i)
		{
			fftwf_execute_dft( this->planfyfor, (fftwf_complex *)&in[i*npy], (fftwf_complex *)&out[i*npy]);
		}
	}
	else
	{
		for (int i=0; i<this->nx;++i)
		{
			fftwf_execute_dft( this->planfyfor, (fftwf_complex *)&in[i*npy], (fftwf_complex *)&out[i*npy]);
		}

		fftwf_execute_dft( this->planfxfor1, (fftwf_complex *)in, (fftwf_complex *)out);
		fftwf_execute_dft( this->planfxfor2, (fftwf_complex *)&in[rixy*nplane], (fftwf_complex *)&out[rixy*nplane]);
	}
	return;
}

void FFT::fftfxybac(std::complex<float>* & in, std::complex<float>* & out)
{
	int npy = this->nplane * this-> ny;
	if(this->xprime)
	{
		for(int i = 0 ; i < this->lixy + 1; ++i)
		{
			fftwf_execute_dft( this->planfybac, (fftwf_complex*)&in[i*npy], (fftwf_complex*)&out[i*npy] );
		}
		for(int i = rixy ; i < this->nx; ++i)
		{
			fftwf_execute_dft( this->planfybac, (fftwf_complex*)&in[i*npy], (fftwf_complex*)&out[i*npy] );
		}

		fftwf_execute_dft( this->planfxbac1, (fftwf_complex *)in, (fftwf_complex *)out);
	}
	else
	{
		fftwf_execute_dft( this->planfxbac1, (fftwf_complex *)in, (fftwf_complex *)out);
		fftwf_execute_dft( this->planfxbac2, (fftwf_complex *)&in[rixy*nplane], (fftwf_complex *)&out[rixy*nplane]);

		for (int i=0; i<this->nx;++i)
		{
			fftwf_execute_dft( this->planfybac, (fftwf_complex*)&in[i*npy], (fftwf_complex*)&out[i*npy] );
		}
	}
	return;
}

void FFT::fftfxyr2c(float* &in, std::complex<float>* & out)
{
	int npy = this->nplane * this-> ny;
	if(this->xprime)
	{
		fftwf_execute_dft_r2c( this->planfxr2c, in, (fftwf_complex *)out);

		for(int i = 0 ; i < this->lixy + 1; ++i)
		{
			fftwf_execute_dft( this->planfyfor, (fftwf_complex *)&out[i*npy], (fftwf_complex *)&out[i*npy]);
		}
	}
	else
	{
		for (int i=0; i<this->nx;++i)
		{
			fftwf_execute_dft_r2c( this->planfyr2c, &in[i*npy], (fftwf_complex*)&out[i*npy] );
		}

		fftwf_execute_dft( this->planfxfor1, (fftwf_complex *)out, (fftwf_complex *)out);
	}
	
	return;
}


void FFT::fftfxyc2r(std::complex<float>* & in, float* & out)
{
	int npy = this->nplane * this-> ny;
	if(this->xprime)
	{
		for(int i = 0 ; i < this->lixy + 1; ++i)
		{
			fftwf_execute_dft( this->planfybac, (fftwf_complex *)&in[i*npy], (fftwf_complex *)&in[i*npy]);
		}

		fftwf_execute_dft_c2r( this->planfxc2r, (fftwf_complex *)in, out);
	}
	else
	{
		fftwf_execute_dft( this->planfxbac1, (fftwf_complex *)in, (fftwf_complex *)in);

		for (int i=0; i<this->nx;++i)
		{
			fftwf_execute_dft_c2r( this->planfyc2r, (fftwf_complex*)&in[i*npy], &out[i*npy] );
		}
	}
	return;
}
#endif
}
