#include "fft.h"

#include "module_base/global_variable.h"
#include "module_base/memory.h"

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
	if(z_auxg!=nullptr) {fftw_free(z_auxg); z_auxg = nullptr;}
	if(z_auxr!=nullptr) {fftw_free(z_auxr); z_auxr = nullptr;}
	d_rspace = nullptr;
#if defined(__CUDA) || defined(__ROCM)
    if (GlobalV::device_flag == "gpu") {
        if (GlobalV::precision_flag == "single") {
            if (c_auxr_3d != nullptr) {
                delmem_cd_op()(gpu_ctx, c_auxr_3d);
                c_auxr_3d = nullptr;
            }
        }
        else {
            if (z_auxr_3d != nullptr) {
                delmem_zd_op()(gpu_ctx, z_auxr_3d);
                z_auxr_3d = nullptr;
            }
        }
    }
#endif
    if (GlobalV::precision_flag == "single") {
        this->cleanfFFT();
        if (c_auxg != nullptr) {
            fftw_free(c_auxg);
            c_auxg = nullptr;
        }
        if (c_auxr != nullptr) {
            fftw_free(c_auxr);
            c_auxr = nullptr;
        }
        s_rspace = nullptr;
    }
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
		z_auxg  = (std::complex<double> *) fftw_malloc(sizeof(fftw_complex) * maxgrids);
		z_auxr  = (std::complex<double> *) fftw_malloc(sizeof(fftw_complex) * maxgrids);
		ModuleBase::Memory::record("FFT::grid", 2 * sizeof(fftw_complex) * maxgrids);
		d_rspace = (double *) z_auxg;
        // auxr_3d = static_cast<std::complex<double> *>(
        //     fftw_malloc(sizeof(fftw_complex) * (this->nx * this->ny * this->nz)));
#if defined(__CUDA) || defined(__ROCM)
        if (GlobalV::device_flag == "gpu") {
            if (GlobalV::precision_flag == "single") {
                resmem_cd_op()(gpu_ctx, this->c_auxr_3d, this->nx * this->ny * this->nz);
            }
            else {
                resmem_zd_op()(gpu_ctx, this->z_auxr_3d, this->nx * this->ny * this->nz);
            }
        }
#endif
        if (GlobalV::precision_flag == "single") {
            c_auxg  = (std::complex<float> *) fftw_malloc(sizeof(fftwf_complex) * maxgrids);
            c_auxr  = (std::complex<float> *) fftw_malloc(sizeof(fftwf_complex) * maxgrids);
			ModuleBase::Memory::record("FFT::grid_s", 2 * sizeof(fftwf_complex) * maxgrids);
            s_rspace = (float *) c_auxg;
        }
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
        if (GlobalV::precision_flag == "single") {
            this->initplanf();
        }
	}
#if defined(__FFTW3_MPI) && defined(__MPI)
	else
	{
		// this->initplan_mpi();
        // if (GlobalV::precision_flag == "single") {
		//     this->initplanf_mpi();
        // }
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
					    (fftw_complex*) z_auxg,  &this->nz,  1,  this->nz,
					    (fftw_complex*) z_auxg,  &this->nz,  1,  this->nz,  FFTW_FORWARD,  FFTW_MEASURE);
	
	this->planzbac = fftw_plan_many_dft(     1,    &this->nz,  this->ns,  
						(fftw_complex*) z_auxg,  &this->nz,  1,  this->nz,
						(fftw_complex*) z_auxg,  &this->nz,  1,  this->nz,  FFTW_BACKWARD,  FFTW_MEASURE);

	//---------------------------------------------------------
	//                              2 D - XY
	//---------------------------------------------------------
	//1D+1D is much faster than 2D FFT!
	//in-place fft is better for c2c and out-of-place fft is better for c2r
	int *embed = nullptr;
	int npy = this->nplane * this->ny;
	if(this->xprime)
	{
		this->planyfor  = fftw_plan_many_dft(  1, &this->ny,	this->nplane,	 (fftw_complex *)z_auxr, 		  embed, nplane,     1,
				(fftw_complex *)z_auxr, 	 embed, nplane,		1,		 FFTW_FORWARD,	FFTW_MEASURE   );
		this->planybac  = fftw_plan_many_dft(  1, &this->ny,	this->nplane,	 (fftw_complex *)z_auxr, 		  embed, nplane,     1,
				(fftw_complex *)z_auxr, 	 embed, nplane,		1,		 FFTW_BACKWARD,	FFTW_MEASURE   );
		if(this->gamma_only)
		{
			this->planxr2c  = fftw_plan_many_dft_r2c(  1, &this->nx, npy,	d_rspace , embed, npy,      1,
				(fftw_complex*)z_auxr, embed, npy,		1,	FFTW_MEASURE   );
			this->planxc2r  = fftw_plan_many_dft_c2r(  1, &this->nx, npy,	(fftw_complex*)z_auxr , embed, npy,      1,
				d_rspace, embed, npy,		1,			FFTW_MEASURE   );
		}
		else
		{
			this->planxfor1  = fftw_plan_many_dft(  1, &this->nx,	npy,	 (fftw_complex *)z_auxr, 		  embed, npy,     1,
				(fftw_complex *)z_auxr, 	 embed, npy,		1,		 FFTW_FORWARD,	FFTW_MEASURE   );
			this->planxbac1  = fftw_plan_many_dft(  1, &this->nx,	npy,	 (fftw_complex *)z_auxr, 		  embed, npy,     1,
				(fftw_complex *)z_auxr, 	 embed, npy,		1,		 FFTW_BACKWARD,	FFTW_MEASURE   );
		}
		
	}
	else
	{
		this->planxfor1  = fftw_plan_many_dft(  1, &this->nx,	this->nplane * (lixy + 1),	 (fftw_complex *)z_auxr, 		  embed, npy,     1,
				(fftw_complex *)z_auxr, 	 embed, npy,		1,		 FFTW_FORWARD,	FFTW_MEASURE   );
		this->planxbac1  = fftw_plan_many_dft(  1, &this->nx,	this->nplane * (lixy + 1),	 (fftw_complex *)z_auxr, 		  embed, npy,     1,
				(fftw_complex *)z_auxr, 	 embed, npy,		1,		 FFTW_BACKWARD,	FFTW_MEASURE   );
		if(this->gamma_only)
		{
			this->planyr2c  = fftw_plan_many_dft_r2c(  1, &this->ny, this->nplane,	d_rspace , embed, this->nplane,      1,
				(fftw_complex*)z_auxr, embed, this->nplane,		1,	FFTW_MEASURE   );
			this->planyc2r  = fftw_plan_many_dft_c2r(  1, &this->ny, this->nplane,	(fftw_complex*)z_auxr , embed, this->nplane,      1,
				d_rspace, embed, this->nplane,		1,			FFTW_MEASURE   );
		}
		else
		{

			this->planxfor2  = fftw_plan_many_dft(  1, &this->nx,	this->nplane * (ny - rixy),	 (fftw_complex *)z_auxr, 		  embed, npy,     1,
					(fftw_complex *)z_auxr, 	 embed, npy,		1,		 FFTW_FORWARD,	FFTW_MEASURE   );
			this->planxbac2  = fftw_plan_many_dft(  1, &this->nx,	this->nplane * (ny - rixy),	 (fftw_complex *)z_auxr, 		  embed, npy,     1,
					(fftw_complex *)z_auxr, 	 embed, npy,		1,		 FFTW_BACKWARD,	FFTW_MEASURE   );
			this->planyfor  = fftw_plan_many_dft(  1, &this->ny, this->nplane,	(fftw_complex*)z_auxr , embed, this->nplane,      1,
				(fftw_complex*)z_auxr, embed, this->nplane,		1,		 FFTW_FORWARD,	FFTW_MEASURE   );
			this->planybac  = fftw_plan_many_dft(  1, &this->ny, this->nplane,	(fftw_complex*)z_auxr , embed, this->nplane,      1,
				(fftw_complex*)z_auxr, embed, this->nplane,		1,		 FFTW_BACKWARD,	FFTW_MEASURE   );
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

#if defined(__CUDA) || defined(__ROCM)
    if (GlobalV::device_flag == "gpu") {
        if (GlobalV::precision_flag == "single") {
        #if defined(__CUDA)
            cufftPlan3d(&c_handle, this->nx, this->ny, this->nz, CUFFT_C2C);
        #elif defined(__ROCM)
            hipfftPlan3d(&c_handle, this->nx, this->ny, this->nz, HIPFFT_C2C);
        #endif
        }
        else {
        #if defined(__CUDA)
            cufftPlan3d(&z_handle, this->nx, this->ny, this->nz, CUFFT_Z2Z);
        #elif defined(__ROCM)
            hipfftPlan3d(&z_handle, this->nx, this->ny, this->nz, HIPFFT_Z2Z);
        #endif
        }
    }
#endif

	destroyp = false;
}

void FFT :: initplanf()
{
	//---------------------------------------------------------
	//                              1 D
	//---------------------------------------------------------

	//               fftw_plan_many_dft(int rank,          const int *n,       int howmany,
	//					                fftw_complex *in,  const int *inembed, int istride, int idist, 
	//					                fftw_complex *out, const int *onembed, int ostride, int odist, int sign, unsigned flags);
	
	this->planfzfor = fftwf_plan_many_dft(     1,    &this->nz,  this->ns,  
					    (fftwf_complex*) c_auxg,  &this->nz,  1,  this->nz,
					    (fftwf_complex*) c_auxg,  &this->nz,  1,  this->nz,  FFTW_FORWARD,  FFTW_MEASURE);
	
	this->planfzbac = fftwf_plan_many_dft(     1,    &this->nz,  this->ns,  
						(fftwf_complex*) c_auxg,  &this->nz,  1,  this->nz,
						(fftwf_complex*) c_auxg,  &this->nz,  1,  this->nz,  FFTW_BACKWARD,  FFTW_MEASURE);
	//---------------------------------------------------------
	//                              2 D
	//---------------------------------------------------------
	
	int *embed = nullptr;
	int npy = this->nplane * this->ny;
	if(this->xprime)
	{
		this->planfyfor  = fftwf_plan_many_dft(  1, &this->ny,	this->nplane,	 (fftwf_complex *)c_auxr, 		  embed, nplane,     1,
				(fftwf_complex *)c_auxr, 	 embed, nplane,		1,		 FFTW_FORWARD,	FFTW_MEASURE   );
		this->planfybac  = fftwf_plan_many_dft(  1, &this->ny,	this->nplane,	 (fftwf_complex *)c_auxr, 		  embed, nplane,     1,
				(fftwf_complex *)c_auxr, 	 embed, nplane,		1,		 FFTW_BACKWARD,	FFTW_MEASURE   );
		if(this->gamma_only)
		{
			this->planfxr2c  = fftwf_plan_many_dft_r2c(  1, &this->nx, npy,	s_rspace , embed, npy,      1,
				(fftwf_complex*)c_auxr, embed, npy,		1,	FFTW_MEASURE   );
			this->planfxc2r  = fftwf_plan_many_dft_c2r(  1, &this->nx, npy,	(fftwf_complex*)c_auxr , embed, npy,      1,
				s_rspace, embed, npy,		1,			FFTW_MEASURE   );
		}
		else
		{
			this->planfxfor1  = fftwf_plan_many_dft(  1, &this->nx,	npy,	 (fftwf_complex *)c_auxr, 		  embed, npy,     1,
				(fftwf_complex *)c_auxr, 	 embed, npy,		1,		 FFTW_FORWARD,	FFTW_MEASURE   );
			this->planfxbac1  = fftwf_plan_many_dft(  1, &this->nx,	npy,	 (fftwf_complex *)c_auxr, 		  embed, npy,     1,
				(fftwf_complex *)c_auxr, 	 embed, npy,		1,		 FFTW_BACKWARD,	FFTW_MEASURE   );
		}
		
	}
	else
	{
		this->planfxfor1  = fftwf_plan_many_dft(  1, &this->nx,	this->nplane * (lixy + 1),	 (fftwf_complex *)c_auxr, 		  embed, npy,     1,
				(fftwf_complex *)c_auxr, 	 embed, npy,		1,		 FFTW_FORWARD,	FFTW_MEASURE   );
		this->planfxbac1  = fftwf_plan_many_dft(  1, &this->nx,	this->nplane * (lixy + 1),	 (fftwf_complex *)c_auxr, 		  embed, npy,     1,
				(fftwf_complex *)c_auxr, 	 embed, npy,		1,		 FFTW_BACKWARD,	FFTW_MEASURE   );
		if(this->gamma_only)
		{
			this->planfyr2c  = fftwf_plan_many_dft_r2c(  1, &this->ny, this->nplane,	s_rspace , embed, this->nplane,      1,
				(fftwf_complex*)c_auxr, embed, this->nplane,		1,	FFTW_MEASURE   );
			this->planfyc2r  = fftwf_plan_many_dft_c2r(  1, &this->ny, this->nplane,	(fftwf_complex*)c_auxr , embed, this->nplane,      1,
				s_rspace, embed, this->nplane,		1,			FFTW_MEASURE   );
		}
		else
		{
			this->planfxfor2  = fftwf_plan_many_dft(  1, &this->nx,	this->nplane * (this->ny - rixy),	 (fftwf_complex *)c_auxr, 		  embed, npy,     1,
				(fftwf_complex *)c_auxr, 	 embed, npy,		1,		 FFTW_FORWARD,	FFTW_MEASURE   );
			this->planfxbac2  = fftwf_plan_many_dft(  1, &this->nx,	this->nplane * (this->ny - rixy),	 (fftwf_complex *)c_auxr, 		  embed, npy,     1,
				(fftwf_complex *)c_auxr, 	 embed, npy,		1,		 FFTW_BACKWARD,	FFTW_MEASURE   );
			this->planfyfor  = fftwf_plan_many_dft(  1, &this->ny, this->nplane,	(fftwf_complex*)c_auxr , embed, this->nplane,      1,
				(fftwf_complex*)c_auxr, embed, this->nplane,		1,		 FFTW_FORWARD,	FFTW_MEASURE   );
			this->planfybac  = fftwf_plan_many_dft(  1, &this->ny, this->nplane,	(fftwf_complex*)c_auxr , embed, this->nplane,      1,
				(fftwf_complex*)c_auxr, embed, this->nplane,		1,		 FFTW_BACKWARD,	FFTW_MEASURE   );
		}
	}
	destroypf = false;
}

// void FFT :: initplan_mpi()
// {

// }

// void FFT :: initplanf_mpi()
// {
	
// }

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
#if defined(__CUDA) || defined(__ROCM)
    if (GlobalV::device_flag == "gpu") {
        if (GlobalV::precision_flag == "single") {
        #if defined(__CUDA)
            cufftDestroy(c_handle);
        #elif defined(__ROCM)
            hipfftDestroy(c_handle);
        #endif
        }
        else {
        #if defined(__CUDA)
            cufftDestroy(z_handle);
        #elif defined(__ROCM)
            hipfftDestroy(z_handle);
        #endif
        }
    }
#endif
	destroyp = true;
}

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

template <>
void FFT::fftzfor(std::complex<float> * in, std::complex<float> * out)
{
    fftwf_execute_dft(this->planfzfor,(fftwf_complex *)in,(fftwf_complex *)out);
}
template <>
void FFT::fftzfor(std::complex<double> * in, std::complex<double> * out)
{
	fftw_execute_dft(this->planzfor,(fftw_complex *)in,(fftw_complex *)out);
}

template <>
void FFT::fftzbac(std::complex<float> * in, std::complex<float> * out)
{
    fftwf_execute_dft(this->planfzbac,(fftwf_complex *)in, (fftwf_complex *)out);
}
template <>
void FFT::fftzbac(std::complex<double> * in, std::complex<double> * out)
{
	fftw_execute_dft(this->planzbac,(fftw_complex *)in, (fftw_complex *)out);
}

template <>
void FFT::fftxyfor(std::complex<float> * in, std::complex<float> * out)
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
}
template <>
void FFT::fftxyfor(std::complex<double> * in, std::complex<double> * out)
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
}

template <>
void FFT::fftxybac(std::complex<float> * in, std::complex<float> * out)
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
}
template <>
void FFT::fftxybac(std::complex<double> * in, std::complex<double> * out)
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
}

template <>
void FFT::fftxyr2c(float * in, std::complex<float> * out)
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
}
template <>
void FFT::fftxyr2c(double * in, std::complex<double> * out)
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
}

template <>
void FFT::fftxyc2r(std::complex<float>* in, float* out)
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
}
template <>
void FFT::fftxyc2r(std::complex<double> * in, double * out)
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
}

#if defined(__CUDA) || defined(__ROCM)
template <>
void FFT::fft3D_forward(const psi::DEVICE_GPU * /*ctx*/, std::complex<float> * in, std::complex<float> * out)
{
#if defined(__CUDA)
    cufftExecC2C(this->c_handle,
          reinterpret_cast<cufftComplex*>(in),
          reinterpret_cast<cufftComplex*>(out),
          CUFFT_FORWARD);
    cudaDeviceSynchronize();
#elif defined(__ROCM)
    hipfftExecC2C(this->c_handle,
          reinterpret_cast<hipfftComplex*>(in),
          reinterpret_cast<hipfftComplex*>(out),
          HIPFFT_FORWARD);
    hipDeviceSynchronize();
#endif
}
template <>
void FFT::fft3D_forward(const psi::DEVICE_GPU * /*ctx*/, std::complex<double> * in, std::complex<double> * out)
{
#if defined(__CUDA)
    cufftExecZ2Z(this->z_handle,
          reinterpret_cast<cufftDoubleComplex*>(in),
          reinterpret_cast<cufftDoubleComplex*>(out),
          CUFFT_FORWARD);
    cudaDeviceSynchronize();
#elif defined(__ROCM)
    hipfftExecZ2Z(this->z_handle,
          reinterpret_cast<hipfftDoubleComplex*>(in),
          reinterpret_cast<hipfftDoubleComplex*>(out),
          HIPFFT_FORWARD);
    hipDeviceSynchronize();
#endif
}

template <>
void FFT::fft3D_backward(const psi::DEVICE_GPU * /*ctx*/, std::complex<float>* in, std::complex<float>* out)
{
#if defined(__CUDA)
    cufftExecC2C(this->c_handle,
             reinterpret_cast<cufftComplex*>(in),
             reinterpret_cast<cufftComplex*>(out),
             CUFFT_INVERSE);
    cudaDeviceSynchronize();
#elif defined(__ROCM)
    hipfftExecC2C(this->c_handle,
             reinterpret_cast<hipfftComplex*>(in),
             reinterpret_cast<hipfftComplex*>(out),
             HIPFFT_BACKWARD);
    hipDeviceSynchronize();
#endif
}
template <>
void FFT::fft3D_backward(const psi::DEVICE_GPU * /*ctx*/, std::complex<double>* in, std::complex<double>* out)
{
#if defined(__CUDA)
    cufftExecZ2Z(this->z_handle,
             reinterpret_cast<cufftDoubleComplex*>(in),
             reinterpret_cast<cufftDoubleComplex*>(out),
             CUFFT_INVERSE);
    cudaDeviceSynchronize();
#elif defined(__ROCM)
    hipfftExecZ2Z(this->z_handle,
             reinterpret_cast<hipfftDoubleComplex*>(in),
             reinterpret_cast<hipfftDoubleComplex*>(out),
             HIPFFT_BACKWARD);
    hipDeviceSynchronize();
#endif
}
#endif


template <>
float * FFT::get_rspace_data() {
    return this->s_rspace;
}
template <>
double * FFT::get_rspace_data() {
    return this->d_rspace;
}

template <>
std::complex<float> * FFT::get_auxr_data() {
    return this->c_auxr;
}
template <>
std::complex<double> * FFT::get_auxr_data() {
    return this->z_auxr;
}

template <>
std::complex<float> * FFT::get_auxg_data() {
    return this->c_auxg;
}
template <>
std::complex<double> * FFT::get_auxg_data() {
    return this->z_auxg;
}

#if defined(__CUDA) || defined(__ROCM)
template <>
std::complex<float> * FFT::get_auxr_3d_data() {
    return this->c_auxr_3d;
}
template <>
std::complex<double> * FFT::get_auxr_3d_data() {
    return this->z_auxr_3d;
}
#endif

}
