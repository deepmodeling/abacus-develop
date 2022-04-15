#include "fft.h"
namespace ModulePW
{

FFT::FFT()
{
	bigny = nx = ny = nz = 2;
	bignxy = nxy = 4;
	ns = 1; 
	nplane = 1;
	mpifft = false; 
	destroyp = true;
	gamma_only = false;
	aux2 = aux1 = NULL;
	r_rspace = NULL;
#ifdef __MIX_PRECISION
	destroypf = true;
	auxf2 = auxf1 = NULL;
	rf_rspace = NULL;
#endif
}

FFT::~FFT()
{
	this->cleanFFT();
	if(aux1!=NULL) fftw_free(aux1);
	if(aux2!=NULL) fftw_free(aux2);
#ifdef __MIX_PRECISION
	this->cleanfFFT();
	if(auxf1!=NULL) fftw_free(auxf1);
	if(auxf2!=NULL) fftw_free(auxf2);
#endif
}

void FFT:: initfft(int nx_in, int bigny_in, int nz_in, int liy_in, int riy_in, int ns_in, int nplane_in, int nproc_in, bool gamma_only_in, bool mpifft_in)
{
	this->gamma_only = gamma_only_in;
	this->nx = nx_in;
	this->bigny = bigny_in;
	if(this->gamma_only) this->ny = int(bigny/2) +1;
	else	this->ny = this->bigny;
	this->nz = nz_in;
	this->ns = ns_in;
	this->liy = liy_in;
	this->riy = riy_in;
	this->nplane = nplane_in;
	this->nproc = nproc_in;
	this->mpifft = mpifft_in;
	this->nxy = this->nx * this-> ny;
	this->bignxy = this->nx * this->bigny;
	this->maxgrids = (this->nz * this->ns > this->bignxy * nplane) ? this->nz * this->ns : this->bignxy * nplane;
	if(!this->mpifft)
	{
		aux1  = (std::complex<double> *) fftw_malloc(sizeof(fftw_complex) * maxgrids);
		aux2  = (std::complex<double> *) fftw_malloc(sizeof(fftw_complex) * maxgrids);
		r_rspace = (double *) aux1;
#ifdef __MIX_PRECISION
		auxf1  = (std::complex<float> *) fftw_malloc(sizeof(fftwf_complex) * maxgrids);
		auxf2  = (std::complex<float> *) fftw_malloc(sizeof(fftwf_complex) * maxgrids);
		rf_rspace = (float *) auxf1;
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
	//                              1 D
	//---------------------------------------------------------

	//               fftw_plan_many_dft(int rank,          const int *n,       int howmany,
	//					                fftw_complex *in,  const int *inembed, int istride, int idist, 
	//					                fftw_complex *out, const int *onembed, int ostride, int odist, int sign, unsigned flags);
	
	//It is better to use out-of-place fft for stride = 1.
	this->planzfor = fftw_plan_many_dft(     1,    &this->nz,  this->ns,  
					    (fftw_complex*) aux1,  &this->nz,  1,  this->nz,
					    (fftw_complex*) aux2,  &this->nz,  1,  this->nz,  FFTW_FORWARD,  FFTW_MEASURE);
	
	this->planzbac = fftw_plan_many_dft(     1,    &this->nz,  this->ns,  
						(fftw_complex*) aux1,  &this->nz,  1,  this->nz,
						(fftw_complex*) aux2,  &this->nz,  1,  this->nz,  FFTW_BACKWARD,  FFTW_MEASURE);

	// this->planzfor = fftw_plan_dft_1d(this->nz,(fftw_complex*) aux1,(fftw_complex*) aux1, FFTW_FORWARD,  FFTW_MEASURE);
	// this->planzbac = fftw_plan_dft_1d(this->nz,(fftw_complex*) aux1,(fftw_complex*) aux1,FFTW_BACKWARD,  FFTW_MEASURE);

	//---------------------------------------------------------
	//                              2 D
	//---------------------------------------------------------
	
	//int nrank[2] = {this->nx,this->bigny};
	int *embed = NULL;
	// It seems 1D+1D is much faster than 2D FFT!
	if(this->gamma_only)
	{
		// int padnpy = this->nplane * this->ny * 2;
		// int rankc[2] = {this->nx, this->padnpy};
		// int rankd[2] = {this->nx, this->padnpy*2};
		// // It seems 1D+1D is much faster than 2D FFT!
		// this->plan2r2c = fftw_plan_many_dft_r2c(   2,   nrank,  this->nplane,  
		// 									r_rspace,   rankd,  this->nplane,   1,
		// 					(fftw_complex*) aux2,   rankc,  this->nplane,   1,  FFTW_MEASURE);

		// this->plan2c2r = fftw_plan_many_dft_c2r(   2,   nrank,  this->nplane,  
		// 					(fftw_complex*) aux2,   rankc,  this->nplane,   1,
		// 									r_rspace,   rankd,  this->nplane,   1,  FFTW_MEASURE);

		// int npy = this->nplane * this->ny;
		int bignpy = this->nplane * this->bigny;
		this->planxfor1  = fftw_plan_many_dft(  1, &this->nx,	this->nplane * (liy + 1),	 (fftw_complex *)aux2, 		  embed, bignpy,     1,
				(fftw_complex *)aux2, 	 embed, bignpy,		1,		 FFTW_FORWARD,	FFTW_MEASURE   );
		this->planxbac1  = fftw_plan_many_dft(  1, &this->nx,	this->nplane * (liy + 1),	 (fftw_complex *)aux2, 		  embed, bignpy,     1,
				(fftw_complex *)aux2, 	 embed, bignpy,		1,		 FFTW_BACKWARD,	FFTW_MEASURE   );
		this->planyr2c  = fftw_plan_many_dft_r2c(  1, &this->bigny, this->nplane,	r_rspace , embed, this->nplane,      1,
			(fftw_complex*)aux1, embed, this->nplane,		1,	FFTW_MEASURE   );
		this->planyc2r  = fftw_plan_many_dft_c2r(  1, &this->bigny, this->nplane,	(fftw_complex*)aux1 , embed, this->nplane,      1,
			r_rspace, embed, this->nplane,		1,			FFTW_MEASURE   );
		
		// int padnpy = npy * 2;
		// this->planxfor  = fftw_plan_many_dft(  1, &this->nx,	npy,	 (fftw_complex *)aux2, 		  embed, padnpy,     1,
		// 		(fftw_complex *)aux2, 	 embed, padnpy,		1,		 FFTW_FORWARD,	FFTW_MEASURE   );
		// this->planxbac  = fftw_plan_many_dft(  1, &this->nx,	npy,	 (fftw_complex *)aux2, 		  embed, padnpy,     1,
		// 		(fftw_complex *)aux2, 	 embed, padnpy,		1,		 FFTW_BACKWARD,	FFTW_MEASURE   );
		// this->planyr2c  = fftw_plan_many_dft_r2c(  1, &this->bigny, this->nplane,	r_rspace , embed, this->nplane*2,      1,
		// 	(fftw_complex*)aux2, embed, this->nplane,		1,	FFTW_MEASURE   );
		// this->planyc2r  = fftw_plan_many_dft_c2r(  1, &this->bigny, this->nplane,	(fftw_complex*)aux2 , embed, this->nplane,      1,
		// 	r_rspace, embed, this->nplane*2,		1,			FFTW_MEASURE   );

	}
	else
	{
		// 	this->plan2for = fftw_plan_many_dft(       2,   nrank,  this->nplane,  
		// 						(fftw_complex*) aux2,   embed,  this->nplane,   1,
		// 						(fftw_complex*) aux2,  embed,  this->nplane,   1,  FFTW_FORWARD,  FFTW_MEASURE);

		// 	this->plan2bac = fftw_plan_many_dft(       2,   nrank,  this->nplane,  
		// 						(fftw_complex*) aux2,   embed,  this->nplane,   1,
		// 						(fftw_complex*) aux2,   embed,  this->nplane,   1,  FFTW_BACKWARD,  FFTW_MEASURE);

		//It is better to use in-place for stride > 1
		int bignpy = this->nplane * this->bigny;
		this->planxfor1  = fftw_plan_many_dft(  1, &this->nx,	this->nplane * (liy + 1),	 (fftw_complex *)aux2, 		  embed, bignpy,     1,
				(fftw_complex *)aux2, 	 embed, bignpy,		1,		 FFTW_FORWARD,	FFTW_MEASURE   );
		this->planxbac1  = fftw_plan_many_dft(  1, &this->nx,	this->nplane * (liy + 1),	 (fftw_complex *)aux2, 		  embed, bignpy,     1,
				(fftw_complex *)aux2, 	 embed, bignpy,		1,		 FFTW_BACKWARD,	FFTW_MEASURE   );
		this->planxfor2  = fftw_plan_many_dft(  1, &this->nx,	this->nplane * (bigny - riy),	 (fftw_complex *)aux2, 		  embed, bignpy,     1,
				(fftw_complex *)aux2, 	 embed, bignpy,		1,		 FFTW_FORWARD,	FFTW_MEASURE   );
		this->planxbac2  = fftw_plan_many_dft(  1, &this->nx,	this->nplane * (bigny - riy),	 (fftw_complex *)aux2, 		  embed, bignpy,     1,
				(fftw_complex *)aux2, 	 embed, bignpy,		1,		 FFTW_BACKWARD,	FFTW_MEASURE   );
		this->planyfor  = fftw_plan_many_dft(  1, &this->bigny, this->nplane,	(fftw_complex*)aux2 , embed, this->nplane,      1,
			(fftw_complex*)aux2, embed, this->nplane,		1,		 FFTW_FORWARD,	FFTW_MEASURE   );
		this->planybac  = fftw_plan_many_dft(  1, &this->bigny, this->nplane,	(fftw_complex*)aux2 , embed, this->nplane,      1,
			(fftw_complex*)aux2, embed, this->nplane,		1,		 FFTW_BACKWARD,	FFTW_MEASURE   );
	}

	

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
	
	//It is better to use out-of-place fft for stride = 1.
	this->planfzfor = fftwf_plan_many_dft(     1,    &this->nz,  this->ns,  
					    (fftwf_complex*) auxf1,  &this->nz,  1,  this->nz,
					    (fftwf_complex*) auxf2,  &this->nz,  1,  this->nz,  FFTW_FORWARD,  FFTW_MEASURE);
	
	this->planfzbac = fftwf_plan_many_dft(     1,    &this->nz,  this->ns,  
						(fftwf_complex*) auxf1,  &this->nz,  1,  this->nz,
						(fftwf_complex*) auxf2,  &this->nz,  1,  this->nz,  FFTW_BACKWARD,  FFTW_MEASURE);
	//---------------------------------------------------------
	//                              2 D
	//---------------------------------------------------------
	
	int *embed = NULL;
	int bignpy = this->nplane * this->bigny;
	this->planfxfor1  = fftwf_plan_many_dft(  1, &this->nx,	this->nplane * (liy + 1),	 (fftwf_complex *)auxf2, 		  embed, bignpy,     1,
			(fftwf_complex *)auxf2, 	 embed, bignpy,		1,		 FFTW_FORWARD,	FFTW_MEASURE   );
	this->planfxbac1  = fftwf_plan_many_dft(  1, &this->nx,	this->nplane * (liy + 1),	 (fftwf_complex *)auxf2, 		  embed, bignpy,     1,
			(fftwf_complex *)auxf2, 	 embed, bignpy,		1,		 FFTW_BACKWARD,	FFTW_MEASURE   );
	if(this->gamma_only)
	{
		this->planfyr2c  = fftwf_plan_many_dft_r2c(  1, &this->bigny, this->nplane,	rf_rspace , embed, this->nplane,      1,
			(fftwf_complex*)auxf1, embed, this->nplane,		1,	FFTW_MEASURE   );
		this->planfyc2r  = fftwf_plan_many_dft_c2r(  1, &this->bigny, this->nplane,	(fftwf_complex*)auxf1 , embed, this->nplane,      1,
			rf_rspace, embed, this->nplane,		1,			FFTW_MEASURE   );
	}
	else
	{
		this->planfxfor2  = fftwf_plan_many_dft(  1, &this->nx,	this->nplane * (bigny - riy),	 (fftwf_complex *)auxf2, 		  embed, bignpy,     1,
			(fftwf_complex *)auxf2, 	 embed, bignpy,		1,		 FFTW_FORWARD,	FFTW_MEASURE   );
		this->planfxbac2  = fftwf_plan_many_dft(  1, &this->nx,	this->nplane * (bigny - riy),	 (fftwf_complex *)auxf2, 		  embed, bignpy,     1,
			(fftwf_complex *)auxf2, 	 embed, bignpy,		1,		 FFTW_BACKWARD,	FFTW_MEASURE   );
		this->planfyfor  = fftwf_plan_many_dft(  1, &this->bigny, this->nplane,	(fftwf_complex*)auxf2 , embed, this->nplane,      1,
			(fftwf_complex*)auxf2, embed, this->nplane,		1,		 FFTW_FORWARD,	FFTW_MEASURE   );
		this->planfybac  = fftwf_plan_many_dft(  1, &this->bigny, this->nplane,	(fftwf_complex*)auxf2 , embed, this->nplane,      1,
			(fftwf_complex*)auxf2, embed, this->nplane,		1,		 FFTW_BACKWARD,	FFTW_MEASURE   );
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
	destroyp = true;
	return;
}
#ifdef __MIX_PRECISION
void FFT:: cleanfFFT()
{
	if(destroypf==true) return;
	fftwf_destroy_plan(planfzfor);
	fftwf_destroy_plan(planfzbac);
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
	destroypf = true;
	return;
}
#endif


void FFT::fftzfor(std::complex<double>* & in, std::complex<double>* & out)
{
	// for(int i = 0 ; i < this->ns ; ++i)
	// {
	// 	fftw_execute_dft(this->planzfor,(fftw_complex *)&in[i*nz],(fftw_complex *)&out[i*nz]);
	// }
	fftw_execute_dft(this->planzfor,(fftw_complex *)in,(fftw_complex *)out);
	return;
}

void FFT::fftzbac(std::complex<double>* & in, std::complex<double>* & out)
{
	// for(int i = 0 ; i < this->ns ; ++i)
	// {
	// 	fftw_execute_dft(this->planzbac,(fftw_complex *)&aux1[i*nz],(fftw_complex *)&aux1[i*nz]);
	// }
	fftw_execute_dft(this->planzbac,(fftw_complex *)in, (fftw_complex *)out);
	return;
}

void FFT::fftxyfor(std::complex<double>* & in, std::complex<double>* & out)
{
	int bignpy = this->nplane * this-> bigny;
	for (int i=0; i<this->nx;++i)
	{
		fftw_execute_dft( this->planyfor, (fftw_complex *)&in[i*bignpy], (fftw_complex *)&out[i*bignpy]);
	}


	fftw_execute_dft( this->planxfor1, (fftw_complex *)in, (fftw_complex *)out);
	fftw_execute_dft( this->planxfor2, (fftw_complex *)&in[riy*nplane], (fftw_complex *)&out[riy*nplane]);
	return;
}

void FFT::fftxybac(std::complex<double>* & in, std::complex<double>* & out)
{
	int bignpy = this->nplane * this-> bigny;
	//x-direction
	fftw_execute_dft( this->planxbac1, (fftw_complex *)in, (fftw_complex *)out);
	fftw_execute_dft( this->planxbac2, (fftw_complex *)&in[riy*nplane], (fftw_complex *)&out[riy*nplane]);

	////y-direction
	for (int i=0; i<this->nx;++i)
	{
		fftw_execute_dft( this->planybac, (fftw_complex*)&in[i*bignpy], (fftw_complex*)&out[i*bignpy] );
	}
	return;
}

void FFT::fftxyr2c(double* &in, std::complex<double>* & out)
{
	int bignpy = this->nplane * this-> bigny;

	for (int i=0; i<this->nx;++i)
	{
		fftw_execute_dft_r2c( this->planyr2c, &in[i*bignpy*2], (fftw_complex*)&out[i*bignpy] );
	}

	fftw_execute_dft( this->planxfor1, (fftw_complex *)out, (fftw_complex *)out);
	return;
}


void FFT::fftxyc2r(std::complex<double>* & in, double* & out)
{
	int bignpy = this->nplane * this-> bigny;
	fftw_execute_dft( this->planxbac1, (fftw_complex *)in, (fftw_complex *)in);
	
	for (int i=0; i<this->nx;++i)
	{
		fftw_execute_dft_c2r( this->planyc2r, (fftw_complex*)&in[i*bignpy], &out[i*bignpy*2] );
	}
	return;
}


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
	int bignpy = this->nplane * this-> bigny;
	for (int i=0; i<this->nx;++i)
	{
		fftwf_execute_dft( this->planfyfor, (fftwf_complex *)&in[i*bignpy], (fftwf_complex *)&out[i*bignpy]);
	}

	fftwf_execute_dft( this->planfxfor1, (fftwf_complex *)in, (fftwf_complex *)out);
	fftwf_execute_dft( this->planfxfor2, (fftwf_complex *)&in[riy*nplane], (fftwf_complex *)&out[riy*nplane]);
	return;
}

void FFT::fftfxybac(std::complex<float>* & in, std::complex<float>* & out)
{
	int bignpy = this->nplane * this-> bigny;
	//x-direction
	fftwf_execute_dft( this->planfxbac1, (fftwf_complex *)in, (fftwf_complex *)out);
	fftwf_execute_dft( this->planfxbac2, (fftwf_complex *)&in[riy*nplane], (fftwf_complex *)&out[riy*nplane]);

	////y-direction
	for (int i=0; i<this->nx;++i)
	{
		fftwf_execute_dft( this->planfybac, (fftwf_complex*)&in[i*bignpy], (fftwf_complex*)&out[i*bignpy] );
	}
	return;
}

void FFT::fftfxyr2c(float* &in, std::complex<float>* & out)
{
	int bignpy = this->nplane * this-> bigny;

	for (int i=0; i<this->nx;++i)
	{
		fftwf_execute_dft_r2c( this->planfyr2c, &in[i*bignpy*2], (fftwf_complex*)&out[i*bignpy] );
	}

	fftwf_execute_dft( this->planfxfor1, (fftwf_complex *)out, (fftwf_complex *)out);
	return;
}


void FFT::fftfxyc2r(std::complex<float>* & in, float* & out)
{
	int bignpy = this->nplane * this-> bigny;
	fftwf_execute_dft( this->planfxbac1, (fftwf_complex *)in, (fftwf_complex *)in);

	for (int i=0; i<this->nx;++i)
	{
		fftwf_execute_dft_c2r( this->planfyc2r, (fftwf_complex*)&in[i*bignpy], &out[i*bignpy*2] );
	}
	return;
}
#endif
}
