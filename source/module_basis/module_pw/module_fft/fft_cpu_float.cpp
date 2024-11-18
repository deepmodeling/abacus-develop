#include "fft_cpu.h"

namespace ModulePW
{
template <>
void FFT_CPU<float>::setupFFT()
{
    unsigned int flag = FFTW_ESTIMATE;
    switch (this->fft_mode)
    {
    case 0:
        flag = FFTW_ESTIMATE;
        break;
    case 1:
        flag = FFTW_MEASURE;
        break;
    case 2:
        flag = FFTW_PATIENT;
        break;
    case 3:
        flag = FFTW_EXHAUSTIVE;
        break;
    default:
        break;
    }
    c_auxg = (std::complex<float>*)fftwf_malloc(sizeof(fftwf_complex) * this->maxgrids); 
    c_auxr = (std::complex<float>*)fftwf_malloc(sizeof(fftwf_complex) * this->maxgrids);
    s_rspace = (float*)c_auxg;
    //---------------------------------------------------------
    //                              1 D
    //---------------------------------------------------------

    //               fftw_plan_many_dft(int rank,          
    //                                  const int *n,       int howmany,
    //					                fftw_complex *in,  const int *inembed, int istride, int idist,
    //					                fftw_complex *out, const int *onembed, int ostride, int odist, int sign, unsigned
    //flags);

    this->planfzfor = fftwf_plan_many_dft(1, 
                                          &this->nz, 
                                          this->ns, 
                                          (fftwf_complex*)c_auxg, 
                                          &this->nz, 
                                          1, 
                                          this->nz,
                                          (fftwf_complex*)c_auxg, 
                                          &this->nz, 
                                          1, 
                                          this->nz, 
                                          FFTW_FORWARD, 
                                          flag);

    this->planfzbac = fftwf_plan_many_dft(1, 
                                          &this->nz, 
                                          this->ns, 
                                          (fftwf_complex*)c_auxg, 
                                          &this->nz, 
                                          1, 
                                          this->nz,
                                          (fftwf_complex*)c_auxg, 
                                          &this->nz, 
                                          1, 
                                          this->nz, 
                                          FFTW_BACKWARD, 
                                          flag);
    //---------------------------------------------------------
    //                              2 D
    //---------------------------------------------------------

    int* embed = nullptr;
    int npy = this->nplane * this->ny;
    if (this->xprime)
    {
        this->planfyfor = fftwf_plan_many_dft(1, 
                                              &this->ny, 
                                              this->nplane, 
                                              (fftwf_complex*)c_auxr, 
                                              embed, 
                                              nplane, 
                                              1,
                                              (fftwf_complex*)c_auxr, 
                                              embed, 
                                              nplane, 
                                              1, 
                                              FFTW_FORWARD, 
                                              flag);
        this->planfybac = fftwf_plan_many_dft(1, 
                                              &this->ny, 
                                              this->nplane, 
                                              (fftwf_complex*)c_auxr, 
                                              embed, 
                                              nplane, 
                                              1,
                                              (fftwf_complex*)c_auxr, 
                                              embed, nplane, 
                                              1, 
                                              FFTW_BACKWARD, 
                                              flag);
        if (this->gamma_only)
        {
            this->planfxr2c = fftwf_plan_many_dft_r2c(1, 
                                                      &this->nx, 
                                                      npy, 
                                                      s_rspace, 
                                                      embed, 
                                                      npy, 
                                                      1,
                                                      (fftwf_complex*)c_auxr, 
                                                      embed, npy, 
                                                      1, 
                                                      flag);
            this->planfxc2r = fftwf_plan_many_dft_c2r(1, 
                                                      &this->nx, 
                                                      npy, 
                                                      (fftwf_complex*)c_auxr, 
                                                      embed, 
                                                      npy,  
                                                      1,
                                                      s_rspace, 
                                                      embed, 
                                                      npy, 
                                                      1, 
                                                      flag);
        }
        else
        {
            this->planfxfor1 = fftwf_plan_many_dft(1, 
                                                   &this->nx, 
                                                   npy, 
                                                   (fftwf_complex*)c_auxr, 
                                                   embed, 
                                                   npy, 
                                                   1,
                                                   (fftwf_complex*)c_auxr, 
                                                   embed, 
                                                   npy, 
                                                   1, 
                                                   FFTW_FORWARD, 
                                                   flag);
            this->planfxbac1 = fftwf_plan_many_dft(1, 
                                                   &this->nx, 
                                                   npy, 
                                                   (fftwf_complex*)c_auxr, 
                                                   embed, 
                                                   npy, 
                                                   1,
                                                   (fftwf_complex*)c_auxr, 
                                                   embed, 
                                                   npy, 
                                                   1, 
                                                   FFTW_BACKWARD, 
                                                   flag);
        }
    }
    else
    {
        this->planfxfor1 = fftwf_plan_many_dft(1, 
                                               &this->nx, 
                                               this->nplane * (lixy + 1), 
                                               (fftwf_complex*)c_auxr, 
                                               embed,
                                               npy, 
                                               1, 
                                               (fftwf_complex*)c_auxr, 
                                               embed, 
                                               npy, 
                                               1, 
                                               FFTW_FORWARD, 
                                               flag);
        this->planfxbac1 = fftwf_plan_many_dft(1, 
                                               &this->nx, 
                                               this->nplane * (lixy + 1), 
                                               (fftwf_complex*)c_auxr, 
                                               embed,
                                               npy, 
                                               1, 
                                               (fftwf_complex*)c_auxr, 
                                               embed, 
                                               npy, 
                                               1, 
                                               FFTW_BACKWARD, 
                                               flag);
        if (this->gamma_only)
        {
            this->planfyr2c = fftwf_plan_many_dft_r2c(1, 
                                                      &this->ny, 
                                                      this->nplane, 
                                                      s_rspace, 
                                                      embed, 
                                                      this->nplane, 
                                                      1,
                                                      (fftwf_complex*)c_auxr, 
                                                      embed, 
                                                      this->nplane, 
                                                      1, 
                                                      flag);
            this->planfyc2r = fftwf_plan_many_dft_c2r(1, 
                                                      &this->ny, 
                                                      this->nplane, 
                                                      (fftwf_complex*)c_auxr, 
                                                      embed,
                                                      this->nplane, 
                                                      1, 
                                                      s_rspace, 
                                                      embed, 
                                                      this->nplane, 
                                                      1, 
                                                      flag);
        }
        else
        {
            this->planfxfor2 = fftwf_plan_many_dft(1, 
                                                  &this->nx, 
                                                  this->nplane * (this->ny - rixy), 
                                                  (fftwf_complex*)c_auxr, 
                                                  embed,
                                                  npy, 
                                                  1, 
                                                  (fftwf_complex*)c_auxr, 
                                                  embed, 
                                                  npy, 
                                                  1, 
                                                  FFTW_FORWARD, 
                                                  flag);
            this->planfxbac2 = fftwf_plan_many_dft(1, 
                                                   &this->nx, 
                                                   this->nplane * (this->ny - rixy), 
                                                   (fftwf_complex*)c_auxr, 
                                                   embed,
                                                   npy, 
                                                   1, 
                                                   (fftwf_complex*)c_auxr, 
                                                   embed, 
                                                   npy, 
                                                   1, 
                                                   FFTW_BACKWARD, 
                                                   flag);
            this->planfyfor = fftwf_plan_many_dft(1, 
                                                  &this->ny, 
                                                  this->nplane, 
                                                  (fftwf_complex*)c_auxr, 
                                                  embed, 
                                                  this->nplane, 
                                                  1,
                                                  (fftwf_complex*)c_auxr, 
                                                  embed, 
                                                  this->nplane, 
                                                  1, 
                                                  FFTW_FORWARD, 
                                                  flag);
            this->planfybac = fftwf_plan_many_dft(1, 
                                                  &this->ny, 
                                                  this->nplane, 
                                                  (fftwf_complex*)c_auxr, 
                                                  embed, 
                                                  this->nplane, 
                                                  1,
                                                  (fftwf_complex*)c_auxr, 
                                                  embed, 
                                                  this->nplane, 
                                                  1, 
                                                  FFTW_BACKWARD, 
                                                  flag);
        }
    }
    return;
}

template <>
void FFT_CPU<float>::clearfft(fftw_plan& plan)
{
    if (plan)
    {
        fftw_destroy_plan(plan);
        plan = nullptr;
    }
}

template <>
void FFT_CPU<float>::cleanFFT()
{
    clearfft(planzfor);
    clearfft(planzbac);
    clearfft(planxfor1);
    clearfft(planxbac1);
    clearfft(planxfor2);
    clearfft(planxbac2);
    clearfft(planyfor);
    clearfft(planybac);
    clearfft(planxr2c);
    clearfft(planxc2r);
    clearfft(planyr2c);
    clearfft(planyc2r);
}


template <>
void FFT_CPU<float>::clear()
{
    this->cleanFFT();
    if (c_auxg != nullptr)
    {
        fftw_free(c_auxg);
        c_auxg = nullptr;
    }
    if (c_auxr != nullptr)
    {
        fftw_free(c_auxr);
        c_auxr = nullptr;
    }
    s_rspace = nullptr;
}


template <>
void FFT_CPU<float>::fftxyfor(std::complex<float>* in, std::complex<float>* out) const
{
    int npy = this->nplane * this->ny;
    if (this->xprime)
    {
        fftwf_execute_dft(this->planfxfor1, (fftwf_complex*)in, (fftwf_complex*)out);

        for (int i = 0; i < this->lixy + 1; ++i)
        {
            fftwf_execute_dft(this->planfyfor, (fftwf_complex*)&in[i * npy], (fftwf_complex*)&out[i * npy]);
        }
        for (int i = rixy; i < this->nx; ++i)
        {
            fftwf_execute_dft(this->planfyfor, (fftwf_complex*)&in[i * npy], (fftwf_complex*)&out[i * npy]);
        }
    }
    else
    {
        for (int i = 0; i < this->nx; ++i)
        {
            fftwf_execute_dft(this->planfyfor, (fftwf_complex*)&in[i * npy], (fftwf_complex*)&out[i * npy]);
        }

        fftwf_execute_dft(this->planfxfor1, (fftwf_complex*)in, (fftwf_complex*)out);
        fftwf_execute_dft(this->planfxfor2, (fftwf_complex*)&in[rixy * nplane], (fftwf_complex*)&out[rixy * nplane]);
    }
}
template <>
void FFT_CPU<float>::fftxybac(std::complex<float>* in,std::complex<float> * out) const
{
    int npy = this->nplane * this->ny;
    if (this->xprime)
    {
        for (int i = 0; i < this->lixy + 1; ++i)
        {
            fftwf_execute_dft(this->planfybac, (fftwf_complex*)&in[i * npy], (fftwf_complex*)&out[i * npy]);
        }
        for (int i = rixy; i < this->nx; ++i)
        {
            fftwf_execute_dft(this->planfybac, (fftwf_complex*)&in[i * npy], (fftwf_complex*)&out[i * npy]);
        }

        fftwf_execute_dft(this->planfxbac1, (fftwf_complex*)in, (fftwf_complex*)out);
    }
    else
    {
        fftwf_execute_dft(this->planfxbac1, (fftwf_complex*)in, (fftwf_complex*)out);
        fftwf_execute_dft(this->planfxbac2, (fftwf_complex*)&in[rixy * nplane], (fftwf_complex*)&out[rixy * nplane]);

        for (int i = 0; i < this->nx; ++i)
        {
            fftwf_execute_dft(this->planfybac, (fftwf_complex*)&in[i * npy], (fftwf_complex*)&out[i * npy]);
        }
    }
}
template <>
void FFT_CPU<float>::fftzfor(std::complex<float>* in, std::complex<float>* out) const
{
    fftwf_execute_dft(this->planfzfor, (fftwf_complex*)in, (fftwf_complex*)out);
}
template <>
void FFT_CPU<float>::fftzbac(std::complex<float>* in, std::complex<float>* out) const
{
    fftwf_execute_dft(this->planfzbac, (fftwf_complex*)in, (fftwf_complex*)out);
}
template <>
void FFT_CPU<float>::fftxyr2c(float* in, std::complex<float>* out) const
{
    int npy = this->nplane * this->ny;
    if (this->xprime)
    {
        fftwf_execute_dft_r2c(this->planfxr2c, in, (fftwf_complex*)out);

        for (int i = 0; i < this->lixy + 1; ++i)
        {
            fftwf_execute_dft(this->planfyfor, (fftwf_complex*)&out[i * npy], (fftwf_complex*)&out[i * npy]);
        }
    }
    else
    {
        for (int i = 0; i < this->nx; ++i)
        {
            fftwf_execute_dft_r2c(this->planfyr2c, &in[i * npy], (fftwf_complex*)&out[i * npy]);
        }

        fftwf_execute_dft(this->planfxfor1, (fftwf_complex*)out, (fftwf_complex*)out);
    }
}
template <>
void FFT_CPU<float>::fftxyc2r(std::complex<float>* in, float* out) const
{
    int npy = this->nplane * this->ny;
    if (this->xprime)
    {
        for (int i = 0; i < this->lixy + 1; ++i)
        {
            fftwf_execute_dft(this->planfybac, (fftwf_complex*)&in[i * npy], (fftwf_complex*)&in[i * npy]);
        }

        fftwf_execute_dft_c2r(this->planfxc2r, (fftwf_complex*)in, out);
    }
    else
    {
        fftwf_execute_dft(this->planfxbac1, (fftwf_complex*)in, (fftwf_complex*)in);

        for (int i = 0; i < this->nx; ++i)
        {
            fftwf_execute_dft_c2r(this->planfyc2r, (fftwf_complex*)&in[i * npy], &out[i * npy]);
        }
    }
}
template <> float* 
FFT_CPU<float>::get_rspace_data() const {return s_rspace;}
template <> std::complex<float>* 
FFT_CPU<float>::get_auxr_data()   const {return c_auxr;}
template <> std::complex<float>* 
FFT_CPU<float>::get_auxg_data()   const {return c_auxg;}
}