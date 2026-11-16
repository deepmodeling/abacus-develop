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
    // Batch the independent x slabs and z planes in one y-transform plan.
    // The array layout is [x][y][z], with z distributed across MPI ranks.
    auto plan_y = [&](int first_x, int count_x, int sign) {
        fftwf_iodim transform = {this->ny, this->nplane, this->nplane};
        fftwf_iodim batches[2] = {{this->nplane, 1, 1}, {count_x, npy, npy}};
        auto* data = reinterpret_cast<fftwf_complex*>(c_auxr + first_x * npy);
        return fftwf_plan_guru_dft(1, &transform, 2, batches, data, data, sign, flag);
    };
    if (this->xprime)
    {
        this->planfyfor = plan_y(0, this->lixy + 1, FFTW_FORWARD);
        this->planfybac = plan_y(0, this->lixy + 1, FFTW_BACKWARD);
        if (!this->gamma_only && this->rixy < this->nx)
        {
            this->planfyfor2 = plan_y(this->rixy, this->nx - this->rixy, FFTW_FORWARD);
            this->planfybac2 = plan_y(this->rixy, this->nx - this->rixy, FFTW_BACKWARD);
        }
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
            this->planfyfor = plan_y(0, this->nx, FFTW_FORWARD);
            this->planfybac = plan_y(0, this->nx, FFTW_BACKWARD);
        }
    }
    return;
}

template <>
void FFT_CPU<float>::clearfft(fftwf_plan& plan)
{
    if (plan)
    {
        fftwf_destroy_plan(plan);
        plan = nullptr;
    }
}

template <>
void FFT_CPU<float>::cleanFFT()
{
    clearfft(planfzfor);
    clearfft(planfzbac);
    clearfft(planfxfor1);
    clearfft(planfxbac1);
    clearfft(planfxfor2);
    clearfft(planfxbac2);
    clearfft(planfyfor);
    clearfft(planfybac);
    clearfft(planfyfor2);
    clearfft(planfybac2);
    clearfft(planfxr2c);
    clearfft(planfxc2r);
    clearfft(planfyr2c);
    clearfft(planfyc2r);
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

        fftwf_execute_dft(this->planfyfor, (fftwf_complex*)in, (fftwf_complex*)out);
        if (this->planfyfor2)
        {
            fftwf_execute_dft(this->planfyfor2, (fftwf_complex*)&in[rixy * npy], (fftwf_complex*)&out[rixy * npy]);
        }
    }
    else
    {
        fftwf_execute_dft(this->planfyfor, (fftwf_complex*)in, (fftwf_complex*)out);

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
        fftwf_execute_dft(this->planfybac, (fftwf_complex*)in, (fftwf_complex*)out);
        if (this->planfybac2)
        {
            fftwf_execute_dft(this->planfybac2, (fftwf_complex*)&in[rixy * npy], (fftwf_complex*)&out[rixy * npy]);
        }

        fftwf_execute_dft(this->planfxbac1, (fftwf_complex*)in, (fftwf_complex*)out);
    }
    else
    {
        fftwf_execute_dft(this->planfxbac1, (fftwf_complex*)in, (fftwf_complex*)out);
        fftwf_execute_dft(this->planfxbac2, (fftwf_complex*)&in[rixy * nplane], (fftwf_complex*)&out[rixy * nplane]);

        fftwf_execute_dft(this->planfybac, (fftwf_complex*)in, (fftwf_complex*)out);
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

        fftwf_execute_dft(this->planfyfor, (fftwf_complex*)out, (fftwf_complex*)out);
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
        fftwf_execute_dft(this->planfybac, (fftwf_complex*)in, (fftwf_complex*)in);

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