#include "fft_cpu.h"
#include "fftw3.h"
namespace ModulePW
{

template <typename FPTYPE>
void FFT_CPU<FPTYPE>::initfft(int nx_in, 
                              int ny_in, 
                              int nz_in, 
                              int lixy_in, 
                              int rixy_in, 
                              int ns_in, 
                              int nplane_in, 
				              int nproc_in, 
                              bool gamma_only_in, 
                              bool xprime_in)
{
    this->gamma_only = gamma_only_in;
    this->xprime = xprime_in;
    this->fftnx = this->nx = nx_in;
    this->fftny = this->ny = ny_in;
    if (this->gamma_only)
    {
        if (xprime) {
            this->fftnx = int(this->nx / 2) + 1;
        } else {
            this->fftny = int(this->ny / 2) + 1;
        }
    }
    this->nz = nz_in;
    this->ns = ns_in;
    this->lixy = lixy_in;
    this->rixy = rixy_in;
    this->nplane = nplane_in;
    this->nproc = nproc_in;
    this->nxy = this->nx * this->ny;
    this->fftnxy = this->fftnx * this->fftny;
    const int nrxx = this->nxy * this->nplane;
    const int nsz = this->nz * this->ns;
    this->maxgrids = (nsz > nrxx) ? nsz : nrxx;
}
template <>
void FFT_CPU<double>::setupFFT()
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
    z_auxg = (std::complex<double>*)fftw_malloc(sizeof(fftw_complex) * this->maxgrids);
    z_auxr = (std::complex<double>*)fftw_malloc(sizeof(fftw_complex) * this->maxgrids);
    d_rspace = (double*)z_auxg;
    this->planzfor = fftw_plan_many_dft(1, 
                                        &this->nz, 
                                        this->ns, 
                                        (fftw_complex*)z_auxg, 
                                        &this->nz, 
                                        1, 
                                        this->nz,
                                        (fftw_complex*)z_auxg, 
                                        &this->nz, 
                                        1, 
                                        this->nz, 
                                        FFTW_FORWARD, 
                                        flag);

    this->planzbac = fftw_plan_many_dft(1, 
                                        &this->nz, 
                                        this->ns, 
                                        (fftw_complex*)z_auxg, 
                                        &this->nz, 
                                        1, 
                                        this->nz,
                                        (fftw_complex*)z_auxg, 
                                        &this->nz, 
                                        1, 
                                        this->nz, 
                                        FFTW_BACKWARD, 
                                        flag);

    //---------------------------------------------------------
    //                              2 D - XY
    //---------------------------------------------------------
    // 1D+1D is much faster than 2D FFT!
    // in-place fft is better for c2c and out-of-place fft is better for c2r
    int* embed = nullptr;
    int npy = this->nplane * this->ny;
    if (this->xprime)
    {
        this->planyfor = fftw_plan_many_dft(1, 
                                            &this->ny, 
                                            this->nplane, 
                                            (fftw_complex*)z_auxr, 
                                            embed,
                                            this->nplane, 
                                            1,
                                            (fftw_complex*)z_auxr, 
                                            embed,
                                            this->nplane, 
                                            1, 
                                            FFTW_FORWARD, 
                                            flag);
        this->planybac = fftw_plan_many_dft(1, 
                                            &this->ny, 
                                            this->nplane, 
                                            (fftw_complex*)z_auxr, 
                                            embed,
                                            this->nplane, 
                                            1,
                                            (fftw_complex*)z_auxr, 
                                            embed,
                                            this->nplane, 
                                            1, 
                                            FFTW_BACKWARD, 
                                            flag);
        if (this->gamma_only)
        {
            this->planxr2c = fftw_plan_many_dft_r2c(1,  
                                                    &this->nx, 
                                                    npy, 
                                                    d_rspace, 
                                                    embed, 
                                                    npy, 
                                                    1, 
                                                    (fftw_complex*)z_auxr,
                                                    embed, 
                                                    npy, 
                                                    1, 
                                                    flag);
            this->planxc2r = fftw_plan_many_dft_c2r(1, 
                                                    &this->nx, 
                                                    npy, 
                                                    (fftw_complex*)z_auxr, 
                                                    embed, 
                                                    npy, 
                                                    1, 
                                                    d_rspace,
                                                    embed, 
                                                    npy, 
                                                    1, 
                                                    flag);
        }
        else
        {
            this->planxfor1 = fftw_plan_many_dft(1, 
                                                 &this->nx, 
                                                 npy, 
                                                 (fftw_complex*)z_auxr, 
                                                 embed, 
                                                 npy, 
                                                 1,
                                                 (fftw_complex*)z_auxr, 
                                                 embed, 
                                                 npy, 
                                                 1, 
                                                 FFTW_FORWARD, 
                                                 flag);
            this->planxbac1 = fftw_plan_many_dft(1, 
                                                 &this->nx, 
                                                 npy, 
                                                 (fftw_complex*)z_auxr, 
                                                 embed, 
                                                 npy, 
                                                 1,
                                                (fftw_complex*)z_auxr, 
                                                embed, 
                                                npy, 
                                                1, 
                                                FFTW_BACKWARD, 
                                                flag);
        }
    }
    else
    {
        this->planxfor1 = fftw_plan_many_dft(1, 
                                             &this->nx, 
                                             this->nplane * (this->lixy + 1), 
                                             (fftw_complex*)z_auxr, 
                                             embed, 
                                             npy,
                                             1, 
                                             (fftw_complex*)z_auxr, 
                                             embed, 
                                             npy, 
                                             1, 
                                             FFTW_FORWARD, 
                                             flag);
        this->planxbac1 = fftw_plan_many_dft(1, 
                                            &this->nx, 
                                            this->nplane * (this->lixy + 1), 
                                            (fftw_complex*)z_auxr, 
                                            embed, 
                                            npy,
                                            1, 
                                            (fftw_complex*)z_auxr, 
                                            embed, 
                                            npy, 
                                            1, 
                                            FFTW_BACKWARD, 
                                            flag);
        if (this->gamma_only)
        {
            this->planyr2c = fftw_plan_many_dft_r2c(1, 
                                                    &this->ny, 
                                                    this->nplane,
                                                    d_rspace, 
                                                    embed, 
                                                    this->nplane, 
                                                    1,
                                                    (fftw_complex*)z_auxr, 
                                                    embed, 
                                                    this->nplane, 
                                                    1, 
                                                    flag);
            this->planyc2r = fftw_plan_many_dft_c2r(1, 
                                                    &this->ny, 
                                                    this->nplane, 
                                                    (fftw_complex*)z_auxr, 
                                                    embed,
                                                    this->nplane, 
                                                    1, 
                                                    d_rspace, 
                                                    embed, 
                                                    this->nplane, 
                                                    1, 
                                                    flag);
        }
        else
        {
            this->planxfor2 = fftw_plan_many_dft(1, 
                                                &this->nx, 
                                                this->nplane * (this->ny - this->rixy), 
                                                (fftw_complex*)z_auxr, 
                                                embed,
                                                npy, 
                                                1, (fftw_complex*)z_auxr, 
                                                embed, 
                                                npy, 
                                                1, 
                                                FFTW_FORWARD, 
                                                flag);
            this->planxbac2 = fftw_plan_many_dft(1, 
                                                 &this->nx, 
                                                 this->nplane * (this->ny - this->rixy), 
                                                 (fftw_complex*)z_auxr, 
                                                 embed,
                                                 npy, 
                                                 1, 
                                                 (fftw_complex*)z_auxr, 
                                                 embed, 
                                                 npy, 
                                                 1, 
                                                 FFTW_BACKWARD, 
                                                 flag);
            this->planyfor = fftw_plan_many_dft(1, 
                                                &this->ny, 
                                                this->nplane, 
                                                (fftw_complex*)z_auxr, 
                                                embed, 
                                                this->nplane,
                                                1, 
                                                (fftw_complex*)z_auxr, 
                                                embed, 
                                                this->nplane, 
                                                1, 
                                                FFTW_FORWARD, 
                                                flag);
            this->planybac = fftw_plan_many_dft(1, 
                                                &this->ny, 
                                                this->nplane, 
                                                (fftw_complex*)z_auxr, 
                                                embed, 
                                                this->nplane,
                                                1, 
                                                (fftw_complex*)z_auxr, 
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
void FFT_CPU<double>::clearfft(fftw_plan& plan)
{
    if (plan)
    {
        fftw_destroy_plan(plan);
        plan = nullptr;
    }
}

template <>
void FFT_CPU<double>::cleanFFT()
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
void FFT_CPU<double>::clear()
{
    this->cleanFFT();
    if (z_auxg != nullptr)
    {
        fftw_free(z_auxg);
        z_auxg = nullptr;
    }
    if (z_auxr != nullptr)
    {
        fftw_free(z_auxr);
        z_auxr = nullptr;
    }
    d_rspace = nullptr;
}

template <>
void FFT_CPU<double>::fftxyfor(std::complex<double>* in, std::complex<double>* out) const
{
    int npy = this->nplane * this->ny;
    if (this->xprime)
    {
        fftw_execute_dft(this->planxfor1, (fftw_complex*)in, (fftw_complex*)out);
        for (int i = 0; i < this->lixy + 1; ++i)
        {
            fftw_execute_dft(this->planyfor, (fftw_complex*)&in[i * npy], (fftw_complex*)&out[i * npy]);
        }
        for (int i = rixy; i < this->nx; ++i)
        {
            fftw_execute_dft(this->planyfor, (fftw_complex*)&in[i * npy], (fftw_complex*)&out[i * npy]);
        }
    }
    else
    {
        for (int i = 0; i < this->nx; ++i)
        {
            fftw_execute_dft(this->planyfor, (fftw_complex*)&in[i * npy], (fftw_complex*)&out[i * npy]);
        }
        fftw_execute_dft(this->planxfor1, (fftw_complex*)in, (fftw_complex*)out);
        fftw_execute_dft(this->planxfor2, (fftw_complex*)&in[rixy * nplane], (fftw_complex*)&out[rixy * nplane]);
    }
}

template <>
void FFT_CPU<double>::fftxybac(std::complex<double>* in,std::complex<double>* out) const
{
    int npy = this->nplane * this->ny;
    if (this->xprime)
    {
        for (int i = 0; i < this->lixy + 1; ++i)
        {
            fftw_execute_dft(this->planybac, (fftw_complex*)&in[i * npy], (fftw_complex*)&out[i * npy]);
        }
        for (int i = rixy; i < this->nx; ++i)
        {
            fftw_execute_dft(this->planybac, (fftw_complex*)&in[i * npy], (fftw_complex*)&out[i * npy]);
        }
        fftw_execute_dft(this->planxbac1, (fftw_complex*)in, (fftw_complex*)out);
    }
    else
    {
        fftw_execute_dft(this->planxbac1, (fftw_complex*)in, (fftw_complex*)out);
        fftw_execute_dft(this->planxbac2, (fftw_complex*)&in[rixy * nplane], (fftw_complex*)&out[rixy * nplane]);
        for (int i = 0; i < this->nx; ++i)
        {
            fftw_execute_dft(this->planybac, (fftw_complex*)&in[i * npy], (fftw_complex*)&out[i * npy]);
        }
    }
}

template <>
void FFT_CPU<double>::fftzfor(std::complex<double>* in, std::complex<double>* out) const
{
    fftw_execute_dft(this->planzfor, (fftw_complex*)in, (fftw_complex*)out);
}

template <>
void FFT_CPU<double>::fftzbac(std::complex<double>* in, std::complex<double>* out) const
{
    fftw_execute_dft(this->planzbac, (fftw_complex*)in, (fftw_complex*)out);
}

template <>
void FFT_CPU<double>::fftxyr2c(double* in, std::complex<double>* out) const
{
    int npy = this->nplane * this->ny;
    if (this->xprime)
    {
        fftw_execute_dft_r2c(this->planxr2c, in, (fftw_complex*)out);
        for (int i = 0; i < this->lixy + 1; ++i)
        {
            fftw_execute_dft(this->planyfor, (fftw_complex*)&out[i * npy], (fftw_complex*)&out[i * npy]);
        }
    }
    else
    {
        for (int i = 0; i < this->nx; ++i)
        {
            fftw_execute_dft_r2c(this->planyr2c, &in[i * npy], (fftw_complex*)&out[i * npy]);
        }
        fftw_execute_dft(this->planxfor1, (fftw_complex*)out, (fftw_complex*)out);
    }
}

template <>
void FFT_CPU<double>::fftxyc2r(std::complex<double> *in,double *out) const
{
    int npy = this->nplane * this->ny;
    if (this->xprime)
    {
        for (int i = 0; i < this->lixy + 1; ++i)
        {
            fftw_execute_dft(this->planybac, (fftw_complex*)&in[i * npy], (fftw_complex*)&in[i * npy]);
        }
        fftw_execute_dft_c2r(this->planxc2r, (fftw_complex*)in, out);
    }
    else
    {
        fftw_execute_dft(this->planxbac1, (fftw_complex*)in, (fftw_complex*)in);
        for (int i = 0; i < this->nx; ++i)
        {
            fftw_execute_dft_c2r(this->planyc2r, (fftw_complex*)&in[i * npy], &out[i * npy]);
        }
    }
}

template <> double* 
FFT_CPU<double>::get_rspace_data() const {return d_rspace;}
template <> std::complex<double>* 
FFT_CPU<double>::get_auxr_data()   const {return z_auxr;}
template <> std::complex<double>* 
FFT_CPU<double>::get_auxg_data()   const {return z_auxg;}

template FFT_CPU<float>::FFT_CPU();
template FFT_CPU<float>::~FFT_CPU();
template FFT_CPU<double>::FFT_CPU();
template FFT_CPU<double>::~FFT_CPU();
}