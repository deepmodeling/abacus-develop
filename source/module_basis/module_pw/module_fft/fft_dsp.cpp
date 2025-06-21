#include "fft_dsp.h"

#include "source_base/global_variable.h"

#include <iostream>
#include <string.h>
#include <vector>
namespace ModulePW
{
template <>
void FFT_DSP<double>::initfft(int nx_in, int ny_in, int nz_in)
{
    this->nx = nx_in;
    this->ny = ny_in;
    this->nz = nz_in;
    cluster_id = GlobalV::MY_RANK;
    nxyz = this->nx * this->ny * this->nz;
}
template <>
void FFT_DSP<double>::setupFFT()
{
    PROBLEM pbm_forward;
    PROBLEM pbm_backward;
    PLAN* ptr_plan_forward;
    PLAN* ptr_plan_backward;
    INT num_thread = 8;
    INT size=0;
    hthread_dat_load(cluster_id, FFT_DAT_DIR);

    // compute the size of and malloc thread
    size = nx * ny * nz * 2 * sizeof(E);
    forward_in = (E*)hthread_malloc((int)cluster_id, size, HT_MEM_RW);

    //init 3d fft problem
    pbm_forward.num_dim = 3;
    pbm_forward.n[0] = nx;
    pbm_forward.n[1] = ny;
    pbm_forward.n[2] = nz;
    pbm_forward.iFFT = 0;
    pbm_forward.in = forward_in;
    pbm_forward.out = forward_in;

    //make ptr plan
    make_plan(&pbm_forward, &ptr_plan_forward, cluster_id, num_thread);
    ptr_plan_forward->in = forward_in;
    ptr_plan_forward->out = forward_in;
    args_for[1] = (unsigned long)ptr_plan_forward;

    // init 3d fft problem
    pbm_backward.num_dim = 3;
    pbm_backward.n[0] = nx;
    pbm_backward.n[1] = ny;
    pbm_backward.n[2] = nz;
    pbm_backward.iFFT = 1;
    pbm_backward.in = forward_in;
    pbm_backward.out = forward_in;

    make_plan(&pbm_backward, &ptr_plan_backward, cluster_id, num_thread);
    ptr_plan_backward->in = forward_in;
    ptr_plan_backward->out = forward_in;
    args_back[1] = (unsigned long)ptr_plan_backward;
}
template <>
void FFT_DSP<double>::resource_handler(const int flag) const
{
    if (flag == 0)
    {
        hthread_barrier_destroy(b_id);
        hthread_group_destroy(thread_id_for);
    }
    else if (flag==1)
    {
        INT num_thread = 8;
        thread_id_for = hthread_group_create(cluster_id, num_thread, NULL, 0, 0, NULL);
        // create b_id for the barrier
        b_id = hthread_barrier_create(cluster_id);
        args_for[0] = b_id;
        args_back[0] = b_id;
    }else{
        ModuleBase::WARNING_QUIT("FFT_DSP", "Error use of fft resource handle");
    }
}
template <>
void FFT_DSP<double>::fft3D_forward(std::complex<double>* in, std::complex<double>* out) const
{
    hthread_group_exec(thread_id_for, "execute_device", 1, 1, args_for);
    hthread_group_wait(thread_id_for);
}

template <>
void FFT_DSP<double>::fft3D_backward(std::complex<double>* in, std::complex<double>* out) const
{
    hthread_group_exec(thread_id_for, "execute_device", 1, 1, args_back);
    hthread_group_wait(thread_id_for);
}
template <>
void FFT_DSP<double>::cleanFFT()
{
    if (ptr_plan_forward != nullptr)
    {
        destroy_plan(ptr_plan_forward);
        ptr_plan_forward = nullptr;
    }
    if (ptr_plan_backward != nullptr)
    {
        destroy_plan(ptr_plan_backward);
        ptr_plan_backward = nullptr;
    }
}

template <>
void FFT_DSP<double>::clear()
{
    this->cleanFFT();
    hthread_free(forward_in);
}

template <>
std::complex<double>* FFT_DSP<double>::get_auxr_3d_data() const
{
    return reinterpret_cast<std::complex<double>*>(this->forward_in);
}
template FFT_DSP<float>::FFT_DSP();
template FFT_DSP<float>::~FFT_DSP();
template FFT_DSP<double>::FFT_DSP();
template FFT_DSP<double>::~FFT_DSP();
} // namespace ModulePW