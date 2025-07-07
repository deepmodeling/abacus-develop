#ifndef KERNEL_H
#define KERNEL_H

#include <torch/torch.h>

class Kernel
{
    // ------------ fill the kernel in reciprocal space ----------
  public:
    Kernel(){};

    void set_para(const int kernel_type_in, const double scaling_in, const double yukawa_alpha_in, const std::string &kernel_file_in)
    {
        this->kernel_type = kernel_type_in;
        this->scaling = scaling_in;
        this->yukawa_alpha = yukawa_alpha_in;
        this->kernel_file = kernel_file_in;
    }

    int kernel_type = 0; // 1: WT, 2: Yukawa
    double scaling = 0.;
    double yukawa_alpha = 0.;
    std::string kernel_file = "none";
    std::vector<torch::Tensor> kernel;

    void fill_kernel(const int fftdim,
                     const int ndata,
                     const torch::Tensor &rho,
                     const double *volume,
                     const std::string *cell,
                     const torch::Device device,
                     const std::vector<torch::Tensor> &fft_gg);
    double wt_kernel(double eta, double tf_weight = 1., double vw_weight = 1.);
    double yukawa_kernel(double eta, double alpha);
    void read_kernel(const int fftdim,
                     const int ndata,
                     const torch::Tensor &rho,
                     const double *volume,
                     const std::string *cell,
                     const torch::Device device,
                     const std::vector<torch::Tensor> &fft_gg);
};
#endif