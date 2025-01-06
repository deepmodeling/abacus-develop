#ifndef DATA_H
#define DATA_H

#include "./input.h"

#include <torch/torch.h>

class Data
{
  public:
    // --------- load the data from .npy files ------
    ~Data();
    
    int nx = 0;
    int nx_tot = 0;

    // =========== data ===========
    torch::Tensor rho;
    torch::Tensor nablaRho;
    torch::Tensor tau_tf;
    // semi-local descriptors
    torch::Tensor gamma;
    torch::Tensor p;
    torch::Tensor q;
    torch::Tensor tanhp;
    torch::Tensor tanhq;
    // non-local descriptors
    std::vector<torch::Tensor> gammanl = {};
    std::vector<torch::Tensor> pnl = {};
    std::vector<torch::Tensor> qnl = {};
    std::vector<torch::Tensor> xi = {};
    std::vector<torch::Tensor> tanhxi = {};
    std::vector<torch::Tensor> tanhxi_nl = {};
    std::vector<torch::Tensor> tanh_pnl = {};
    std::vector<torch::Tensor> tanh_qnl = {};
    std::vector<torch::Tensor> tanhp_nl = {};
    std::vector<torch::Tensor> tanhq_nl = {};
    // target
    torch::Tensor enhancement;
    torch::Tensor pauli;
    torch::Tensor enhancement_mean;
    torch::Tensor tau_mean; // mean Pauli energy
    torch::Tensor pauli_mean;

    // =========== label ===========
    bool load_gamma = false;
    bool load_p = false;
    bool load_q = false;
    bool load_tanhp = false;
    bool load_tanhq = false;
    bool* load_gammanl = nullptr;
    bool* load_pnl = nullptr;
    bool* load_qnl = nullptr;
    bool* load_xi = nullptr;
    bool* load_tanhxi = nullptr;
    bool* load_tanhxi_nl = nullptr;
    bool* load_tanh_pnl = nullptr;
    bool* load_tanh_qnl = nullptr;
    bool* load_tanhp_nl = nullptr;
    bool* load_tanhq_nl = nullptr;

    void load_data(Input &input, const int ndata, std::string *dir, const torch::Device device);
    torch::Tensor get_data(std::string parameter, const int ikernel);

  private:
    void init_label(Input &input);
    void init_data(const int nkernel, const int ndata, const int fftdim, const torch::Device device);
    void load_data_(Input &input, const int ndata, const int fftdim, std::string *dir);
    
    const double cTF = 3.0/10.0 * std::pow(3*std::pow(M_PI, 2.0), 2.0/3.0) * 2; // 10/3*(3*pi^2)^{2/3}, multiply by 2 to convert unit from Hartree to Ry, finally in Ry*Bohr^(-2)

  public:
    void loadTensor(std::string file,
                    std::vector<long unsigned int> cshape,
                    bool fortran_order,
                    std::vector<double> &container,
                    const int index,
                    const int fftdim,
                    torch::Tensor &data);
    // -------- dump Tensor into .npy files ---------
    void dumpTensor(const torch::Tensor &data, std::string filename, int nx);
    std::string file_name(std::string parameter, const int kernel_type, const double kernel_scaling);
};
#endif