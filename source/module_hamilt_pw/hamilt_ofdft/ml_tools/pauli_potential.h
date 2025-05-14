#ifndef PAULI_POTENTIAL_H
#define PAULI_POTENTIAL_H

#include <torch/torch.h>
#include "./input.h"
#include "./data.h"
#include "./kernel.h"
#include "./grid.h"

class PauliPotential{

public:
    void init(const Input &input, const int ninput, const std::vector<std::string> &descriptor_type, const std::vector<int> &kernel_index);

    int fftdim = 0;
    int istru = 0;

    std::map<std::string, std::vector<int>> descriptor2kernel = {};
    std::map<std::string, std::vector<int>> descriptor2index = {};
    
    // semi-local descriptors
    bool ml_gamma = false;
    bool ml_p = false;
    bool ml_q = false;
    bool ml_tanhp = false;
    bool ml_tanhq = false;
    // non-local descriptors
    bool ml_gammanl = false;
    bool ml_pnl = false;
    bool ml_qnl = false;
    bool ml_xi = false;
    bool ml_tanhxi = false;
    bool ml_tanhxi_nl = false;
    bool ml_tanh_pnl = false;
    bool ml_tanh_qnl = false;
    bool ml_tanhp_nl = false;
    bool ml_tanhq_nl = false;

    double chi_p = 1.;
    double chi_q = 1.;
    double* chi_xi = nullptr;
    double* chi_pnl = nullptr;
    double* chi_qnl = nullptr;

    torch::Tensor get_potential(
        const int istru,
        const Data &data,
        const torch::Tensor &F,
        const torch::Tensor &gradient,
        const Kernel *kernels,
        const Grid &grid
    );
private:

    torch::Tensor potGammaTerm(
        const torch::Tensor &gamma,
        const torch::Tensor &gradient
    );
    torch::Tensor potPTerm1(
        const torch::Tensor &p,
        const torch::Tensor &gradient
    );
    torch::Tensor potQTerm1(
        const torch::Tensor &q,
        const torch::Tensor &gradient
    );
    torch::Tensor potGammanlTerm(
        const torch::Tensor &rho,
        const torch::Tensor &gamma,
        const Kernel *kernels,
        // const torch::Tensor &kernel,
        const torch::Tensor &tauTF,
        const torch::Tensor &gradient
    );
    torch::Tensor potPPnlTerm(
        const torch::Tensor &rho,
        const torch::Tensor &nablaRho,
        const torch::Tensor &p,
        const Kernel *kernels,
        // const torch::Tensor &kernel,
        const torch::Tensor &tauTF,
        const torch::Tensor &gradient,
        const std::vector<torch::Tensor> &grid
    );
    torch::Tensor potQQnlTerm(
        const torch::Tensor &rho,
        const torch::Tensor &q,
        const Kernel *kernels,
        // const torch::Tensor &kernel,
        const torch::Tensor &tauTF,
        const torch::Tensor &gradient,
        const torch::Tensor &gg    
    );


    torch::Tensor potXiTerm1(
        const torch::Tensor &rho,
        const std::vector<torch::Tensor> &xi,
        const torch::Tensor &gradient
    );
    torch::Tensor potTanhxiTerm1(
        const torch::Tensor &rho,
        const std::vector<torch::Tensor> &xi,
        const std::vector<torch::Tensor> &tanhxi,
        const torch::Tensor &gradient
    );
    torch::Tensor potTanhpTerm1(
        const torch::Tensor &p,
        const torch::Tensor &tanhp,
        const torch::Tensor &gradient
    );
    torch::Tensor potTanhqTerm1(
        const torch::Tensor &q,
        const torch::Tensor &tanhq,
        const torch::Tensor &gradient
    );
    torch::Tensor potXinlTerm(
        const torch::Tensor &rho,
        const Kernel *kernels,
        // const torch::Tensor &kernel,
        const torch::Tensor &tauTF,
        const torch::Tensor &gradient
    );
    torch::Tensor potTanhxinlTerm(
        const torch::Tensor &rho,
        const std::vector<torch::Tensor> &tanhxi,
        const Kernel *kernels,
        // const torch::Tensor &kernel,
        const torch::Tensor &tauTF,
        const torch::Tensor &gradient
    );
    torch::Tensor potTanhxi_nlTerm(
        const torch::Tensor &rho,
        const std::vector<torch::Tensor> &xi,
        const std::vector<torch::Tensor> &tanhxi,
        const Kernel *kernels,
        // const torch::Tensor &kernel,
        const torch::Tensor &tauTF,
        const torch::Tensor &gradient
    );
    torch::Tensor potTanhpTanh_pnlTerm(
        const torch::Tensor &rho,
        const torch::Tensor &nablaRho,
        const torch::Tensor &p,
        const torch::Tensor &tanhp,
        const std::vector<torch::Tensor> &tanh_pnl,
        const Kernel *kernels,
        // const torch::Tensor &kernel,
        const torch::Tensor &tauTF,
        const torch::Tensor &gradient,
        const std::vector<torch::Tensor> &grid
    );
    torch::Tensor potTanhqTanh_qnlTerm(
        const torch::Tensor &rho,
        const torch::Tensor &q,
        const torch::Tensor &tanhq,
        const std::vector<torch::Tensor> &tanh_qnl,
        const Kernel *kernels,
        // const torch::Tensor &kernel,
        const torch::Tensor &tauTF,
        const torch::Tensor &gradient,
        const torch::Tensor &gg
    );
    torch::Tensor potTanhpTanhp_nlTerm(
        const torch::Tensor &rho,
        const torch::Tensor &nablaRho,
        const torch::Tensor &p,
        const torch::Tensor &tanhp,
        const Kernel *kernels,
        // const torch::Tensor &kernel,
        const torch::Tensor &tauTF,
        const torch::Tensor &gradient,
        const std::vector<torch::Tensor> &grid
    );
    torch::Tensor potTanhqTanhq_nlTerm(
        const torch::Tensor &rho,
        const torch::Tensor &q,
        const torch::Tensor &tanhq,
        const Kernel *kernels,
        // const torch::Tensor &kernel,
        const torch::Tensor &tauTF,
        const torch::Tensor &gradient,
        const torch::Tensor &gg
    );

    // Tools for getting potential
    torch::Tensor divergence(
        const torch::Tensor &input,
        const std::vector<torch::Tensor> &grid
    );
    torch::Tensor Laplacian(
        const torch::Tensor &input,
        const torch::Tensor &gg
    );
    torch::Tensor dtanh(
        const torch::Tensor &tanhx,
        const double chi
    );

    const double cTF = 3.0/10.0 * std::pow(3*std::pow(M_PI, 2.0), 2.0/3.0) * 2; // 10/3*(3*pi^2)^{2/3}, multiply by 2 to convert unit from Hartree to Ry, finally in Ry*Bohr^(-2)
    const double pqcoef = 1.0 / (4.0 * std::pow(3*std::pow(M_PI, 2.0), 2.0/3.0)); // coefficient of p and q
};
#endif