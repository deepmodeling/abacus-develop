#include "./pauli_potential.h"

void PauliPotential::init(const Input &input,
                     const int ninput,
                     const std::vector<std::string> &descriptor_type,
                     const std::vector<int> &kernel_index)
{
    this->fftdim = input.fftdim;
    std::cout << "descriptor_type    " << descriptor_type << std::endl;
    std::cout << "kernel_index    " << kernel_index << std::endl;

    // this->descriptor_type = descriptor_type;
    // this->kernel_index = kernel_index;

    this->chi_xi = input.chi_xi;
    this->chi_p = input.chi_p;
    this->chi_q = input.chi_q;
    this->chi_pnl = input.chi_pnl;
    this->chi_qnl = input.chi_qnl;

    this->descriptor2kernel = {{"gamma", {}},
                               {"p", {}},
                               {"q", {}},
                               {"tanhp", {}},
                               {"tanhq", {}},
                               {"gammanl", {}},
                               {"pnl", {}},
                               {"qnl", {}},
                               {"xi", {}},
                               {"tanhxi", {}},
                               {"tanhxi_nl", {}},
                               {"tanh_pnl", {}},
                               {"tanh_qnl", {}},
                               {"tanhp_nl", {}},
                               {"tanhq_nl", {}}};
    this->descriptor2index = this->descriptor2kernel;

    for (int i = 0; i < ninput; ++i)
    {
        this->descriptor2kernel[descriptor_type[i]].push_back(kernel_index[i]);
        std::cout << "this->descriptor2kernel[descriptor_type[i]]    " << this->descriptor2kernel[descriptor_type[i]]
                  << std::endl;
        this->descriptor2index[descriptor_type[i]].push_back(i);
        std::cout << "this->descriptor2index[descriptor_type[i]]    " << this->descriptor2index[descriptor_type[i]]
                  << std::endl;
    }
    std::cout << "descriptor2index    " << descriptor2index << std::endl;
    std::cout << "descriptor2kernel    " << descriptor2kernel << std::endl;

    this->ml_gamma = this->descriptor2index["gamma"].size() > 0;
    this->ml_p = this->descriptor2index["p"].size() > 0;
    this->ml_q = this->descriptor2index["q"].size() > 0;
    this->ml_tanhp = this->descriptor2index["tanhp"].size() > 0;
    this->ml_tanhq = this->descriptor2index["tanhq"].size() > 0;
    this->ml_gammanl = this->descriptor2index["gammanl"].size() > 0;
    this->ml_pnl = this->descriptor2index["pnl"].size() > 0;
    this->ml_qnl = this->descriptor2index["qnl"].size() > 0;
    this->ml_xi = this->descriptor2index["xi"].size() > 0;
    this->ml_tanhxi = this->descriptor2index["tanhxi"].size() > 0;
    this->ml_tanhxi_nl = this->descriptor2index["tanhxi_nl"].size() > 0;
    this->ml_tanh_pnl = this->descriptor2index["tanh_pnl"].size() > 0;
    this->ml_tanh_qnl = this->descriptor2index["tanh_qnl"].size() > 0;
    this->ml_tanhp_nl = this->descriptor2index["tanhp_nl"].size() > 0;
    this->ml_tanhq_nl = this->descriptor2index["tanhq_nl"].size() > 0;
}

torch::Tensor PauliPotential::get_potential(const int istru,
                                       const Data &data,
                                       const torch::Tensor &F,
                                       const torch::Tensor &gradient,
                                       const Kernel *kernels,
                                       const Grid &grid)
{
    // Input::print("get potential begin");
    this->istru = istru;
    torch::Tensor potential = 5. / 3. * F;

    // semi-local potential terms
    if (this->ml_gamma) {
        potential += this->potGammaTerm(data.gamma[istru], gradient);
    }
    if (this->ml_p) {
        potential += this->potPTerm1(data.p[istru], gradient);
    }
    if (this->ml_q) {
        potential += this->potQTerm1(data.q[istru], gradient);
    }
    if (this->ml_xi) {
        potential += this->potXiTerm1(data.rho[istru], data.xi, gradient);
    }
    if (this->ml_tanhxi) {
        potential += this->potTanhxiTerm1(data.rho[istru], data.xi, data.tanhxi, gradient);
    }
    if (this->ml_tanhp) {
        potential += this->potTanhpTerm1(data.p[istru], data.tanhp[istru], gradient);
    }
    if (this->ml_tanhq) {
        potential += this->potTanhqTerm1(data.q[istru], data.tanhq[istru], gradient);
    }
    potential *= data.tau_tf[istru] / data.rho[istru];

    // non-local potential terms
    if (this->ml_gammanl) {
        potential += this->potGammanlTerm(data.rho[istru], data.gamma[istru], kernels, data.tau_tf[istru], gradient);
    }
    if (this->ml_p || this->ml_pnl) {
        potential += this->potPPnlTerm(data.rho[istru],
                                       data.nablaRho[istru],
                                       data.p[istru],
                                       kernels,
                                       data.tau_tf[istru],
                                       gradient,
                                       grid.fft_grid[istru]);
    }
    if (this->ml_q || this->ml_qnl) {
        potential += this->potQQnlTerm(data.rho[istru],
                                       data.q[istru],
                                       kernels,
                                       data.tau_tf[istru],
                                       gradient,
                                       grid.fft_gg[istru]);
    }
    if (this->ml_xi) {
        potential += this->potXinlTerm(data.rho[istru], kernels, data.tau_tf[istru], gradient);
    }
    if (this->ml_tanhxi) {
        potential += this->potTanhxinlTerm(data.rho[istru], data.tanhxi, kernels, data.tau_tf[istru], gradient);
    }
    if (this->ml_tanhxi_nl) {
        potential
            += this->potTanhxi_nlTerm(data.rho[istru], data.xi, data.tanhxi, kernels, data.tau_tf[istru], gradient);
    }
    if ((this->ml_tanhp || this->ml_tanhp_nl) && !this->ml_tanh_pnl) {
        potential += this->potTanhpTanhp_nlTerm(data.rho[istru],
                                                data.nablaRho[istru],
                                                data.p[istru],
                                                data.tanhp[istru],
                                                kernels,
                                                data.tau_tf[istru],
                                                gradient,
                                                grid.fft_grid[istru]);
    }
    if ((this->ml_tanhq || this->ml_tanhq_nl) && !this->ml_tanh_qnl) {
        potential += this->potTanhqTanhq_nlTerm(data.rho[istru],
                                                data.q[istru],
                                                data.tanhq[istru],
                                                kernels,
                                                data.tau_tf[istru],
                                                gradient,
                                                grid.fft_gg[istru]);
    }
    if (this->ml_tanh_pnl) {
        potential += this->potTanhpTanh_pnlTerm(data.rho[istru],
                                                data.nablaRho[istru],
                                                data.p[istru],
                                                data.tanhp[istru],
                                                data.tanh_pnl,
                                                kernels,
                                                data.tau_tf[istru],
                                                gradient,
                                                grid.fft_grid[istru]);
    }
    if (this->ml_tanh_qnl) {
        potential += this->potTanhqTanh_qnlTerm(data.rho[istru],
                                                data.q[istru],
                                                data.tanhq[istru],
                                                data.tanh_qnl,
                                                kernels,
                                                data.tau_tf[istru],
                                                gradient,
                                                grid.fft_gg[istru]);
    }

    // Input::print("get potential done");
    return potential;
}

torch::Tensor PauliPotential::potGammaTerm(const torch::Tensor &gamma, const torch::Tensor &gradient)
{
    // std::cout << "potGammaTerm" << std::endl;
    return 1. / 3. * gamma
           * gradient.index({"...", this->descriptor2index["gamma"][0]})
                 .reshape({this->fftdim, this->fftdim, this->fftdim});
}

torch::Tensor PauliPotential::potPTerm1(const torch::Tensor &p, const torch::Tensor &gradient)
{
    // std::cout << "potPTerm1" << std::endl;
    return -8. / 3. * p
           * gradient.index({"...", this->descriptor2index["p"][0]})
                 .reshape({this->fftdim, this->fftdim, this->fftdim});
}

torch::Tensor PauliPotential::potQTerm1(const torch::Tensor &q, const torch::Tensor &gradient)
{
    // std::cout << "potQTerm1" << std::endl;
    return -5. / 3. * q
           * gradient.index({"...", this->descriptor2index["q"][0]})
                 .reshape({this->fftdim, this->fftdim, this->fftdim});
}

torch::Tensor PauliPotential::potGammanlTerm(const torch::Tensor &rho,
                                        const torch::Tensor &gamma,
                                        const Kernel *kernels,
                                        // const torch::Tensor &kernel,
                                        const torch::Tensor &tauTF,
                                        const torch::Tensor &gradient)
{
    // std::cout << "potGmmamnlTerm" << std::endl;
    torch::Tensor result = torch::zeros_like(gamma);
    for (int ik = 0; ik < this->descriptor2kernel["gammanl"].size(); ++ik)
    {
        int d2k = this->descriptor2kernel["gammanl"][ik];
        int d2i = this->descriptor2index["gammanl"][ik];
        result += 1. / 3. * gamma / rho
                  * torch::real(
                      torch::fft::ifftn(torch::fft::fftn(gradient.index({"...", d2i})
                                                             .reshape({this->fftdim, this->fftdim, this->fftdim})
                                                         * tauTF)
                                        * kernels[d2k].kernel[istru]));
    }
    return result;
}

torch::Tensor PauliPotential::potPPnlTerm(const torch::Tensor &rho,
                                     const torch::Tensor &nablaRho,
                                     const torch::Tensor &p,
                                     const Kernel *kernels,
                                     //  const torch::Tensor &kernel,
                                     const torch::Tensor &tauTF,
                                     const torch::Tensor &gradient,
                                     const std::vector<torch::Tensor> &grid)
{
    // std::cout << "potPPnlTerm" << std::endl;
    torch::Tensor dFdpnl_nl = torch::zeros_like(p);
    for (int ik = 0; ik < this->descriptor2index["pnl"].size(); ++ik)
    {
        int d2k = this->descriptor2kernel["pnl"][ik];
        int d2i = this->descriptor2index["pnl"][ik];
        dFdpnl_nl
            += torch::real(torch::fft::ifftn(torch::fft::fftn(gradient.index({"...", d2i})
                                                                  .reshape({this->fftdim, this->fftdim, this->fftdim})
                                                              * tauTF)
                                             * kernels[d2k].kernel[istru]));
    }

    torch::Tensor temp = torch::zeros_like(nablaRho);
    for (int i = 0; i < 3; ++i)
    {
        temp[i] = (this->ml_p) ? -3. / 20.
                                     * gradient.index({"...", this->descriptor2index["p"][0]})
                                           .reshape({this->fftdim, this->fftdim, this->fftdim})
                                     * nablaRho[i] / rho * /*Ha to Ry*/ 2.
                               : torch::zeros_like(nablaRho[i]);
        if (this->ml_pnl) {
            temp[i] += -this->pqcoef * 2. * nablaRho[i] / torch::pow(rho, 8. / 3.) * dFdpnl_nl;
        }
    }
    // std::cout << torch::slice(temp[0][0][0], 0, 0, 10);
    torch::Tensor result = this->divergence(temp, grid);

    if (this->ml_pnl) {
        result += -8. / 3. * p / rho * dFdpnl_nl;
    }
    // std::cout << torch::slice(result[0][0], 0, 20) << std::endl;

    // std::cout << "potPPnlTerm done" << std::endl;
    return result;
}

torch::Tensor PauliPotential::potQQnlTerm(const torch::Tensor &rho,
                                     const torch::Tensor &q,
                                     const Kernel *kernels,
                                     //  const torch::Tensor &kernel,
                                     const torch::Tensor &tauTF,
                                     const torch::Tensor &gradient,
                                     const torch::Tensor &gg)
{
    // std::cout << "potQQnlTerm" << std::endl;
    torch::Tensor dFdqnl_nl = torch::zeros_like(q);
    for (int ik = 0; ik < this->descriptor2index["qnl"].size(); ++ik)
    {
        int d2k = this->descriptor2kernel["qnl"][ik];
        int d2i = this->descriptor2index["qnl"][ik];
        dFdqnl_nl
            = torch::real(torch::fft::ifftn(torch::fft::fftn(gradient.index({"...", d2i})
                                                                 .reshape({this->fftdim, this->fftdim, this->fftdim})
                                                             * tauTF)
                                            * kernels[d2k].kernel[istru]));
    }

    torch::Tensor temp = (this->ml_q) ? 3. / 40.
                                            * gradient.index({"...", this->descriptor2index["q"][0]})
                                                  .reshape({this->fftdim, this->fftdim, this->fftdim})
                                            * /*Ha2Ry*/ 2.
                                      : torch::zeros_like(q);
    if (this->ml_qnl) {
        temp += this->pqcoef / torch::pow(rho, 5. / 3.) * dFdqnl_nl;
    }
    torch::Tensor result = this->Laplacian(temp, gg);

    if (this->ml_qnl) {
        result += -5. / 3. * q / rho * dFdqnl_nl;
    }

    // std::cout << "potQQnlTerm done" << std::endl;
    return result;
}

torch::Tensor PauliPotential::potXiTerm1(const torch::Tensor &rho, const std::vector<torch::Tensor> &xi, const torch::Tensor &gradient)
{
    torch::Tensor result = torch::zeros_like(rho);
    for (int ik = 0; ik < this->descriptor2kernel["xi"].size(); ++ik)
    {
        int d2k = this->descriptor2kernel["xi"][ik];
        int d2i = this->descriptor2index["xi"][ik];
        result += -1. / 3. * xi[d2k][istru]
                  * gradient.index({"...", d2i})
                        .reshape({this->fftdim, this->fftdim, this->fftdim});
    }
    return result;
}

torch::Tensor PauliPotential::potTanhxiTerm1(const torch::Tensor &rho,
                                        const std::vector<torch::Tensor> &xi,
                                        const std::vector<torch::Tensor> &tanhxi,
                                        const torch::Tensor &gradient)
{
    torch::Tensor result = torch::zeros_like(rho);
    for (int ik = 0; ik < this->descriptor2kernel["tanhxi"].size(); ++ik)
    {
        int d2k = this->descriptor2kernel["tanhxi"][ik];
        int d2i = this->descriptor2index["tanhxi"][ik];
        result += -1. / 3. * xi[d2k][istru] * this->dtanh(tanhxi[d2k][istru], this->chi_xi[d2k])
                  * gradient.index({"...", d2i}).reshape({this->fftdim, this->fftdim, this->fftdim});
    }
    return result;
}

torch::Tensor PauliPotential::potTanhpTerm1(const torch::Tensor &p,
                                       const torch::Tensor &tanhp,
                                       const torch::Tensor &gradient)
{
    return -8. / 3. * p * this->dtanh(tanhp, this->chi_p)
           * gradient.index({"...", this->descriptor2index["tanhp"][0]})
                 .reshape({this->fftdim, this->fftdim, this->fftdim});
}

torch::Tensor PauliPotential::potTanhqTerm1(const torch::Tensor &q,
                                       const torch::Tensor &tanhq,
                                       const torch::Tensor &gradient)
{
    return -5. / 3. * q * this->dtanh(tanhq, this->chi_q)
           * gradient.index({"...", this->descriptor2index["tanhq"][0]})
                 .reshape({this->fftdim, this->fftdim, this->fftdim});
}

torch::Tensor PauliPotential::potXinlTerm(const torch::Tensor &rho,
                                     const Kernel *kernels,
                                     //  const torch::Tensor &kernel,
                                     const torch::Tensor &tauTF,
                                     const torch::Tensor &gradient)
{
    torch::Tensor result = torch::zeros_like(rho);
    for (int ik = 0; ik < this->descriptor2kernel["xi"].size(); ++ik)
    {
        int d2k = this->descriptor2kernel["xi"][ik];
        int d2i = this->descriptor2index["xi"][ik];
        result += 1. / 3. * torch::pow(rho, -2. / 3.)
                  * torch::real(
                      torch::fft::ifftn(torch::fft::fftn(gradient.index({"...", d2i})
                                                             .reshape({this->fftdim, this->fftdim, this->fftdim})
                                                         * tauTF * torch::pow(rho, -1. / 3.))
                                        * kernels[d2k].kernel[istru]));
    }
    return result;
}

torch::Tensor PauliPotential::potTanhxinlTerm(const torch::Tensor &rho,
                                         const std::vector<torch::Tensor> &tanhxi,
                                         const Kernel *kernels,
                                         //  const torch::Tensor &kernel,
                                         const torch::Tensor &tauTF,
                                         const torch::Tensor &gradient)
{
    torch::Tensor result = torch::zeros_like(rho);

    for (int ik = 0; ik < this->descriptor2kernel["tanhxi"].size(); ++ik)
    {
        int d2k = this->descriptor2kernel["tanhxi"][ik];
        int d2i = this->descriptor2index["tanhxi"][ik];
        result += 1. / 3. * torch::pow(rho, -2. / 3.)
                  * torch::real(torch::fft::ifftn(
                      torch::fft::fftn(gradient.index({"...", d2i}).reshape({this->fftdim, this->fftdim, this->fftdim})
                                       * this->dtanh(tanhxi[d2k][istru], this->chi_xi[d2k]) * tauTF
                                       * torch::pow(rho, -1. / 3.))
                      * kernels[d2k].kernel[istru]));
    }
    return result;
}

torch::Tensor PauliPotential::potTanhxi_nlTerm(const torch::Tensor &rho,
                                          const std::vector<torch::Tensor> &xi,
                                          const std::vector<torch::Tensor> &tanhxi,
                                          const Kernel *kernels,
                                          //   const torch::Tensor &kernel,
                                          const torch::Tensor &tauTF,
                                          const torch::Tensor &gradient)
{
    torch::Tensor result = torch::zeros_like(rho);
    for (int ik = 0; ik < this->descriptor2kernel["tanhxi_nl"].size(); ++ik)
    {
        int d2k = this->descriptor2kernel["tanhxi_nl"][ik];
        int d2i = this->descriptor2index["tanhxi_nl"][ik];
        torch::Tensor dFdxi
            = torch::real(torch::fft::ifftn(
                  torch::fft::fftn(tauTF
                                   * gradient.index({"...", d2i}).reshape({this->fftdim, this->fftdim, this->fftdim}))
                  * kernels[d2k].kernel[istru]))
              * this->dtanh(tanhxi[d2k][istru], this->chi_xi[d2k]) * torch::pow(rho, -1. / 3.);
        // std::cout << "tanhxi\n" << torch::slice(tanhxi[d2k][istru][0][0], 0, 0, 10) << std::endl;
        // std::cout << "gradient\n" << torch::slice(gradient, 0, 0, 10) << std::endl;
        // std::cout << "kernel\n" << torch::slice(kernels[d2k].kernel[istru][0][0], 0, 0, 10) << std::endl;
        result += 1. / 3. * torch::pow(rho, -2. / 3.)
                  * (-xi[d2k][istru] * dFdxi
                     + torch::real(torch::fft::ifftn(torch::fft::fftn(dFdxi) * kernels[d2k].kernel[istru])));
    }
    return result;
}

torch::Tensor PauliPotential::potTanhpTanh_pnlTerm(const torch::Tensor &rho,
                                              const torch::Tensor &nablaRho,
                                              const torch::Tensor &p,
                                              const torch::Tensor &tanhp,
                                              const std::vector<torch::Tensor> &tanh_pnl,
                                              const Kernel *kernels,
                                              //   const torch::Tensor &kernel,
                                              const torch::Tensor &tauTF,
                                              const torch::Tensor &gradient,
                                              const std::vector<torch::Tensor> &grid)
{
    torch::Tensor dFdpnl_nl = torch::zeros_like(tanhp);
    for (int ik = 0; ik < this->descriptor2kernel["tanh_pnl"].size(); ++ik)
    {
        int d2k = this->descriptor2kernel["tanh_pnl"][ik];
        int d2i = this->descriptor2index["tanh_pnl"][ik];
        dFdpnl_nl += torch::real(torch::fft::ifftn(
            torch::fft::fftn(gradient.index({"...", d2i})
                                 .reshape({this->fftdim, this->fftdim, this->fftdim})
                             * this->dtanh(tanh_pnl[d2k][istru], this->chi_pnl[d2k])
                             * tauTF)
            * kernels[d2k].kernel[istru]));
    }

    torch::Tensor temp = torch::zeros_like(nablaRho);
    for (int i = 0; i < 3; ++i)
    {
        temp[i] = (this->ml_tanhp) ? -3. / 20.
                                         * gradient.index({"...", this->descriptor2index["tanhp"][0]})
                                               .reshape({this->fftdim, this->fftdim, this->fftdim})
                                         * this->dtanh(tanhp, this->chi_p) * nablaRho[i] / rho * /*Ha to Ry*/ 2.
                                   : torch::zeros_like(nablaRho[i]);
        if (this->ml_tanh_pnl) {
            temp[i] += -this->pqcoef * 2. * nablaRho[i] / torch::pow(rho, 8. / 3.) * dFdpnl_nl;
        }
    }
    torch::Tensor result = this->divergence(temp, grid);

    if (this->ml_tanh_pnl) {
        result += -8. / 3. * p / rho * dFdpnl_nl;
    }

    return result;
}

torch::Tensor PauliPotential::potTanhqTanh_qnlTerm(const torch::Tensor &rho,
                                              const torch::Tensor &q,
                                              const torch::Tensor &tanhq,
                                              const std::vector<torch::Tensor> &tanh_qnl,
                                              const Kernel *kernels,
                                              //   const torch::Tensor &kernel,
                                              const torch::Tensor &tauTF,
                                              const torch::Tensor &gradient,
                                              const torch::Tensor &gg)
{
    torch::Tensor dFdqnl_nl = torch::zeros_like(tanhq);
    for (int ik = 0; ik < this->descriptor2kernel["tanh_qnl"].size(); ++ik)
    {
        int d2k = this->descriptor2kernel["tanh_qnl"][ik];
        int d2i = this->descriptor2index["tanh_qnl"][ik];
        dFdqnl_nl += torch::real(torch::fft::ifftn(
            torch::fft::fftn(gradient.index({"...", d2i})
                                 .reshape({this->fftdim, this->fftdim, this->fftdim})
                             * this->dtanh(tanh_qnl[d2k][istru], this->chi_qnl[d2k])
                             * tauTF)
            * kernels[d2k].kernel[istru]));
    }

    torch::Tensor temp = (this->ml_tanhq) ? 3. / 40.
                                                * gradient.index({"...", this->descriptor2index["tanhq"][0]})
                                                      .reshape({this->fftdim, this->fftdim, this->fftdim})
                                                * this->dtanh(tanhq, this->chi_q) * /*Ha2Ry*/ 2.
                                          : torch::zeros_like(q);
    if (this->ml_tanh_qnl) {
        temp += this->pqcoef / torch::pow(rho, 5. / 3.) * dFdqnl_nl;
    }
    torch::Tensor result = this->Laplacian(temp, gg);

    if (this->ml_tanh_qnl) {
        result += -5. / 3. * q / rho * dFdqnl_nl;
    }

    return result;
}

torch::Tensor PauliPotential::potTanhpTanhp_nlTerm(const torch::Tensor &rho,
                                              const torch::Tensor &nablaRho,
                                              const torch::Tensor &p,
                                              const torch::Tensor &tanhp,
                                              const Kernel *kernels,
                                              //   const torch::Tensor &kernel,
                                              const torch::Tensor &tauTF,
                                              const torch::Tensor &gradient,
                                              const std::vector<torch::Tensor> &grid)
{
    torch::Tensor dFdpnl_nl = torch::zeros_like(tanhp);
    for (int ik = 0; ik < this->descriptor2kernel["tanhp_nl"].size(); ++ik)
    {
        int d2k = this->descriptor2kernel["tanhp_nl"][ik];
        int d2i = this->descriptor2index["tanhp_nl"][ik];
        dFdpnl_nl += torch::real(torch::fft::ifftn(
                         torch::fft::fftn(gradient.index({"...", d2i})
                                              .reshape({this->fftdim, this->fftdim, this->fftdim})
                                          * tauTF)
                         * kernels[d2k].kernel[istru]))
                     * this->dtanh(tanhp, this->chi_p);
    }

    torch::Tensor temp = torch::zeros_like(nablaRho);
    for (int i = 0; i < 3; ++i)
    {
        temp[i] = (this->ml_tanhp) ? -3. / 20.
                                         * gradient.index({"...", this->descriptor2index["tanhp"][0]})
                                               .reshape({this->fftdim, this->fftdim, this->fftdim})
                                         * this->dtanh(tanhp, this->chi_p) * nablaRho[i] / rho * /*Ha to Ry*/ 2.
                                   : torch::zeros_like(nablaRho[i]);
        if (this->ml_tanhp_nl) {
            temp[i] += -this->pqcoef * 2. * nablaRho[i] / torch::pow(rho, 8. / 3.) * dFdpnl_nl;
        }
    }
    torch::Tensor result = this->divergence(temp, grid);

    if (this->ml_tanhp_nl) {
        result += -8. / 3. * p / rho * dFdpnl_nl;
    }

    return result;
}

torch::Tensor PauliPotential::potTanhqTanhq_nlTerm(const torch::Tensor &rho,
                                              const torch::Tensor &q,
                                              const torch::Tensor &tanhq,
                                              const Kernel *kernels,
                                              //   const torch::Tensor &kernel,
                                              const torch::Tensor &tauTF,
                                              const torch::Tensor &gradient,
                                              const torch::Tensor &gg)
{
    torch::Tensor dFdqnl_nl = torch::zeros_like(tanhq);
    for (int ik = 0; ik < this->descriptor2kernel["tanhq_nl"].size(); ++ik)
    {
        int d2k = this->descriptor2kernel["tanhq_nl"][ik];
        int d2i = this->descriptor2index["tanhq_nl"][ik];
        dFdqnl_nl += torch::real(torch::fft::ifftn(
                         torch::fft::fftn(gradient.index({"...", d2i})
                                              .reshape({this->fftdim, this->fftdim, this->fftdim})
                                          * tauTF)
                         * kernels[d2k].kernel[istru]))
                     * this->dtanh(tanhq, this->chi_q);
    }

    torch::Tensor temp = (this->ml_tanhq) ? 3. / 40.
                                                * gradient.index({"...", this->descriptor2index["tanhq"][0]})
                                                      .reshape({this->fftdim, this->fftdim, this->fftdim})
                                                * this->dtanh(tanhq, this->chi_q) * /*Ha2Ry*/ 2.
                                          : torch::zeros_like(q);
    if (this->ml_tanhq_nl) {
        temp += this->pqcoef / torch::pow(rho, 5. / 3.) * dFdqnl_nl;
    }
    torch::Tensor result = this->Laplacian(temp, gg);

    if (this->ml_tanhq_nl) {
        result += -5. / 3. * q / rho * dFdqnl_nl;
    }

    return result;
}

torch::Tensor PauliPotential::divergence(const torch::Tensor &input, const std::vector<torch::Tensor> &grid)
{
    torch::Tensor result = torch::zeros_like(input[0]);
    // torch::Tensor img = torch::tensor({1.0j});
    // for (int i = 0; i < 3; ++i)
    // {
    //     result += torch::real(torch::fft::ifftn(torch::fft::fftn(input[i]) * grid[i] * img));
    // }
    for (int i = 0; i < 3; ++i)
    {
        result -= torch::imag(torch::fft::ifftn(torch::fft::fftn(input[i]) * grid[i]));
    }
    return result;
}

torch::Tensor PauliPotential::Laplacian(const torch::Tensor &input, const torch::Tensor &gg)
{
    return torch::real(torch::fft::ifftn(torch::fft::fftn(input) * -gg));
}

torch::Tensor PauliPotential::dtanh(const torch::Tensor &tanhx, const double chi)
{
    return (torch::ones_like(tanhx) - tanhx * tanhx) * chi;
    // return (1. - tanhx * tanhx) * chi;
}