#include "./kernel.h"

void Kernel::fill_kernel(const int fftdim,
                         const int ndata,
                         const torch::Tensor &rho,
                         const double *volume,
                         const std::string *cell,
                         const torch::Device device,
                         const std::vector<torch::Tensor> &fft_gg)
{
    double rho0 = 0.;
    double tkF = 0.;
    double eta = 0.;
    this->kernel = std::vector<torch::Tensor>(ndata);

    // read in the kernel
    if (this->kernel_type == 3 || this->kernel_type == 4) // 3 for TKK Al, and 4 for TKK Si
    {
        this->read_kernel(fftdim, ndata, rho, volume, cell, device, fft_gg);
    }
    else
    {
        for (int id = 0; id < ndata; ++id)
        {
            rho0 = torch::sum(rho[id]).item<double>() / std::pow(fftdim, 3);
            std::cout << "There are " << rho0 * volume[id] << " electrons in " << cell[id] << " strcture." << std::endl;
            tkF = 2. * std::pow(3. * std::pow(M_PI, 2) * rho0, 1. / 3.);
            this->kernel[id] = torch::zeros({fftdim, fftdim, fftdim}).to(device);
            for (int ix = 0; ix < fftdim; ++ix)
            {
                for (int iy = 0; iy < fftdim; ++iy)
                {
                    for (int iz = 0; iz < fftdim; ++iz)
                    {
                        eta = sqrt(fft_gg[id][ix][iy][iz].item<double>()) / tkF;
                        eta = eta * this->scaling;
                        if (this->kernel_type == 1)
                        {
                            this->kernel[id][ix][iy][iz] = this->wt_kernel(eta);
                            // this->kernel[id][ix][iy][iz] = std::pow(1. / this->scaling, 3) * this->wt_kernel(eta);
                            // if (ix == 26 && iy == 0 && iz == 26)
                            // {
                            //     std::cout << "kernel1    " << this->kernel[id][ix][iy][iz] << std::endl;
                            //     std::cout << "eta        " << eta << std::endl;
                            //     std::cout << "tkF        " << tkF << std::endl;
                            // }
                        }
                        else if (this->kernel_type == 2)
                        {
                            this->kernel[id][ix][iy][iz] = this->yukawa_kernel(eta, this->yukawa_alpha);
                        }
                    }
                }
            }
        }
    }
    std::cout << "Fill kernel done" << std::endl;
}

double Kernel::wt_kernel(double eta, double tf_weight, double vw_weight)
{
    if (eta < 0.) 
    {
        return 0.;
    }
    // limit for small eta
    else if (eta < 1e-10)
    {
        return 1. - tf_weight + eta * eta * (1./3. - 3. * vw_weight);
    }
    // around the singularity
    else if (std::abs(eta - 1.) < 1e-10)
    {
        return 2. - tf_weight - 3. * vw_weight + 20. * (eta - 1);
    }
    // Taylor expansion for high eta
    else if (eta > 3.65)
    {
        double eta2 = eta * eta;
        double invEta2 = 1. / eta2;
        double LindG = 3. * (1. - vw_weight) * eta2 
                        -tf_weight-0.6 
                        + invEta2 * (-0.13714285714285712 
                        + invEta2 * (-6.39999999999999875E-2
                        + invEta2 * (-3.77825602968460128E-2
                        + invEta2 * (-2.51824061652633074E-2
                        + invEta2 * (-1.80879839616166146E-2
                        + invEta2 * (-1.36715733124818332E-2
                        + invEta2 * (-1.07236045520990083E-2
                        + invEta2 * (-8.65192783339199453E-3 
                        + invEta2 * (-7.1372762502456763E-3
                        + invEta2 * (-5.9945117538835746E-3
                        + invEta2 * (-5.10997527675418131E-3 
                        + invEta2 * (-4.41060829979912465E-3 
                        + invEta2 * (-3.84763737842981233E-3 
                        + invEta2 * (-3.38745061493813488E-3 
                        + invEta2 * (-3.00624946457977689E-3)))))))))))))));
        return LindG;
    }
    else
    {
        return 1. / (0.5 + 0.25 * (1. - eta * eta) / eta * std::log((1 + eta)/std::abs(1 - eta)))
                 - 3. * vw_weight * eta * eta - tf_weight;
    }
}

double Kernel::yukawa_kernel(double eta, double alpha)
{
    return (eta == 0 && alpha == 0) ? 0. : M_PI / (eta * eta + alpha * alpha / 4.);
}

// Read kernel from file
void Kernel::read_kernel(const int fftdim,
                         const int ndata,
                         const torch::Tensor &rho,
                         const double *volume,
                         const std::string *cell,
                         const torch::Device device,
                         const std::vector<torch::Tensor> &fft_gg)
{
    std::ifstream ifs(kernel_file.c_str(), std::ios::in);

    if (!ifs)
    {
        std::cout << " Can't find the kernel file " << kernel_file << std::endl;
        exit(0);
    }

    std::cout << "Read WT kernel from " << kernel_file << std::endl;

    int kineType = 0;
    double kF_in = 0.;
    double rho0_in = 0.;
    int nq_in = 0;
    double maxEta_in = 0.;

    ifs >> kineType;
    ifs >> kF_in;
    ifs >> rho0_in;
    ifs >> nq_in;

    double *eta_in = new double[nq_in];
    double *w0_in = new double[nq_in];

    for (int iq = 0; iq < nq_in; ++iq)
    {
        ifs >> eta_in[iq] >> w0_in[iq];
    }

    maxEta_in = eta_in[nq_in-1];

    double rho0 = 0.;
    double tkF = 0.;
    double eta = 0.;
    for (int id = 0; id < ndata; ++id)
    {
        rho0 = torch::sum(rho[id]).item<double>() / std::pow(fftdim, 3);
        std::cout << "There are " << rho0 * volume[id] << " electrons in " << cell[id] << " strcture." << std::endl;
        tkF = 2. * std::pow(3. * std::pow(M_PI, 2) * rho0, 1. / 3.);
        this->kernel[id] = torch::zeros({fftdim, fftdim, fftdim}).to(device);

        double eta = 0.;
        double maxEta = 0.;
        int ind1 = 0;
        int ind2 = 0;
        int ind_mid = 0;
        double fac1 = 0.;
        double fac2 = 0.;

        for (int ix = 0; ix < fftdim; ++ix)
        {
            for (int iy = 0; iy < fftdim; ++iy)
            {
                for (int iz = 0; iz < fftdim; ++iz)
                {
                    eta = sqrt(fft_gg[id][ix][iy][iz].item<double>()) / tkF;
                    eta = eta * this->scaling;
                    maxEta = std::max(eta, maxEta);

                    if (eta <= eta_in[0]) {
                        this->kernel[id][ix][iy][iz] = w0_in[0];
                    } else if (eta > maxEta_in) {
                        this->kernel[id][ix][iy][iz] = w0_in[nq_in-1];
                    } else
                    {
                        ind1 = 1;
                        ind2 = nq_in;
                        while (ind1 < ind2 - 1)
                        {
                            ind_mid = (ind1 + ind2)/2;
                            if (eta > eta_in[ind_mid])
                            {
                                ind1 = ind_mid;
                            }
                            else
                            {
                                ind2 = ind_mid;
                            }
                        }
                        fac1 = (eta_in[ind2] - eta)/(eta_in[ind2] - eta_in[ind1]);
                        fac2 = (eta - eta_in[ind1])/(eta_in[ind2] - eta_in[ind1]);
                        this->kernel[id][ix][iy][iz] = fac1 * w0_in[ind1] + fac2 * w0_in[ind2];
                        // this->kernel[id][ix][iy][iz] *= std::pow(1. / this->scaling, 3);
                    }
                }
            }
        }
        if (maxEta > maxEta_in) { std::cout << "Please increase the maximal eta value in KEDF kernel file" << std::endl;
}
    }


    delete[] eta_in;
    delete[] w0_in;
}