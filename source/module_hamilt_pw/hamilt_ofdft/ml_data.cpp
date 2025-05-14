#ifdef __MLKEDF

#include "ml_data.h"

#include "npy.hpp"
#include "module_elecstate/module_charge/symmetry_rho.h"

void ML_data::set_para(
    const int &nx,
    const double &nelec, 
    const double &tf_weight, 
    const double &vw_weight,
    const double &chi_p,
    const double &chi_q,
    const std::vector<double> &chi_xi,
    const std::vector<double> &chi_pnl,
    const std::vector<double> &chi_qnl,
    const int &nkernel,
    const std::vector<int> &kernel_type,
    const std::vector<double> &kernel_scaling,
    const std::vector<double> &yukawa_alpha,
    const std::vector<std::string> &kernel_file,
    const double &omega,
    ModulePW::PW_Basis *pw_rho
)
{
    this->nx = nx;
    this->nkernel = nkernel;
    this->chi_p = chi_p;
    this->chi_q = chi_q;
    this->chi_xi = chi_xi;
    this->chi_pnl = chi_pnl;
    this->chi_qnl = chi_qnl;

    this->kernel_type = kernel_type;
    this->kernel_scaling = kernel_scaling;
    this->yukawa_alpha = yukawa_alpha;
    this->kernel_file = kernel_file;
    std::cout << "nkernel = " << nkernel << std::endl;

    if (PARAM.inp.of_wt_rho0 != 0)
    {
        this->rho0 = PARAM.inp.of_wt_rho0;
    }
    else
    {
        this->rho0 = 1./omega * nelec;
    }

    this->kF = std::pow(3. * std::pow(ModuleBase::PI, 2) * this->rho0, 1./3.);
    this->tkF = 2. * this->kF;

    this->kernel = new double*[this->nkernel];
    for (int ik = 0; ik < this->nkernel; ++ik)
    {
        // delete[] this->kernel[ik];
        this->kernel[ik] = new double[pw_rho->npw];
        if (this->kernel_type[ik] == 3 || this->kernel_type[ik] == 4) // 3 for TKK Al, and 4 for TKK Si
        {
            this->read_kernel(this->kernel_file[ik], this->kernel_scaling[ik], pw_rho, this->kernel[ik]);
        }
        else
        {
            double eta = 0.;
            for (int ip = 0; ip < pw_rho->npw; ++ip)
            {
                eta = sqrt(pw_rho->gg[ip]) * pw_rho->tpiba / this->tkF * this->kernel_scaling[ik];
                if (this->kernel_type[ik] == 1)
                {
                    this->kernel[ik][ip] = this->MLkernel(eta, tf_weight, vw_weight);
                }
                else if (this->kernel_type[ik] == 2)
                {
                    this->kernel[ik][ip] = this->MLkernel_yukawa(eta, this->yukawa_alpha[ik]);
                }
            }
        }
    }
}

void ML_data::generateTrainData_WT(
    const double * const *prho, 
    KEDF_WT &wt, 
    KEDF_TF &tf, 
    ModulePW::PW_Basis *pw_rho,
    const double* veff    
)
{
    std::vector<std::vector<double>> nablaRho(3, std::vector<double>(this->nx, 0.));

    this->generate_descriptor(prho, pw_rho, nablaRho);

    std::vector<double> container(this->nx);
    const long unsigned cshape[] = {(long unsigned) this->nx}; // shape of container and containernl
    
    // enhancement factor of Pauli potential
    if (PARAM.inp.of_kinetic == "wt")
    {
        this->getF_WT(wt, tf, prho, pw_rho, container);
        npy::SaveArrayAsNumpy("enhancement.npy", false, 1, cshape, container);

        // Pauli potential
        this->getPauli_WT(wt, tf, prho, pw_rho, container);
        npy::SaveArrayAsNumpy("pauli.npy", false, 1, cshape, container);
    }

    for (int ir = 0; ir < this->nx; ++ir)
    {
        container[ir] = veff[ir];
    }
    npy::SaveArrayAsNumpy("veff.npy", false, 1, cshape, container);
}

void ML_data::generateTrainData_KS(
    psi::Psi<std::complex<double>> *psi,
    elecstate::ElecState *pelec,
    ModulePW::PW_Basis_K *pw_psi,
    ModulePW::PW_Basis *pw_rho,
    UnitCell& ucell,
    const double* veff
)
{
    std::vector<std::vector<double>> nablaRho(3, std::vector<double>(this->nx, 0.));

    this->generate_descriptor(pelec->charge->rho, pw_rho, nablaRho);

    std::vector<double> container(this->nx);
    std::vector<double> containernl(this->nx);

    const long unsigned cshape[] = {(long unsigned) this->nx}; // shape of container and containernl
    // enhancement factor of Pauli energy, and Pauli potential
    this->getF_KS1(psi, pelec, pw_psi, pw_rho, ucell, nablaRho, container, containernl);

    Symmetry_rho srho;

    Charge* ptempRho = new Charge();
    ptempRho->nspin = PARAM.inp.nspin;
    ptempRho->nrxx = this->nx;
    ptempRho->rho_core = pelec->charge->rho_core;
    ptempRho->rho = new double*[1];
    ptempRho->rho[0] = new double[this->nx];
    ptempRho->rhog = new std::complex<double>*[1];
    ptempRho->rhog[0] = new std::complex<double>[pw_rho->npw];

    for (int ir = 0; ir < this->nx; ++ir){
        ptempRho->rho[0][ir] = container[ir];
    }
    srho.begin(0, *ptempRho, pw_rho, ucell.symm);
    for (int ir = 0; ir < this->nx; ++ir){
        container[ir] = ptempRho->rho[0][ir];
    }

    for (int ir = 0; ir < this->nx; ++ir){
        ptempRho->rho[0][ir] = containernl[ir];
    }
    srho.begin(0, *ptempRho, pw_rho, ucell.symm);
    for (int ir = 0; ir < this->nx; ++ir){
        containernl[ir] = ptempRho->rho[0][ir];
    }

    npy::SaveArrayAsNumpy("enhancement.npy", false, 1, cshape, container);
    npy::SaveArrayAsNumpy("pauli.npy", false, 1, cshape, containernl);

    // enhancement factor of Pauli energy, and Pauli potential
    this->getF_KS2(psi, pelec, pw_psi, pw_rho, ucell, container, containernl);

    for (int ir = 0; ir < this->nx; ++ir){
        ptempRho->rho[0][ir] = container[ir];
    }
    srho.begin(0, *ptempRho, pw_rho, ucell.symm);
    for (int ir = 0; ir < this->nx; ++ir){
        container[ir] = ptempRho->rho[0][ir];
    }

    for (int ir = 0; ir < this->nx; ++ir){
        ptempRho->rho[0][ir] = containernl[ir];
    }
    srho.begin(0, *ptempRho, pw_rho, ucell.symm);
    for (int ir = 0; ir < this->nx; ++ir){
        containernl[ir] = ptempRho->rho[0][ir];
    }

    npy::SaveArrayAsNumpy("enhancement2.npy", false, 1, cshape, container);
    npy::SaveArrayAsNumpy("pauli2.npy", false, 1, cshape, containernl);

    for (int ir = 0; ir < this->nx; ++ir)
    {
        container[ir] = veff[ir];
    }
    npy::SaveArrayAsNumpy("veff.npy", false, 1, cshape, container);

    delete ptempRho;
}

void ML_data::generate_descriptor(
    const double * const *prho, 
    ModulePW::PW_Basis *pw_rho,
    std::vector<std::vector<double>> &nablaRho
)
{
    // container which will contain gamma, p, q in turn
    std::vector<double> container(this->nx);
    std::vector<double> new_container(this->nx);
    // container contains gammanl, pnl, qnl in turn
    std::vector<double> containernl(this->nx);
    std::vector<double> new_containernl(this->nx);

    const long unsigned cshape[] = {(long unsigned) this->nx}; // shape of container and containernl

    // rho
    std::vector<double> rho(this->nx);
    for (int ir = 0; ir < this->nx; ++ir){
        rho[ir] = prho[0][ir];
    }
    npy::SaveArrayAsNumpy("rho.npy", false, 1, cshape, rho);

    // gamma
    this->getGamma(prho, container);
    npy::SaveArrayAsNumpy("gamma.npy", false, 1, cshape, container);

    for (int ik = 0; ik < this->nkernel; ++ik)
    {
        int ktype = this->kernel_type[ik];
        double kscaling = this->kernel_scaling[ik];

        // gamma_nl
        this->getGammanl(ik, container, pw_rho, containernl);
        npy::SaveArrayAsNumpy(this->file_name("gammanl", ktype, kscaling), false, 1, cshape, containernl);

        // xi = gamma_nl/gamma
        this->getXi(container, containernl, new_container);
        npy::SaveArrayAsNumpy(this->file_name("xi", ktype, kscaling), false, 1, cshape, new_container);

        // tanhxi = tanh(xi)
        this->getTanhXi(ik, container, containernl, new_container);
        npy::SaveArrayAsNumpy(this->file_name("tanhxi", ktype, kscaling), false, 1, cshape, new_container);

        // (tanhxi)_nl
        this->getTanhXi_nl(ik, new_container, pw_rho, new_containernl);
        npy::SaveArrayAsNumpy(this->file_name("tanhxi_nl", ktype, kscaling), false, 1, cshape, new_containernl);
    }

    // nabla rho
    this->getNablaRho(prho, pw_rho, nablaRho);
    npy::SaveArrayAsNumpy("nablaRhox.npy", false, 1, cshape, nablaRho[0]);
    npy::SaveArrayAsNumpy("nablaRhoy.npy", false, 1, cshape, nablaRho[1]);
    npy::SaveArrayAsNumpy("nablaRhoz.npy", false, 1, cshape, nablaRho[2]);

    // p
    this->getP(prho, pw_rho, nablaRho, container);
    npy::SaveArrayAsNumpy("p.npy", false, 1, cshape, container);

    for (int ik = 0; ik < this->nkernel; ++ik)
    {
        int ktype = this->kernel_type[ik];
        double kscaling = this->kernel_scaling[ik];

        // p_nl
        this->getPnl(ik, container, pw_rho, containernl);
        npy::SaveArrayAsNumpy(this->file_name("pnl", ktype, kscaling), false, 1, cshape, containernl);

        // tanh(p_nl)
        this->getTanh_Pnl(ik, containernl, new_containernl);
        npy::SaveArrayAsNumpy(this->file_name("tanh_pnl", ktype, kscaling), false, 1, cshape, new_containernl);

        // tanh(p)
        this->getTanhP(container, new_container);
        npy::SaveArrayAsNumpy("tanhp.npy", false, 1, cshape, new_container);

        // tanh(p)_nl
        this->getTanhP_nl(ik, new_container, pw_rho, new_containernl);
        npy::SaveArrayAsNumpy(this->file_name("tanhp_nl", ktype, kscaling), false, 1, cshape, new_containernl);
    }

    // q
    this->getQ(prho, pw_rho, container);
    npy::SaveArrayAsNumpy("q.npy", false, 1, cshape, container);

    for (int ik = 0; ik < this->nkernel; ++ik)
    {
        int ktype = this->kernel_type[ik];
        double kscaling = this->kernel_scaling[ik];

        // q_nl
        this->getQnl(ik, container, pw_rho, containernl);
        npy::SaveArrayAsNumpy(this->file_name("qnl", ktype, kscaling), false, 1, cshape, containernl);

        // tanh(q_nl)
        this->getTanh_Qnl(ik, containernl, new_containernl);
        npy::SaveArrayAsNumpy(this->file_name("tanh_qnl", ktype, kscaling), false, 1, cshape, new_containernl);

        // tanh(q)
        this->getTanhQ(container, new_container);
        npy::SaveArrayAsNumpy("tanhq.npy", false, 1, cshape, new_container);

        // tanh(q)_nl
        this->getTanhQ_nl(ik, new_container, pw_rho, new_containernl);
        npy::SaveArrayAsNumpy(this->file_name("tanhq_nl", ktype, kscaling), false, 1, cshape, new_containernl);
    }
}

double ML_data::MLkernel(double eta, double tf_weight, double vw_weight)
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

double ML_data::MLkernel_yukawa(double eta, double alpha)
{
    return (eta == 0 && alpha == 0) ? 0. : M_PI / (eta * eta + alpha * alpha / 4.);
}

// Read kernel from file
void ML_data::read_kernel(const std::string &fileName, const double& scaling, ModulePW::PW_Basis *pw_rho, double* kernel_)
{
    std::ifstream ifs(fileName.c_str(), std::ios::in);

    GlobalV::ofs_running << "Read WT kernel from " << fileName << std::endl;
    if (!ifs) ModuleBase::WARNING_QUIT("ml_data.cpp", "The kernel file not found");

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

    double eta = 0.;
    double maxEta = 0.;
    int ind1 = 0;
    int ind2 = 0;
    int ind_mid = 0;
    double fac1 = 0.;
    double fac2 = 0.;
    for (int ig = 0; ig < pw_rho->npw; ++ig)
    {
        eta = sqrt(pw_rho->gg[ig]) * pw_rho->tpiba / this->tkF;
        eta = eta * scaling;
        maxEta = std::max(eta, maxEta);

        if (eta <= eta_in[0])
            kernel_[ig] = w0_in[0];
        else if (eta > maxEta_in)
            kernel_[ig] = w0_in[nq_in-1];
        else
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
            kernel_[ig] = fac1 * w0_in[ind1] + fac2 * w0_in[ind2];
        }
    }

    if (maxEta > maxEta_in) ModuleBase::WARNING("kedf_wt.cpp", "Please increase the maximal eta value in KEDF kernel file");

    delete[] eta_in;
    delete[] w0_in;
    ModuleBase::GlobalFunc::DONE(GlobalV::ofs_running, "FILL WT KERNEL");
}

void ML_data::multiKernel(const int ikernel, double *pinput, ModulePW::PW_Basis *pw_rho, double *routput)
{
    std::complex<double> *recipOutput = new std::complex<double>[pw_rho->npw];

    pw_rho->real2recip(pinput, recipOutput);
    for (int ip = 0; ip < pw_rho->npw; ++ip)
    {
        recipOutput[ip] *= this->kernel[ikernel][ip];
    }
    pw_rho->recip2real(recipOutput, routput);

    delete[] recipOutput;
}

void ML_data::Laplacian(double * pinput, ModulePW::PW_Basis *pw_rho, double * routput)
{
    std::complex<double> *recipContainer = new std::complex<double>[pw_rho->npw];

    pw_rho->real2recip(pinput, recipContainer);
    for (int ip = 0; ip < pw_rho->npw; ++ip)
    {
        recipContainer[ip] *= - pw_rho->gg[ip] * pw_rho->tpiba2;
    }
    pw_rho->recip2real(recipContainer, routput);

    delete[] recipContainer;
}

void ML_data::divergence(double ** pinput, ModulePW::PW_Basis *pw_rho, double * routput)
{
    std::complex<double> *recipContainer = new std::complex<double>[pw_rho->npw];
    std::complex<double> img(0.0, 1.0);
    ModuleBase::GlobalFunc::ZEROS(routput, this->nx);
    for (int i = 0; i < 3; ++i)
    {
        pw_rho->real2recip(pinput[i], recipContainer);
        for (int ip = 0; ip < pw_rho->npw; ++ip)
        {
            recipContainer[ip] = img * pw_rho->gcar[ip][i] * pw_rho->tpiba * recipContainer[ip];
        }
        pw_rho->recip2real(recipContainer, routput, true);
    }

    delete[] recipContainer;
}

void ML_data::loadVector(std::string filename, std::vector<double> &data)
{
    std::vector<long unsigned int> cshape = {(long unsigned) this->nx};
    bool fortran_order = false;
    npy::LoadArrayFromNumpy(filename, cshape, fortran_order, data);
}

void ML_data::dumpVector(std::string filename, const std::vector<double> &data)
{
    const long unsigned cshape[] = {(long unsigned) this->nx}; // shape
    npy::SaveArrayAsNumpy(filename, false, 1, cshape, data);
}

void ML_data::tanh(std::vector<double> &pinput, std::vector<double> &routput, double chi)
{
    for (int i = 0; i < this->nx; ++i)
    {
        routput[i] = std::tanh(pinput[i] * chi);
    }
}

double ML_data::dtanh(double tanhx, double chi)
{
    return (1. - tanhx * tanhx) * chi;
}

std::string ML_data::file_name(std::string parameter, const int kernel_type, const double kernel_scaling)
{
    std::stringstream ss;
    ss << "./" << parameter << "_" << kernel_type << "_" << kernel_scaling << ".npy";
    return ss.str();
}
#endif