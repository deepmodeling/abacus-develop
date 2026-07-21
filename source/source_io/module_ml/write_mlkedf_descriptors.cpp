#ifdef __MLALGO

#include "write_mlkedf_descriptors.h"

#include "npy.hpp"
#include "source_estate/module_charge/symmetry_rho.h"

namespace ModuleIO
{
void Write_MLKEDF_Descriptors::generateTrainData_KS(
    const std::string& out_dir,
    psi::Psi<std::complex<double>> *psi,
    elecstate::ElecState *pelec,
    ModulePW::PW_Basis_K *pw_psi,
    ModulePW::PW_Basis *pw_rho,
    UnitCell& ucell,
    const double* veff,
    const int nrxx
)
{
    if (nrxx <= 0)
    {
        ModuleBase::WARNING_QUIT("Write_MLKEDF_Descriptors::generateTrainData_KS", "nrxx must be greater than 0");
    }

    std::vector<std::vector<double>> drho(3, std::vector<double>(nrxx, 0.));

    this->generate_descriptor(out_dir, pelec->charge->rho, pw_rho, drho, nrxx);

    std::vector<double> enhancement(nrxx);
    std::vector<double> pauli(nrxx);

    this->cal_tool->getF_KS(psi, pelec, pw_psi, pw_rho, ucell, drho, enhancement, pauli);

    Symmetry_rho srho;

    std::vector<double> rho_vec(nrxx);
    std::vector<std::complex<double>> rhog_vec(pw_rho->npw);
    double* rho_ptr = rho_vec.data();
    std::complex<double>* rhog_ptr = rhog_vec.data();

    std::copy(enhancement.begin(), enhancement.end(), rho_vec.begin());
    srho.begin(0, &rho_ptr, &rhog_ptr, pw_rho->npw, nullptr, pw_rho, ucell.symm);
    std::copy(rho_vec.begin(), rho_vec.end(), enhancement.begin());

    std::copy(pauli.begin(), pauli.end(), rho_vec.begin());
    srho.begin(0, &rho_ptr, &rhog_ptr, pw_rho->npw, nullptr, pw_rho, ucell.symm);
    std::copy(rho_vec.begin(), rho_vec.end(), pauli.begin());


    // output data in .npy format
    const long unsigned cshape[] = {(long unsigned) nrxx};
    npy::SaveArrayAsNumpy(out_dir + "/enhancement.npy", false, 1, cshape, enhancement);
    npy::SaveArrayAsNumpy(out_dir + "/pauli.npy", false, 1, cshape, pauli);

    for (int ir = 0; ir < nrxx; ++ir)
    {
        enhancement[ir] = veff[ir];
    }
    npy::SaveArrayAsNumpy(out_dir + "/veff.npy", false, 1, cshape, enhancement);
}

void Write_MLKEDF_Descriptors::generateTrainData_KS(
    const std::string& out_dir,
    psi::Psi<std::complex<float>> *psi,
    elecstate::ElecState *pelec,
    ModulePW::PW_Basis_K *pw_psi,
    ModulePW::PW_Basis *pw_rho,
    UnitCell& ucell,
    const double* veff,
    const int nrxx
)
{
    psi::Psi<std::complex<double>, base_device::DEVICE_CPU> psi_double(*psi);

    this->generateTrainData_KS(out_dir, &psi_double, pelec, pw_psi, pw_rho, ucell, veff, nrxx);
}

#if ((defined __CUDA) || (defined __ROCM))
void Write_MLKEDF_Descriptors::generateTrainData_KS(
    const std::string& out_dir,
    psi::Psi<std::complex<double>, base_device::DEVICE_GPU>* psi,
    elecstate::ElecState *pelec,
    ModulePW::PW_Basis_K *pw_psi,
    ModulePW::PW_Basis *pw_rho,
    UnitCell& ucell,
    const double* veff,
    const int nrxx
)
{
    psi::Psi<std::complex<double>, base_device::DEVICE_CPU> psi_cpu(*psi);

    this->generateTrainData_KS(out_dir, &psi_cpu, pelec, pw_psi, pw_rho, ucell, veff, nrxx);
}

void Write_MLKEDF_Descriptors::generateTrainData_KS(
    const std::string& dir,
    psi::Psi<std::complex<float>, base_device::DEVICE_GPU>* psi,
    elecstate::ElecState *pelec,
    ModulePW::PW_Basis_K *pw_psi,
    ModulePW::PW_Basis *pw_rho,
    UnitCell& ucell,
    const double *veff,
    const int nrxx
)
{
    psi::Psi<std::complex<double>, base_device::DEVICE_CPU> psi_cpu_double(*psi);

    this->generateTrainData_KS(dir, &psi_cpu_double, pelec, pw_psi, pw_rho, ucell, veff, nrxx);
}
#endif

void Write_MLKEDF_Descriptors::generate_descriptor(
    const std::string& out_dir,
    const double * const *prho, 
    ModulePW::PW_Basis *pw_rho,
    std::vector<std::vector<double>> &nablaRho,
    const int nrxx
)
{
    // container which will contain gamma, p, q in turn
    std::vector<double> container(nrxx);
    std::vector<double> new_container(nrxx);
    // container contains gammanl, pnl, qnl in turn
    std::vector<double> containernl(nrxx);
    std::vector<double> new_containernl(nrxx);

    const long unsigned cshape[] = {(long unsigned) nrxx};

    // rho
    std::vector<double> rho(nrxx);
    for (int ir = 0; ir < nrxx; ++ir){
        rho[ir] = prho[0][ir];
    }
    npy::SaveArrayAsNumpy(out_dir + "/rho.npy", false, 1, cshape, rho);

    // gamma
    this->cal_tool->getGamma(prho, container);
    npy::SaveArrayAsNumpy(out_dir + "/gamma.npy", false, 1, cshape, container);

    for (int ik = 0; ik < this->cal_tool->nkernel; ++ik)
    {
        int ktype = this->cal_tool->kernel_type[ik];
        double kscaling = this->cal_tool->kernel_scaling[ik];

        // gamma_nl
        this->cal_tool->getGammanl(ik, container, pw_rho, containernl);
        npy::SaveArrayAsNumpy(this->file_name(out_dir, "gammanl", ktype, kscaling), false, 1, cshape, containernl);

        // xi = gamma_nl/gamma
        this->cal_tool->getXi(container, containernl, new_container);
        npy::SaveArrayAsNumpy(this->file_name(out_dir, "xi", ktype, kscaling), false, 1, cshape, new_container);

        // tanhxi = tanh(xi)
        this->cal_tool->getTanhXi(ik, container, containernl, new_container);
        npy::SaveArrayAsNumpy(this->file_name(out_dir, "tanhxi", ktype, kscaling), false, 1, cshape, new_container);

        // (tanhxi)_nl
        this->cal_tool->getTanhXi_nl(ik, new_container, pw_rho, new_containernl);
        npy::SaveArrayAsNumpy(this->file_name(out_dir, "tanhxi_nl", ktype, kscaling), false, 1, cshape, new_containernl);
    }

    // nabla rho
    this->cal_tool->getNablaRho(prho, pw_rho, nablaRho);
    npy::SaveArrayAsNumpy(out_dir + "/nablaRhox.npy", false, 1, cshape, nablaRho[0]);
    npy::SaveArrayAsNumpy(out_dir + "/nablaRhoy.npy", false, 1, cshape, nablaRho[1]);
    npy::SaveArrayAsNumpy(out_dir + "/nablaRhoz.npy", false, 1, cshape, nablaRho[2]);

    // p
    this->cal_tool->getP(prho, pw_rho, nablaRho, container);
    npy::SaveArrayAsNumpy(out_dir + "/p.npy", false, 1, cshape, container);

    for (int ik = 0; ik < this->cal_tool->nkernel; ++ik)
    {
        int ktype = this->cal_tool->kernel_type[ik];
        double kscaling = this->cal_tool->kernel_scaling[ik];

        // p_nl
        this->cal_tool->getPnl(ik, container, pw_rho, containernl);
        npy::SaveArrayAsNumpy(this->file_name(out_dir, "pnl", ktype, kscaling), false, 1, cshape, containernl);

        // tanh(p_nl)
        this->cal_tool->getTanh_Pnl(ik, containernl, new_containernl);
        npy::SaveArrayAsNumpy(this->file_name(out_dir, "tanh_pnl", ktype, kscaling), false, 1, cshape, new_containernl);

        // tanh(p)
        this->cal_tool->getTanhP(container, new_container);
        npy::SaveArrayAsNumpy(out_dir + "/tanhp.npy", false, 1, cshape, new_container);

        // tanh(p)_nl
        this->cal_tool->getTanhP_nl(ik, new_container, pw_rho, new_containernl);
        npy::SaveArrayAsNumpy(this->file_name(out_dir, "tanhp_nl", ktype, kscaling), false, 1, cshape, new_containernl);
    }

    // q
    this->cal_tool->getQ(prho, pw_rho, container);
    npy::SaveArrayAsNumpy(out_dir + "/q.npy", false, 1, cshape, container);

    for (int ik = 0; ik < this->cal_tool->nkernel; ++ik)
    {
        int ktype = this->cal_tool->kernel_type[ik];
        double kscaling = this->cal_tool->kernel_scaling[ik];

        // q_nl
        this->cal_tool->getQnl(ik, container, pw_rho, containernl);
        npy::SaveArrayAsNumpy(this->file_name(out_dir, "qnl", ktype, kscaling), false, 1, cshape, containernl);

        // tanh(q_nl)
        this->cal_tool->getTanh_Qnl(ik, containernl, new_containernl);
        npy::SaveArrayAsNumpy(this->file_name(out_dir, "tanh_qnl", ktype, kscaling), false, 1, cshape, new_containernl);

        // tanh(q)
        this->cal_tool->getTanhQ(container, new_container);
        npy::SaveArrayAsNumpy(out_dir + "/tanhq.npy", false, 1, cshape, new_container);

        // tanh(q)_nl
        this->cal_tool->getTanhQ_nl(ik, new_container, pw_rho, new_containernl);
        npy::SaveArrayAsNumpy(this->file_name(out_dir, "tanhq_nl", ktype, kscaling), false, 1, cshape, new_containernl);
    }
}

std::string Write_MLKEDF_Descriptors::file_name(
    const std::string& out_dir,
    std::string parameter,
    const int kernel_type,
    const double kernel_scaling
)
{
    std::stringstream ss;
    ss << out_dir << "/" << parameter << "_" << kernel_type << "_" << kernel_scaling << ".npy";
    return ss.str();
}

}

#endif
