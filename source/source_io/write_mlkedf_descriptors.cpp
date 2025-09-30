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
    const double* veff
)
{
    std::vector<std::vector<double>> nablaRho(3, std::vector<double>(this->cal_tool->nx, 0.));

    this->generate_descriptor(out_dir, pelec->charge->rho, pw_rho, nablaRho);

    std::vector<double> container(this->cal_tool->nx);
    std::vector<double> containernl(this->cal_tool->nx);

    const long unsigned cshape[] = {(long unsigned) this->cal_tool->nx}; // shape of container and containernl
    // enhancement factor of Pauli energy, and Pauli potential
    this->cal_tool->getF_KS(psi, pelec, pw_psi, pw_rho, ucell, nablaRho, container, containernl);

    Symmetry_rho srho;

    Charge* ptempRho = new Charge();
    ptempRho->nspin = PARAM.inp.nspin;
    ptempRho->nrxx = this->cal_tool->nx;
    ptempRho->rho_core = pelec->charge->rho_core;
    ptempRho->rho = new double*[1];
    ptempRho->rho[0] = new double[this->cal_tool->nx];
    ptempRho->rhog = new std::complex<double>*[1];
    ptempRho->rhog[0] = new std::complex<double>[pw_rho->npw];

    for (int ir = 0; ir < this->cal_tool->nx; ++ir){
        ptempRho->rho[0][ir] = container[ir];
    }
    srho.begin(0, *ptempRho, pw_rho, ucell.symm);
    for (int ir = 0; ir < this->cal_tool->nx; ++ir){
        container[ir] = ptempRho->rho[0][ir];
    }

    for (int ir = 0; ir < this->cal_tool->nx; ++ir){
        ptempRho->rho[0][ir] = containernl[ir];
    }
    srho.begin(0, *ptempRho, pw_rho, ucell.symm);
    for (int ir = 0; ir < this->cal_tool->nx; ++ir){
        containernl[ir] = ptempRho->rho[0][ir];
    }

    npy::SaveArrayAsNumpy(out_dir + "/enhancement.npy", false, 1, cshape, container);
    npy::SaveArrayAsNumpy(out_dir + "/pauli.npy", false, 1, cshape, containernl);

    for (int ir = 0; ir < this->cal_tool->nx; ++ir)
    {
        container[ir] = veff[ir];
    }
    npy::SaveArrayAsNumpy(out_dir + "/veff.npy", false, 1, cshape, container);

    delete ptempRho;
}

void Write_MLKEDF_Descriptors::generateTrainData_KS(
    const std::string& out_dir,
    psi::Psi<std::complex<float>> *psi,
    elecstate::ElecState *pelec,
    ModulePW::PW_Basis_K *pw_psi,
    ModulePW::PW_Basis *pw_rho,
    UnitCell& ucell,
    const double* veff
)
{
    psi::Psi<std::complex<double>, base_device::DEVICE_CPU> psi_double(*psi);

    this->generateTrainData_KS(out_dir, &psi_double, pelec, pw_psi, pw_rho, ucell, veff);
}

#if ((defined __CUDA) || (defined __ROCM))
void Write_MLKEDF_Descriptors::generateTrainData_KS(
    const std::string& out_dir,
    psi::Psi<std::complex<double>, base_device::DEVICE_GPU>* psi,
    elecstate::ElecState *pelec,
    ModulePW::PW_Basis_K *pw_psi,
    ModulePW::PW_Basis *pw_rho,
    UnitCell& ucell,
    const double* veff
)
{
    psi::Psi<std::complex<double>, base_device::DEVICE_CPU> psi_cpu(*psi);

    this->generateTrainData_KS(out_dir, &psi_cpu, pelec, pw_psi, pw_rho, ucell, veff);
}

void Write_MLKEDF_Descriptors::generateTrainData_KS(
    const std::string& dir,
    psi::Psi<std::complex<float>, base_device::DEVICE_GPU>* psi,
    elecstate::ElecState *pelec,
    ModulePW::PW_Basis_K *pw_psi,
    ModulePW::PW_Basis *pw_rho,
    UnitCell& ucell,
    const double *veff
)
{
    psi::Psi<std::complex<double>, base_device::DEVICE_CPU> psi_cpu_double(*psi);

    this->generateTrainData_KS(dir, &psi_cpu_double, pelec, pw_psi, pw_rho, ucell, veff);
}
#endif

void Write_MLKEDF_Descriptors::generate_descriptor(
    const std::string& out_dir,
    const double * const *prho, 
    ModulePW::PW_Basis *pw_rho,
    std::vector<std::vector<double>> &nablaRho
)
{
    // container which will contain gamma, p, q in turn
    std::vector<double> container(this->cal_tool->nx);
    std::vector<double> new_container(this->cal_tool->nx);
    // container contains gammanl, pnl, qnl in turn
    std::vector<double> containernl(this->cal_tool->nx);
    std::vector<double> new_containernl(this->cal_tool->nx);

    const long unsigned cshape[] = {(long unsigned) this->cal_tool->nx}; // shape of container and containernl

    // rho
    std::vector<double> rho(this->cal_tool->nx);
    for (int ir = 0; ir < this->cal_tool->nx; ++ir){
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
    npy::SaveArrayAsNumpy("p.npy", false, 1, cshape, container);

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