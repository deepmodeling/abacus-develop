#include "conv_coulomb_pot_k.h"

#include "../module_base/constants.h"
#include "../module_base/math_erf_complex.h"
#include "../module_basis/module_ao/ORB_atomic_lm.h"
#include "../module_hamilt_pw/hamilt_pwdft/global.h"
#include "module_parameter/parameter.h"
namespace Conv_Coulomb_Pot_K
{

std::vector<double> cal_psi_ccp(const std::vector<double>& psif, const double hybrid_alpha)
{
    std::vector<double> psik2_ccp(psif.size());
    for (size_t ik = 0; ik < psif.size(); ++ik)
        psik2_ccp[ik] = hybrid_alpha * ModuleBase::FOUR_PI * psif[ik];
    return psik2_ccp;
}

std::vector<double> cal_psi_ccp_cam(const std::vector<double>& psif,
                                    const std::vector<double>& k_radial,
                                    const double omega,
                                    const double hybrid_alpha,
                                    const double hybrid_beta)
{
    std::vector<double> psik2_ccp(psif.size());
    for (size_t ik = 0; ik < psif.size(); ++ik)
    {
        double fock_part = 1 - std::exp(-(k_radial[ik] * k_radial[ik]) / (4 * omega * omega));
        psik2_ccp[ik] = ModuleBase::FOUR_PI * psif[ik] * (hybrid_alpha + hybrid_beta * fock_part);
    }
    return psik2_ccp;
}

// rongshi add 2022-07-27
// Sphere truction -- Spencer
std::vector<double> cal_psi_hf(const std::vector<double>& psif,
                               const std::vector<double>& k_radial,
                               const int Rcut_type,
                               const double hybrid_alpha,
                               const double Rc)
{
    std::vector<double> psik2_ccp(psif.size());
    double Rw = 1.092 * Rc;
    double gamma = 5.0 / Rc;
    double dr = 0.01;
    int nr = static_cast<int>(5 * Rc / dr);
    for (size_t ik = 0; ik < psif.size(); ++ik)
    {
        if (Rcut_type == 0)
        {
            psik2_ccp[ik] = hybrid_alpha * ModuleBase::FOUR_PI * psif[ik] * (1 - std::cos(k_radial[ik] * Rc));
        }
        else if (Rcut_type == 1)
        {
            double sum_r = 0;
            for (size_t ir = 0; ir != nr; ++ir)
            {
                double r = dr * ir;
                sum_r += (std::erfc(gamma * r)
                          + 0.5 * std::erfc((std::log(r) - std::log(Rc)) / std::log(1.092)) * std::erf(gamma * r))
                         * std::sin(r * k_radial[ik]);
            }
            psik2_ccp[ik] = hybrid_alpha * ModuleBase::FOUR_PI * psif[ik] * dr * sum_r * k_radial[ik];
        }
    }
    return psik2_ccp;
}

std::vector<double> cal_psi_erfc(const std::vector<double>& psif,
                                 const std::vector<double>& k_radial,
                                 const double hse_omega,
                                 const double hybrid_beta)
{
    std::vector<double> psik2_ccp(psif.size());
    for (size_t ik = 0; ik < psif.size(); ++ik)
        psik2_ccp[ik] = hybrid_beta * ModuleBase::FOUR_PI * psif[ik]
                        * (1 - std::exp(-(k_radial[ik] * k_radial[ik]) / (4 * hse_omega * hse_omega)));
    return psik2_ccp;
}

std::vector<double> cal_psi_cam(const std::vector<double>& psif,
                                const std::vector<double>& k_radial,
                                const double omega,
                                const double hybrid_alpha,
                                const double hybrid_beta,
                                const int Rcut_type,
                                const double Rc)
{
    double eps = 1e-14;
    std::vector<double> psik2_ccp(psif.size());
    for (size_t ik = 0; ik < psif.size(); ++ik)
    {
        double coulomb_part = 1 - std::cos(k_radial[ik] * Rc);
        double temp0 = std::cos(k_radial[ik] * Rc) * std::erfc(omega * Rc);
        double temp1 = std::exp(-(k_radial[ik] * k_radial[ik]) / (4 * omega * omega));
        std::complex<double> temp2 = std::complex<double>(0, 0);
        std::complex<double> temp3 = std::complex<double>(0, 0);
        if (temp1 >= eps)
        {
            temp2 = ModuleBase::ErrorFunc::erf(0.5 * (ModuleBase::IMAG_UNIT * k_radial[ik] + 2 * omega * omega * Rc)
                                               / omega);
            temp3 = ModuleBase::NEG_IMAG_UNIT
                    * ModuleBase::ErrorFunc::erfi(0.5 * k_radial[ik] / omega + ModuleBase::IMAG_UNIT * omega * Rc);
        }
        std::complex<double> fock_part = -0.5 * (-2 + 2 * temp0 + temp1 * (temp2 + temp3));
        psik2_ccp[ik] = ModuleBase::FOUR_PI * psif[ik] * (hybrid_alpha * coulomb_part + hybrid_beta * fock_part.real());
    }
    return psik2_ccp;
}

template <>
Numerical_Orbital_Lm cal_orbs_ccp<Numerical_Orbital_Lm>(const Numerical_Orbital_Lm& orbs,
                                                        const Ccp_Type& ccp_type,
                                                        const std::map<std::string, double>& parameter,
                                                        const double rmesh_times)
{
    std::vector<double> psik2_ccp;
    switch (ccp_type)
    {
    case Ccp_Type::Ccp:
        psik2_ccp = cal_psi_ccp(orbs.get_psif(), parameter.at("hybrid_alpha"));
        break;
    case Ccp_Type::Hf:
        psik2_ccp = cal_psi_hf(orbs.get_psif(),
                               orbs.get_k_radial(),
                               parameter.at("Rcut_type"),
                               parameter.at("hybrid_alpha"),
                               parameter.at("hf_Rcut"));
        break;
    case Ccp_Type::Erfc:
        psik2_ccp = cal_psi_erfc(orbs.get_psif(),
                                 orbs.get_k_radial(),
                                 parameter.at("hse_omega"),
                                 parameter.at("hybrid_beta"));
        break;
    case Ccp_Type::Cam:
        psik2_ccp = cal_psi_cam(orbs.get_psif(),
                                orbs.get_k_radial(),
                                parameter.at("hse_omega"),
                                parameter.at("hybrid_alpha"),
                                parameter.at("hybrid_beta"),
                                parameter.at("Rcut_type"),
                                parameter.at("hf_Rcut"));
        break;
    case Ccp_Type::Ccp_Cam:
        psik2_ccp = cal_psi_ccp_cam(orbs.get_psif(),
                                    orbs.get_k_radial(),
                                    parameter.at("hse_omega"),
                                    parameter.at("hybrid_alpha"),
                                    parameter.at("hybrid_beta"));
        break;
    default:
        throw(ModuleBase::GlobalFunc::TO_STRING(__FILE__) + " line " + ModuleBase::GlobalFunc::TO_STRING(__LINE__));
        break;
    }

    const double dr = orbs.get_rab().back();
    const int Nr = (static_cast<int>(orbs.getNr() * rmesh_times)) | 1;
    std::vector<double> rab(Nr);
    for (size_t ir = 0; ir < std::min(orbs.getNr(), Nr); ++ir)
        rab[ir] = orbs.getRab(ir);
    for (size_t ir = orbs.getNr(); ir < Nr; ++ir)
        rab[ir] = dr;
    std::vector<double> r_radial(Nr);
    for (size_t ir = 0; ir < std::min(orbs.getNr(), Nr); ++ir)
        r_radial[ir] = orbs.getRadial(ir);
    for (size_t ir = orbs.getNr(); ir < Nr; ++ir)
        r_radial[ir] = orbs.get_r_radial().back() + (ir - orbs.getNr() + 1) * dr;

    Numerical_Orbital_Lm orbs_ccp;
    orbs_ccp.set_orbital_info(orbs.getLabel(),
                              orbs.getType(),
                              orbs.getL(),
                              orbs.getChi(),
                              Nr,
                              ModuleBase::GlobalFunc::VECTOR_TO_PTR(rab),
                              ModuleBase::GlobalFunc::VECTOR_TO_PTR(r_radial),
                              Numerical_Orbital_Lm::Psi_Type::Psik2,
                              ModuleBase::GlobalFunc::VECTOR_TO_PTR(psik2_ccp),
                              orbs.getNk(),
                              orbs.getDk(),
                              orbs.getDruniform(),
                              false,
                              true,
                              PARAM.inp.cal_force);
    return orbs_ccp;
}

template <>
double get_rmesh_proportion(const Numerical_Orbital_Lm& orbs, const double psi_threshold)
{
    for (int ir = orbs.getNr() - 1; ir >= 0; --ir)
    {
        if (std::abs(orbs.getPsi(ir)) >= psi_threshold)
            return static_cast<double>(ir) / orbs.getNr();
    }
    return 0.0;
}

} // namespace Conv_Coulomb_Pot_K
