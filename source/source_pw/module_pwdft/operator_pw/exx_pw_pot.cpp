#include "op_exx_pw.h"
#include "source_io/module_parameter/parameter.h"
#include "source_pw/module_pwdft/global.h"

namespace hamilt
{
template <typename Real, typename Device>
void get_exx_potential(const K_Vectors* kv,
                       const ModulePW::PW_Basis_K* wfcpw,
                       ModulePW::PW_Basis* rhopw_dev,
                       Real* pot,
                       double tpiba,
                       bool gamma_extrapolation,
                       double ucell_omega,
                       int ik,
                       int iq,
                       bool is_stress)
{
    using setmem_real_cpu_op = base_device::memory::set_memory_op<Real, base_device::DEVICE_CPU>;
    using syncmem_real_c2d_op = base_device::memory::synchronize_memory_op<Real, base_device::DEVICE_CPU, Device>;

    Real nqs_half1 = 0.5 * kv->nmp[0];
    Real nqs_half2 = 0.5 * kv->nmp[1];
    Real nqs_half3 = 0.5 * kv->nmp[2];

    Real* pot_cpu = nullptr;
    int nks = wfcpw->nks, npw = rhopw_dev->npw;
    double tpiba2 = tpiba * tpiba;
    pot_cpu = new Real[npw];
    // fill zero
    setmem_real_cpu_op()(pot_cpu, 0, npw);

    // calculate Fock pot
    auto param_fock = GlobalC::exx_info.info_global.coulomb_param[Conv_Coulomb_Pot_K::Coulomb_Type::Fock];
    for (int i = 0; i < param_fock.size(); i++)
    {
        auto param = param_fock[i];
        double exx_div = OperatorEXXPW<std::complex<Real>, Device>::fock_div[i];
        double alpha = std::stod(param["alpha"]);
        const ModuleBase::Vector3<double> k_c = wfcpw->kvec_c[ik];
        const ModuleBase::Vector3<double> k_d = wfcpw->kvec_d[ik];
        const ModuleBase::Vector3<double> q_c = wfcpw->kvec_c[iq];
        const ModuleBase::Vector3<double> q_d = wfcpw->kvec_d[iq];

#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
        for (int ig = 0; ig < rhopw_dev->npw; ig++)
        {
            const ModuleBase::Vector3<double> g_d = rhopw_dev->gdirect[ig];
            const ModuleBase::Vector3<double> kqg_d = k_d - q_d + g_d;
            // For gamma_extrapolation (https://doi.org/10.1103/PhysRevB.79.205114)
            // 7/8 of the points in the grid are "activated" and 1/8 are disabled.
            // grid_factor is designed for the 7/8 of the grid to function like all of the points
            Real grid_factor = 1;
            double extrapolate_grid = 8.0 / 7.0;
            if (gamma_extrapolation)
            {
                // if isint(kqg_d[0] * nqs_half1) && isint(kqg_d[1] * nqs_half2) && isint(kqg_d[2] * nqs_half3)
                auto isint = [](double x) {
                    double epsilon = 1e-6; // this follows the isint judgement in q-e
                    return std::abs(x - std::round(x)) < epsilon;
                };
                if (isint(kqg_d[0] * nqs_half1) && isint(kqg_d[1] * nqs_half2) && isint(kqg_d[2] * nqs_half3))
                {
                    grid_factor = 0;
                }
                else
                {
                    grid_factor = extrapolate_grid;
                }
            }

            const int nk_fac = PARAM.inp.nspin == 2 ? 2 : 1;
            const int nk = nks / nk_fac;

            Real gg = (k_c - q_c + rhopw_dev->gcar[ig]).norm2() * tpiba2;
            // if (kqgcar2 > 1e-12) // vasp uses 1/40 of the smallest (k spacing)**2
            if (gg >= 1e-8)
            {
                Real fac = -ModuleBase::FOUR_PI * ModuleBase::e2 / gg;
                pot_cpu[ig] += fac * grid_factor * alpha;
            }
            // }
            else
            {
                pot_cpu[ig] += exx_div * alpha;
            }
            // assert(is_finite(density_recip[ig]));
        }
    }

    // calculate erfc pot
    auto param_erfc = GlobalC::exx_info.info_global.coulomb_param[Conv_Coulomb_Pot_K::Coulomb_Type::Erfc];
    for (int i = 0; i < param_erfc.size(); i++)
    {
        auto param = param_erfc[i];
        double erfc_omega = std::stod(param["omega"]);
        double erfc_omega2 = erfc_omega * erfc_omega;
        double alpha = std::stod(param["alpha"]);
        // double exx_div = OperatorEXXPW<std::complex<Real>, Device>::erfc_div[i];
        double exx_div = exx_divergence(Conv_Coulomb_Pot_K::Coulomb_Type::Erfc,
                                          erfc_omega,
                                          kv,
                                          wfcpw,
                                          rhopw_dev,
                                          tpiba,
                                          gamma_extrapolation,
                                          ucell_omega);
        const ModuleBase::Vector3<double> k_c = wfcpw->kvec_c[ik];
        const ModuleBase::Vector3<double> k_d = wfcpw->kvec_d[ik];
        const ModuleBase::Vector3<double> q_c = wfcpw->kvec_c[iq];
        const ModuleBase::Vector3<double> q_d = wfcpw->kvec_d[iq];

#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
        for (int ig = 0; ig < rhopw_dev->npw; ig++)
        {
            const ModuleBase::Vector3<double> g_d = rhopw_dev->gdirect[ig];
            const ModuleBase::Vector3<double> kqg_d = k_d - q_d + g_d;
            // For gamma_extrapolation (https://doi.org/10.1103/PhysRevB.79.205114)
            // 7/8 of the points in the grid are "activated" and 1/8 are disabled.
            // grid_factor is designed for the 7/8 of the grid to function like all of the points
            Real grid_factor = 1;
            double extrapolate_grid = 8.0 / 7.0;
            if (gamma_extrapolation)
            {
                // if isint(kqg_d[0] * nqs_half1) && isint(kqg_d[1] * nqs_half2) && isint(kqg_d[2] * nqs_half3)
                auto isint = [](double x) {
                    double epsilon = 1e-6; // this follows the isint judgement in q-e
                    return std::abs(x - std::round(x)) < epsilon;
                };
                if (isint(kqg_d[0] * nqs_half1) && isint(kqg_d[1] * nqs_half2) && isint(kqg_d[2] * nqs_half3))
                {
                    grid_factor = 0;
                }
                else
                {
                    grid_factor = extrapolate_grid;
                }
            }

            const int nk_fac = PARAM.inp.nspin == 2 ? 2 : 1;
            const int nk = nks / nk_fac;
            // const int ig_kq = ik * nks * npw + iq * npw + ig;

            Real gg = (k_c - q_c + rhopw_dev->gcar[ig]).norm2() * tpiba2;
            // if (kqgcar2 > 1e-12) // vasp uses 1/40 of the smallest (k spacing)**2
            if (gg >= 1e-8)
            {
                Real fac = -ModuleBase::FOUR_PI * ModuleBase::e2 / gg;
                pot_cpu[ig] += fac * (1.0 - std::exp(-gg / 4.0 / erfc_omega2)) * grid_factor * alpha;
            }
            // }
            else
            {
                // if (PARAM.inp.dft_functional == "hse")
                if (!gamma_extrapolation)
                {
                    if (is_stress)
                        pot_cpu[ig] += (- ModuleBase::PI * ModuleBase::e2 / erfc_omega2) * alpha;
                    else
                        pot_cpu[ig] += (exx_div - ModuleBase::PI * ModuleBase::e2 / erfc_omega2) * alpha;
                }
                else
                {
                    pot_cpu[ig] += exx_div * alpha;
                }
            }
            // assert(is_finite(density_recip[ig]));
        }
    }

    // copy the potential to the device memory
#ifdef __CUDA
    cudaError_t err = cudaHostRegister(pot_cpu, sizeof(Real) * npw, cudaHostRegisterPortable);
    if (err != cudaSuccess) {
        throw std::runtime_error("failed to register potential CPU memory operations");
    }
#endif
    syncmem_real_c2d_op()(pot, pot_cpu, rhopw_dev->npw);
#ifdef __CUDA
    cudaHostUnregister(pot_cpu);
#endif

    delete pot_cpu;
}

template <typename Real, typename Device>
void get_exx_stress_potential(const K_Vectors* kv,
                              const ModulePW::PW_Basis_K* wfcpw,
                              ModulePW::PW_Basis* rhopw_dev,
                              Real* pot,
                              double tpiba,
                              bool gamma_extrapolation,
                              double ucell_omega,
                              int ik,
                              int iq)
{
    using setmem_real_cpu_op = base_device::memory::set_memory_op<Real, base_device::DEVICE_CPU>;
    using syncmem_real_c2d_op = base_device::memory::synchronize_memory_op<Real, base_device::DEVICE_CPU, Device>;

    Real nqs_half1 = 0.5 * kv->nmp[0];
    Real nqs_half2 = 0.5 * kv->nmp[1];
    Real nqs_half3 = 0.5 * kv->nmp[2];

    Real* pot_cpu = nullptr;
    int nks = wfcpw->nks, npw = rhopw_dev->npw;
    double tpiba2 = tpiba * tpiba;
    pot_cpu = new Real[npw];
    // fill zero
    setmem_real_cpu_op()(pot_cpu, 0, npw);

    // calculate Fock pot
    auto param_fock = GlobalC::exx_info.info_global.coulomb_param[Conv_Coulomb_Pot_K::Coulomb_Type::Fock];
    for (auto param: param_fock)
    {
        // double exx_div = exx_divergence(Conv_Coulomb_Pot_K::Coulomb_Type::Fock,
        //                                 0.0,
        //                                 kv,
        //                                 wfcpw,
        //                                 rhopw_dev,
        //                                 tpiba,
        //                                 gamma_extrapolation,
        //                                 ucell_omega);
        double alpha = std::stod(param["alpha"]);

        const ModuleBase::Vector3<double> k_c = wfcpw->kvec_c[ik];
        const ModuleBase::Vector3<double> k_d = wfcpw->kvec_d[ik];
        const ModuleBase::Vector3<double> q_c = wfcpw->kvec_c[iq];
        const ModuleBase::Vector3<double> q_d = wfcpw->kvec_d[iq];

#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
        for (int ig = 0; ig < rhopw_dev->npw; ig++)
        {
            const ModuleBase::Vector3<double> g_d = rhopw_dev->gdirect[ig];
            const ModuleBase::Vector3<double> kqg_d = k_d - q_d + g_d;
            // For gamma_extrapolation (https://doi.org/10.1103/PhysRevB.79.205114)
            // 7/8 of the points in the grid are "activated" and 1/8 are disabled.
            // grid_factor is designed for the 7/8 of the grid to function like all of the points
            Real grid_factor = 1;
            double extrapolate_grid = 8.0 / 7.0;
            if (gamma_extrapolation)
            {
                // if isint(kqg_d[0] * nqs_half1) && isint(kqg_d[1] * nqs_half2) && isint(kqg_d[2] * nqs_half3)
                auto isint = [](double x) {
                    double epsilon = 1e-6; // this follows the isint judgement in q-e
                    return std::abs(x - std::round(x)) < epsilon;
                };
                if (isint(kqg_d[0] * nqs_half1) && isint(kqg_d[1] * nqs_half2) && isint(kqg_d[2] * nqs_half3))
                {
                    grid_factor = 0;
                }
                else
                {
                    grid_factor = extrapolate_grid;
                }
            }

            const int nk_fac = PARAM.inp.nspin == 2 ? 2 : 1;
            const int nk = nks / nk_fac;
            // const int ig_kq = ik * nks * npw + iq * npw + ig;

            Real gg = (k_c - q_c + rhopw_dev->gcar[ig]).norm2() * tpiba2;
            // if (kqgcar2 > 1e-12) // vasp uses 1/40 of the smallest (k spacing)**2
            if (gg >= 1e-8)
            {
                Real fac = -ModuleBase::FOUR_PI * ModuleBase::e2 / gg;
                pot_cpu[ig] += 1.0 / gg * grid_factor * alpha;
            }
        }
    }

    // calculate erfc pot
    auto param_erfc = GlobalC::exx_info.info_global.coulomb_param[Conv_Coulomb_Pot_K::Coulomb_Type::Erfc];
    for (auto param: param_erfc)
    {
        double erfc_omega = std::stod(param["omega"]);
        double erfc_omega2 = erfc_omega * erfc_omega;
        double alpha = std::stod(param["alpha"]);
        // double exx_div = exx_divergence(Conv_Coulomb_Pot_K::Coulomb_Type::Erfc,
        //                                 erfc_omega,
        //                                 kv,
        //                                 wfcpw,
        //                                 rhopw_dev,
        //                                 tpiba,
        //                                 gamma_extrapolation,
        //                                 ucell_omega);

        const ModuleBase::Vector3<double> k_c = wfcpw->kvec_c[ik];
        const ModuleBase::Vector3<double> k_d = wfcpw->kvec_d[ik];
        const ModuleBase::Vector3<double> q_c = wfcpw->kvec_c[iq];
        const ModuleBase::Vector3<double> q_d = wfcpw->kvec_d[iq];

#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
        for (int ig = 0; ig < rhopw_dev->npw; ig++)
        {
            const ModuleBase::Vector3<double> g_d = rhopw_dev->gdirect[ig];
            const ModuleBase::Vector3<double> kqg_d = k_d - q_d + g_d;
            // For gamma_extrapolation (https://doi.org/10.1103/PhysRevB.79.205114)
            // 7/8 of the points in the grid are "activated" and 1/8 are disabled.
            // grid_factor is designed for the 7/8 of the grid to function like all of the points
            Real grid_factor = 1;
            double extrapolate_grid = 8.0 / 7.0;
            if (gamma_extrapolation)
            {
                // if isint(kqg_d[0] * nqs_half1) && isint(kqg_d[1] * nqs_half2) && isint(kqg_d[2] * nqs_half3)
                auto isint = [](double x) {
                    double epsilon = 1e-6; // this follows the isint judgement in q-e
                    return std::abs(x - std::round(x)) < epsilon;
                };
                if (isint(kqg_d[0] * nqs_half1) && isint(kqg_d[1] * nqs_half2) && isint(kqg_d[2] * nqs_half3))
                {
                    grid_factor = 0;
                }
                else
                {
                    grid_factor = extrapolate_grid;
                }
            }

            const int nk_fac = PARAM.inp.nspin == 2 ? 2 : 1;
            const int nk = nks / nk_fac;
            // const int ig_kq = ik * nks * npw + iq * npw + ig;

            Real gg = (k_c - q_c + rhopw_dev->gcar[ig]).norm2() * tpiba2;
            // if (kqgcar2 > 1e-12) // vasp uses 1/40 of the smallest (k spacing)**2
            if (gg >= 1e-8)
            {
                Real fac = -ModuleBase::FOUR_PI * ModuleBase::e2 / gg;
                pot_cpu[ig] += (1.0 - (1.0 + gg / 4.0 / erfc_omega2) * std::exp(-gg / 4.0 / erfc_omega2))
                               / (1.0 - std::exp(-gg / 4.0 / erfc_omega2)) / gg * grid_factor * alpha;
            }
            // }
            else
            {
                // if (PARAM.inp.dft_functional == "hse")
                if (!gamma_extrapolation)
                {
                    pot_cpu[ig] += 1.0 / 4.0 / erfc_omega2 * alpha;
                }
            }
            // assert(is_finite(density_recip[ig]));
        }
    }

    // copy the potential to the device memory
#ifdef __CUDA
    cudaError_t err = cudaHostRegister(pot_cpu, sizeof(Real) * npw, cudaHostRegisterPortable);
    if (err != cudaSuccess) {
        throw std::runtime_error("failed to register potential CPU memory operations");
    }
#endif
    syncmem_real_c2d_op()(pot, pot_cpu, rhopw_dev->npw);
#ifdef __CUDA
    cudaHostUnregister(pot_cpu);
#endif

    delete pot_cpu;
}

double exx_divergence(Conv_Coulomb_Pot_K::Coulomb_Type coulomb_type,
                      double erfc_omega,
                      const K_Vectors* kv,
                      const ModulePW::PW_Basis_K* wfcpw,
                      ModulePW::PW_Basis* rhopw_dev,
                      double tpiba,
                      bool gamma_extrapolation,
                      double ucell_omega)
{
    double exx_div = 0;

    double nqs_half1 = 0.5 * kv->nmp[0];
    double nqs_half2 = 0.5 * kv->nmp[1];
    double nqs_half3 = 0.5 * kv->nmp[2];

    int nk_fac = PARAM.inp.nspin == 2 ? 2 : 1;

    // here we follow the exx_divergence subroutine in q-e (PW/src/exx_base.f90)
    double alpha = 10.0 / wfcpw->gk_ecut;
    double tpiba2 = tpiba * tpiba;
    double div = 0;

    // this is the \sum_q F(q) part
    // temporarily for all k points, should be replaced to q points later
    for (int ik = 0; ik < wfcpw->nks; ik++)
    {
        const ModuleBase::Vector3<double> k_c = wfcpw->kvec_c[ik];
        const ModuleBase::Vector3<double> k_d = wfcpw->kvec_d[ik];
#ifdef _OPENMP
#pragma omp parallel for reduction(+ : div)
#endif
        for (int ig = 0; ig < rhopw_dev->npw; ig++)
        {
            const ModuleBase::Vector3<double> q_c = k_c + rhopw_dev->gcar[ig];
            const ModuleBase::Vector3<double> q_d = k_d + rhopw_dev->gdirect[ig];
            double qq = q_c.norm2();
            // For gamma_extrapolation (https://doi.org/10.1103/PhysRevB.79.205114)
            // 7/8 of the points in the grid are "activated" and 1/8 are disabled.
            // grid_factor is designed for the 7/8 of the grid to function like all of the points
            double grid_factor = 1;
            double extrapolate_grid = 8.0 / 7.0;
            if (gamma_extrapolation)
            {
                auto isint = [](double x) {
                    double epsilon = 1e-6; // this follows the isint judgement in q-e
                    return std::abs(x - std::round(x)) < epsilon;
                };
                if (isint(q_d[0] * nqs_half1) && isint(q_d[1] * nqs_half2) && isint(q_d[2] * nqs_half3))
                {
                    grid_factor = 0;
                }
                else
                {
                    grid_factor = extrapolate_grid;
                }
            }

            if (qq <= 1e-8)
                continue;
            // else if (PARAM.inp.dft_functional == "hse")
            else if (coulomb_type == Conv_Coulomb_Pot_K::Coulomb_Type::Erfc)
            {
                double omega = erfc_omega;
                double omega2 = omega * omega;
                div += std::exp(-alpha * qq) / qq * (1.0 - std::exp(-qq * tpiba2 / 4.0 / omega2)) * grid_factor;
            }
            else
            {
                div += std::exp(-alpha * qq) / qq * grid_factor;
            }
        }
    }

    Parallel_Reduce::reduce_pool(div);
    // std::cout << "EXX div: " << div << std::endl;

    // if (PARAM.inp.dft_functional == "hse")
    if (!gamma_extrapolation)
    {
        if (coulomb_type == Conv_Coulomb_Pot_K::Coulomb_Type::Erfc)
        {
            double omega = erfc_omega;
            div += tpiba2 / 4.0 / omega / omega; // compensate for the finite value when qq = 0
        }
        else
        {
            div -= alpha;
        }
    }

    div *= ModuleBase::e2 * ModuleBase::FOUR_PI / tpiba2 / wfcpw->nks;
    //    std::cout << "div: " << div << std::endl;

    // numerically value the mean value of F(q) in the reciprocal space
    // This means we need to calculate the average of F(q) in the first brillouin zone
    alpha /= tpiba2;
    int nqq = 100000;
    double dq = 5.0 / std::sqrt(alpha) / nqq;
    double aa = 0.0;
    // if (PARAM.inp.dft_functional == "hse")
    if (coulomb_type == Conv_Coulomb_Pot_K::Coulomb_Type::Erfc)
    {
        double omega = erfc_omega;
        double omega2 = omega * omega;
#ifdef _OPENMP
#pragma omp parallel for reduction(+ : aa)
#endif
        for (int i = 0; i < nqq; i++)
        {
            double q = dq * (i + 0.5);
            aa -= exp(-alpha * q * q) * exp(-q * q / 4.0 / omega2) * dq;
        }
    }
    aa *= 8 / ModuleBase::FOUR_PI;
    aa += 1.0 / std::sqrt(alpha * ModuleBase::PI);

    div -= ModuleBase::e2 * ucell_omega * aa;
    exx_div = div * wfcpw->nks / nk_fac;
    //    exx_div = 0;
    //    std::cout << "EXX divergence: " << exx_div << std::endl;

    return exx_div;
}
template class OperatorEXXPW<std::complex<float>, base_device::DEVICE_CPU>;
template class OperatorEXXPW<std::complex<double>, base_device::DEVICE_CPU>;
template void get_exx_potential<float, base_device::DEVICE_CPU>(const K_Vectors*,
                                                                const ModulePW::PW_Basis_K*,
                                                                ModulePW::PW_Basis*,
                                                                float*,
                                                                double,
                                                                bool,
                                                                double,
                                                                int,
                                                                int,
                                                                        bool);
template void get_exx_potential<double, base_device::DEVICE_CPU>(const K_Vectors*,
                                                                 const ModulePW::PW_Basis_K*,
                                                                 ModulePW::PW_Basis*,
                                                                 double*,
                                                                 double,
                                                                 bool,
                                                                 double,
                                                                 int,
                                                                 int,
                                                                        bool);
template void get_exx_stress_potential<float, base_device::DEVICE_CPU>(const K_Vectors*,
                                                                       const ModulePW::PW_Basis_K*,
                                                                       ModulePW::PW_Basis*,
                                                                       float*,
                                                                       double,
                                                                       bool,
                                                                       double,
                                                                       int,
                                                                       int);
template void get_exx_stress_potential<double, base_device::DEVICE_CPU>(const K_Vectors*,
                                                                        const ModulePW::PW_Basis_K*,
                                                                        ModulePW::PW_Basis*,
                                                                        double*,
                                                                        double,
                                                                        bool,
                                                                        double,
                                                                        int,
                                                                        int);
#if ((defined __CUDA) || (defined __ROCM))
template class OperatorEXXPW<std::complex<float>, base_device::DEVICE_GPU>;
template class OperatorEXXPW<std::complex<double>, base_device::DEVICE_GPU>;
template void get_exx_potential<float, base_device::DEVICE_GPU>(const K_Vectors*,
                                                                const ModulePW::PW_Basis_K*,
                                                                ModulePW::PW_Basis*,
                                                                float*,
                                                                double,
                                                                bool,
                                                                double,
                                                                int,
                                                                int,
                                                                bool);
template void get_exx_potential<double, base_device::DEVICE_GPU>(const K_Vectors*,
                                                                 const ModulePW::PW_Basis_K*,
                                                                 ModulePW::PW_Basis*,
                                                                 double*,
                                                                 double,
                                                                 bool,
                                                                 double,
                                                                 int,
                                                                 int,
                                                                 bool);
template void get_exx_stress_potential<float, base_device::DEVICE_GPU>(const K_Vectors*,
                                                                       const ModulePW::PW_Basis_K*,
                                                                       ModulePW::PW_Basis*,
                                                                       float*,
                                                                       double,
                                                                       bool,
                                                                       double,
                                                                       int,
                                                                       int);
template void get_exx_stress_potential<double, base_device::DEVICE_GPU>(const K_Vectors*,
                                                                        const ModulePW::PW_Basis_K*,
                                                                        ModulePW::PW_Basis*,
                                                                        double*,
                                                                        double,
                                                                        bool,
                                                                        double,
                                                                        int,
                                                                        int);
#endif
} // namespace hamilt