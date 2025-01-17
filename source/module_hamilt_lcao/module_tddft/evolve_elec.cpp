#include "evolve_elec.h"

#include "evolve_psi.h"
#include "module_base/parallel_reduce.h"
#include "module_base/timer.h"
#include "module_elecstate/module_charge/symmetry_rho.h"
#include "module_hamilt_lcao/hamilt_lcaodft/hamilt_lcao.h"
#include "module_hamilt_lcao/module_dftu/dftu.h"
#include "module_hamilt_pw/hamilt_pwdft/global.h"

namespace module_tddft
{
template <typename Device>
Evolve_elec<Device>::Evolve_elec(){};
template <typename Device>
Evolve_elec<Device>::~Evolve_elec(){};

template <typename Device>
Device* Evolve_elec<Device>::ctx = {};
template <typename Device>
base_device::DEVICE_CPU* Evolve_elec<Device>::cpu_ctx = {};
template <typename Device>
ct::DeviceType Evolve_elec<Device>::ct_device_type = ct::DeviceTypeToEnum<Device>::value;

// this routine only serves for TDDFT using LCAO basis set
template <typename Device>
void Evolve_elec<Device>::solve_psi(const int& istep,
                                    const int nband,
                                    const int nlocal,
                                    hamilt::Hamilt<std::complex<double>>* phm,
                                    Parallel_Orbitals& para_orb,
                                    psi::Psi<std::complex<double>>* psi,
                                    psi::Psi<std::complex<double>>* psi_laststep,
                                    std::complex<double>** Hk_laststep,
                                    std::complex<double>** Sk_laststep,
                                    ModuleBase::matrix& ekb,
                                    int htype,
                                    int propagator,
                                    const int& nks)
{
    ModuleBase::TITLE("Evolve_elec", "solve_psi");
    ModuleBase::timer::tick("Evolve_elec", "solve_psi");

    const int print_matrix = 0;
    // const bool use_tensor = true;
    const bool use_tensor = false;
    const bool use_lapack = true;
    // const bool use_lapack = false;

    for (int ik = 0; ik < nks; ik++)
    {
        phm->updateHk(ik);

        ModuleBase::timer::tick("Efficiency", "evolve_k");
        psi->fix_k(ik);
        psi_laststep->fix_k(ik);
        if (htype == 0)
        {
            evolve_psi(nband,
                       nlocal,
                       &(para_orb),
                       phm,
                       psi[0].get_pointer(),
                       psi_laststep[0].get_pointer(),
                       nullptr,
                       nullptr,
                       &(ekb(ik, 0)),
                       htype,
                       propagator,
                       print_matrix);
        }
        else if (htype == 1)
        {
            if (!use_tensor)
            {
                evolve_psi(nband,
                           nlocal,
                           &(para_orb),
                           phm,
                           psi[0].get_pointer(),
                           psi_laststep[0].get_pointer(),
                           Hk_laststep[ik],
                           Sk_laststep[ik],
                           &(ekb(ik, 0)),
                           htype,
                           propagator,
                           print_matrix);
                // std::cout << "Print ekb: " << std::endl;
                // ekb.print(std::cout);
            }
            else
            {
                // std::cout << "nband = " << nband << std::endl;
                // std::cout << "psi->get_nbands() = " << psi->get_nbands() << std::endl;
                // std::cout << "nlocal = " << nlocal << std::endl;
                // std::cout << "psi->get_nbasis() = " << psi->get_nbasis() << std::endl;
                // std::cout << "ekb.nr = " << ekb.nr << std::endl;
                // std::cout << "ekb.nc = " << ekb.nc << std::endl;

                // Create Tensor for psi_k, psi_k_laststep, H_laststep, S_laststep, ekb
                ct::Tensor psi_k_tensor(ct::DataType::DT_COMPLEX_DOUBLE,
                                        ct_device_type,
                                        ct::TensorShape({psi->get_nbands(), psi->get_nbasis()}));
                ct::Tensor psi_k_laststep_tensor(ct::DataType::DT_COMPLEX_DOUBLE,
                                                 ct_device_type,
                                                 ct::TensorShape({psi->get_nbands(), psi->get_nbasis()}));
                ct::Tensor H_laststep_tensor(ct::DataType::DT_COMPLEX_DOUBLE,
                                             ct_device_type,
                                             ct::TensorShape({para_orb.nloc}));
                ct::Tensor S_laststep_tensor(ct::DataType::DT_COMPLEX_DOUBLE,
                                             ct_device_type,
                                             ct::TensorShape({para_orb.nloc}));
                ct::Tensor ekb_tensor(ct::DataType::DT_DOUBLE, ct_device_type, ct::TensorShape({nband}));

                // Syncronize data from CPU to Device
                syncmem_complex_h2d_op()(ctx,
                                         cpu_ctx,
                                         psi_k_tensor.data<std::complex<double>>(),
                                         psi[0].get_pointer(),
                                         psi->get_nbands() * psi->get_nbasis());
                syncmem_complex_h2d_op()(ctx,
                                         cpu_ctx,
                                         psi_k_laststep_tensor.data<std::complex<double>>(),
                                         psi_laststep[0].get_pointer(),
                                         psi->get_nbands() * psi->get_nbasis());
                syncmem_complex_h2d_op()(ctx,
                                         cpu_ctx,
                                         H_laststep_tensor.data<std::complex<double>>(),
                                         Hk_laststep[ik],
                                         para_orb.nloc);
                syncmem_complex_h2d_op()(ctx,
                                         cpu_ctx,
                                         S_laststep_tensor.data<std::complex<double>>(),
                                         Sk_laststep[ik],
                                         para_orb.nloc);
                syncmem_double_h2d_op()(ctx, cpu_ctx, ekb_tensor.data<double>(), &(ekb(ik, 0)), nband);

                evolve_psi_tensor<Device>(nband,
                                          nlocal,
                                          &(para_orb),
                                          phm,
                                          psi_k_tensor,
                                          psi_k_laststep_tensor,
                                          H_laststep_tensor,
                                          S_laststep_tensor,
                                          ekb_tensor,
                                          htype,
                                          propagator,
                                          print_matrix,
                                          use_lapack);

                // Syncronize data from Device to CPU
                syncmem_complex_d2h_op()(cpu_ctx,
                                         ctx,
                                         psi[0].get_pointer(),
                                         psi_k_tensor.data<std::complex<double>>(),
                                         psi->get_nbands() * psi->get_nbasis());
                syncmem_complex_d2h_op()(cpu_ctx,
                                         ctx,
                                         psi_laststep[0].get_pointer(),
                                         psi_k_laststep_tensor.data<std::complex<double>>(),
                                         psi->get_nbands() * psi->get_nbasis());
                syncmem_complex_d2h_op()(cpu_ctx,
                                         ctx,
                                         Hk_laststep[ik],
                                         H_laststep_tensor.data<std::complex<double>>(),
                                         para_orb.nloc);
                syncmem_complex_d2h_op()(cpu_ctx,
                                         ctx,
                                         Sk_laststep[ik],
                                         S_laststep_tensor.data<std::complex<double>>(),
                                         para_orb.nloc);
                syncmem_double_d2h_op()(cpu_ctx, ctx, &(ekb(ik, 0)), ekb_tensor.data<double>(), nband);

                // std::cout << "Print ekb tensor: " << std::endl;
                // ekb.print(std::cout);
            }
        }
        else
        {
            std::cout << "method of htype is wrong" << std::endl;
        }

        ModuleBase::timer::tick("Efficiency", "evolve_k");
    } // end k

    ModuleBase::timer::tick("Evolve_elec", "solve_psi");
    return;
}

template class Evolve_elec<base_device::DEVICE_CPU>;
#if ((defined __CUDA) /* || (defined __ROCM) */)
template class Evolve_elec<base_device::DEVICE_GPU>;
#endif
} // namespace module_tddft