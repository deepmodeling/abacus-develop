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
ct::DeviceType Evolve_elec<Device>::ct_device_type = ct::DeviceTypeToEnum<Device>::value;

// this routine only serves for TDDFT using LCAO basis set
template <typename Device>
void Evolve_elec<Device>::solve_psi(const int& istep,
                                    const int nband,
                                    const int nlocal,
                                    const int& nks,
                                    hamilt::Hamilt<std::complex<double>>* phm,
                                    Parallel_Orbitals& para_orb,
                                    psi::Psi<std::complex<double>>* psi,
                                    psi::Psi<std::complex<double>>* psi_laststep,
                                    std::complex<double>** Hk_laststep,
                                    std::complex<double>** Sk_laststep,
                                    ModuleBase::matrix& ekb,
                                    std::ofstream& ofs_running,
                                    const int htype,
                                    const int propagator,
                                    const bool use_tensor,
                                    const bool use_lapack)
{
    ModuleBase::TITLE("Evolve_elec", "solve_psi");
    ModuleBase::timer::tick("Evolve_elec", "solve_psi");

    // Control the print of matrix to running_md.log
    const int print_matrix = 0;

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
                       ofs_running,
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
                           ofs_running,
                           print_matrix);
                // std::cout << "Print ekb: " << std::endl;
                // ekb.print(std::cout);
            }
            else
            {
                const int len_psi_k_1 = use_lapack ? nband : psi->get_nbands();
                const int len_psi_k_2 = use_lapack ? nlocal : psi->get_nbasis();
                const int len_HS_laststep = use_lapack ? nlocal * nlocal : para_orb.nloc;

                // Create Tensor for psi_k, psi_k_laststep, H_laststep, S_laststep, ekb
                ct::Tensor psi_k_tensor(ct::DataType::DT_COMPLEX_DOUBLE,
                                        ct_device_type,
                                        ct::TensorShape({len_psi_k_1, len_psi_k_2}));
                ct::Tensor psi_k_laststep_tensor(ct::DataType::DT_COMPLEX_DOUBLE,
                                                 ct_device_type,
                                                 ct::TensorShape({len_psi_k_1, len_psi_k_2}));
                ct::Tensor H_laststep_tensor(ct::DataType::DT_COMPLEX_DOUBLE,
                                             ct_device_type,
                                             ct::TensorShape({len_HS_laststep}));
                ct::Tensor S_laststep_tensor(ct::DataType::DT_COMPLEX_DOUBLE,
                                             ct_device_type,
                                             ct::TensorShape({len_HS_laststep}));
                ct::Tensor ekb_tensor(ct::DataType::DT_DOUBLE, ct_device_type, ct::TensorShape({nband}));

                // Global psi
                ModuleESolver::Matrix_g<std::complex<double>> psi_g;
                ModuleESolver::Matrix_g<std::complex<double>> psi_laststep_g;

                if (use_lapack)
                {
                    // Need to gather the psi to the root process on CPU
                    // H_laststep and S_laststep are already gathered in esolver_ks_lcao_tddft.cpp
#ifdef __MPI
                    // Access the rank of the calling process in the communicator
                    int myid = 0;
                    int root_proc = 0;
                    MPI_Comm_rank(MPI_COMM_WORLD, &myid);

                    // Gather psi to the root process
                    gatherPsi(myid, root_proc, psi[0].get_pointer(), para_orb, psi_g);
                    gatherPsi(myid, root_proc, psi_laststep[0].get_pointer(), para_orb, psi_laststep_g);

                    // Syncronize data from CPU to Device
                    syncmem_complex_h2d_op()(psi_k_tensor.data<std::complex<double>>(),
                                             psi_g.p.get(),
                                             len_psi_k_1 * len_psi_k_2);
                    syncmem_complex_h2d_op()(psi_k_laststep_tensor.data<std::complex<double>>(),
                                             psi_laststep_g.p.get(),
                                             len_psi_k_1 * len_psi_k_2);
#endif
                }
                else
                {
                    // Syncronize data from CPU to Device
                    syncmem_complex_h2d_op()(psi_k_tensor.data<std::complex<double>>(),
                                             psi[0].get_pointer(),
                                             len_psi_k_1 * len_psi_k_2);
                    syncmem_complex_h2d_op()(psi_k_laststep_tensor.data<std::complex<double>>(),
                                             psi_laststep[0].get_pointer(),
                                             len_psi_k_1 * len_psi_k_2);
                }

                syncmem_complex_h2d_op()(H_laststep_tensor.data<std::complex<double>>(),
                                         Hk_laststep[ik],
                                         len_HS_laststep);
                syncmem_complex_h2d_op()(S_laststep_tensor.data<std::complex<double>>(),
                                         Sk_laststep[ik],
                                         len_HS_laststep);
                syncmem_double_h2d_op()(ekb_tensor.data<double>(), &(ekb(ik, 0)), nband);

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
                                          ofs_running,
                                          print_matrix,
                                          use_lapack);

                // Need to distribute global psi back to all processes
                if (use_lapack)
                {
#ifdef __MPI
                    // Syncronize data from Device to CPU
                    syncmem_complex_d2h_op()(psi_g.p.get(),
                                             psi_k_tensor.data<std::complex<double>>(),
                                             len_psi_k_1 * len_psi_k_2);
                    syncmem_complex_d2h_op()(psi_laststep_g.p.get(),
                                             psi_k_laststep_tensor.data<std::complex<double>>(),
                                             len_psi_k_1 * len_psi_k_2);

                    // Distribute psi to all processes
                    distributePsi(para_orb, psi[0].get_pointer(), psi_g);
                    distributePsi(para_orb, psi_laststep[0].get_pointer(), psi_laststep_g);
#endif
                }
                else
                {
                    // Syncronize data from Device to CPU
                    syncmem_complex_d2h_op()(psi[0].get_pointer(),
                                             psi_k_tensor.data<std::complex<double>>(),
                                             len_psi_k_1 * len_psi_k_2);
                    syncmem_complex_d2h_op()(psi_laststep[0].get_pointer(),
                                             psi_k_laststep_tensor.data<std::complex<double>>(),
                                             len_psi_k_1 * len_psi_k_2);
                }
                syncmem_complex_d2h_op()(Hk_laststep[ik],
                                         H_laststep_tensor.data<std::complex<double>>(),
                                         len_HS_laststep);
                syncmem_complex_d2h_op()(Sk_laststep[ik],
                                         S_laststep_tensor.data<std::complex<double>>(),
                                         len_HS_laststep);
                syncmem_double_d2h_op()(&(ekb(ik, 0)), ekb_tensor.data<double>(), nband);

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