#include "evolve_elec.h"

#include "evolve_psi.h"
#include "source_base/parallel_reduce.h"
#include "source_base/timer.h"
#include "source_estate/module_charge/symmetry_rho.h"
#include "source_lcao/hamilt_lcao.h"
#include "source_lcao/module_dftu/dftu.h"

namespace module_rt
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
                                    ct::Tensor& Hk_laststep,
                                    ct::Tensor& Sk_laststep,
                                    ModuleBase::matrix& ekb,
                                    std::ofstream& ofs_running,
                                    const int propagator,
                                    const bool use_tensor,
                                    const bool use_lapack)
{
    ModuleBase::TITLE("Evolve_elec", "solve_psi");
    ModuleBase::timer::start("Evolve_elec", "solve_psi");

    // Control the print of matrix to running_md.log
    const int print_matrix = 0;

    // Multi-GPU support
    CublasMpResources cublas_res;
#ifdef __CUBLASMP
    init_cublasmp_resources(cublas_res, MPI_COMM_WORLD, para_orb.desc);
#endif

    for (int ik = 0; ik < nks; ik++)
    {
        phm->updateHk(ik);

        ModuleBase::timer::start("TD_Efficiency", "evolve_k");
        psi->fix_k(ik);
        psi_laststep->fix_k(ik);

        if (!use_tensor)
        {
            const int len_HS_laststep = use_lapack ? nlocal * nlocal : para_orb.nloc;
            evolve_psi(nband,
                       nlocal,
                       &(para_orb),
                       phm,
                       psi[0].get_pointer(),
                       psi_laststep[0].get_pointer(),
                       Hk_laststep.data<std::complex<double>>() + ik * len_HS_laststep,
                       Sk_laststep.data<std::complex<double>>() + ik * len_HS_laststep,
                       &(ekb(ik, 0)),
                       propagator,
                       ofs_running,
                       print_matrix);
            // GlobalV::ofs_running << "Print ekb: " << std::endl;
            // ekb.print(GlobalV::ofs_running);
        }
        else
        {
            ModuleBase::timer::start("TD_Efficiency", "host_device_comm");

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
            module_rt::Matrix_g<std::complex<double>> psi_g;
            module_rt::Matrix_g<std::complex<double>> psi_laststep_g;

            // Prepare host pointers for psi and psi_laststep
            std::complex<double>* p_psi_host = nullptr;
            std::complex<double>* p_psi_last_host = nullptr;

            if (use_lapack)
            {
#ifdef __MPI
                int myid = 0;
                const int root_proc = 0;
                int num_procs = 1;
                MPI_Comm_rank(MPI_COMM_WORLD, &myid);
                MPI_Comm_size(MPI_COMM_WORLD, &num_procs);

                if (num_procs == 1)
                {
                    // Single process: directly point to local data without gather
                    p_psi_host = psi[0].get_pointer();
                    p_psi_last_host = psi_laststep[0].get_pointer();
                }
                else
                {
                    // Multiple processes: gather data to the root process (myid == 0) and point to the gathered data
                    gatherPsi(myid, root_proc, psi[0].get_pointer(), para_orb, psi_g);
                    gatherPsi(myid, root_proc, psi_laststep[0].get_pointer(), para_orb, psi_laststep_g);

                    if (myid == root_proc)
                    {
                        p_psi_host = psi_g.p.get();
                        p_psi_last_host = psi_laststep_g.p.get();
                    }
                }

                // Only the root process (myid == 0) performs the copy
                if (myid == root_proc)
                {
                    syncmem_complex_h2d_op()(psi_k_tensor.data<std::complex<double>>(),
                                             p_psi_host,
                                             len_psi_k_1 * len_psi_k_2);
                    syncmem_complex_h2d_op()(psi_k_laststep_tensor.data<std::complex<double>>(),
                                             p_psi_last_host,
                                             len_psi_k_1 * len_psi_k_2);
                }
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
                                     Hk_laststep.data<std::complex<double>>() + ik * len_HS_laststep,
                                     len_HS_laststep);
            syncmem_complex_h2d_op()(S_laststep_tensor.data<std::complex<double>>(),
                                     Sk_laststep.data<std::complex<double>>() + ik * len_HS_laststep,
                                     len_HS_laststep);
            syncmem_double_h2d_op()(ekb_tensor.data<double>(), &(ekb(ik, 0)), nband);

            ModuleBase::timer::end("TD_Efficiency", "host_device_comm");

            evolve_psi_tensor<Device>(nband,
                                      nlocal,
                                      &(para_orb),
                                      phm,
                                      psi_k_tensor,
                                      psi_k_laststep_tensor,
                                      H_laststep_tensor,
                                      S_laststep_tensor,
                                      ekb_tensor,
                                      propagator,
                                      ofs_running,
                                      print_matrix,
                                      use_lapack,
                                      cublas_res);

            ModuleBase::timer::start("TD_Efficiency", "host_device_comm");
            // Need to distribute global psi back to all processes
            if (use_lapack)
            {
#ifdef __MPI
                int myid = 0;
                int num_procs = 1;
                MPI_Comm_rank(MPI_COMM_WORLD, &myid);
                MPI_Comm_size(MPI_COMM_WORLD, &num_procs);

                if (myid == 0)
                {
                    syncmem_complex_d2h_op()(p_psi_host,
                                             psi_k_tensor.data<std::complex<double>>(),
                                             len_psi_k_1 * len_psi_k_2);
                    syncmem_complex_d2h_op()(p_psi_last_host,
                                             psi_k_laststep_tensor.data<std::complex<double>>(),
                                             len_psi_k_1 * len_psi_k_2);
                }

                // If it's multi-process, distribute back; if it's single-process, the data is already in psi[0]
                if (num_procs > 1)
                {
                    distributePsi(para_orb, psi[0].get_pointer(), psi_g);
                    distributePsi(para_orb, psi_laststep[0].get_pointer(), psi_laststep_g);
                }
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
            syncmem_complex_d2h_op()(Hk_laststep.data<std::complex<double>>() + ik * len_HS_laststep,
                                     H_laststep_tensor.data<std::complex<double>>(),
                                     len_HS_laststep);
            syncmem_complex_d2h_op()(Sk_laststep.data<std::complex<double>>() + ik * len_HS_laststep,
                                     S_laststep_tensor.data<std::complex<double>>(),
                                     len_HS_laststep);
            syncmem_double_d2h_op()(&(ekb(ik, 0)), ekb_tensor.data<double>(), nband);

#ifdef __MPI
            const int root_proc = 0;
            if (use_lapack)
            {
                // Synchronize ekb to all MPI processes
                MPI_Bcast(&(ekb(ik, 0)), nband, MPI_DOUBLE, root_proc, MPI_COMM_WORLD);
            }
#endif

            ModuleBase::timer::end("TD_Efficiency", "host_device_comm");

            // GlobalV::ofs_running << "Print ekb: " << std::endl;
            // ekb.print(GlobalV::ofs_running);
        }

        ModuleBase::timer::end("TD_Efficiency", "evolve_k");
    } // end k

#ifdef __CUBLASMP
    finalize_cublasmp_resources(cublas_res);
#endif

    ModuleBase::timer::end("Evolve_elec", "solve_psi");
    return;
}

template class Evolve_elec<base_device::DEVICE_CPU>;
#if ((defined __CUDA) /* || (defined __ROCM) */)
template class Evolve_elec<base_device::DEVICE_GPU>;
#endif
} // namespace module_rt
