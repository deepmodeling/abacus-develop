#pragma once

#include <memory>
#include <vector>
#include "module_hamilt_lcao/module_hcontainer/hcontainer.h"
#include "gint.h"
#include "gint_info.h"
#include "module_hamilt_lcao/module_gint/temp_gint/kernel/cuda_mem_wrapper.h"

namespace ModuleGint
{

class Gint_rho_gpu: public Gint
{
    public:
    Gint_rho_gpu(
        const std::vector<HContainer<double>*>& dm_vec,
        const int nspin,
        double **rho,
        bool is_dm_symm = true)
        : dm_vec_(dm_vec), nspin_(nspin), rho_(rho), is_dm_symm_(is_dm_symm) {}
    
    void cal_gint();

    private:
    void init_dm_gint_();

    void cal_rho_();

    void transfer_cpu_to_gpu_();

    void transfer_gpu_to_cpu_();

    // input
    const std::vector<HContainer<double>*> dm_vec_;
    const int nspin_;

    // if true, it means the DMR matrix is symmetric,
    // which leads to faster computations compared to the asymmetric case.
    const bool is_dm_symm_;

    // output
    double **rho_;

    // Intermediate variables
    std::vector<HContainer<double>> dm_gint_vec_;

    std::vector<CudaMemWrapper<double>> dm_gint_d_vec_;
    std::vector<CudaMemWrapper<double>> rho_d_vec_;
};

}