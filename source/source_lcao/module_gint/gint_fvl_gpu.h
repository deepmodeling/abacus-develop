#pragma once

#include <memory>
#include <vector>
#include "source_lcao/module_hcontainer/hcontainer.h"
#include "source_base/matrix.h"
#include "gint.h"
#include "gint_info.h"
#include "source_lcao/module_gint/kernel/cuda_mem_wrapper.h"

namespace ModuleGint
{

class Gint_fvl_gpu : public Gint
{
    public:
    Gint_fvl_gpu(
        const int nspin,
        const std::vector<const double*>& vr_eff,
        const std::vector<HContainer<double>*>& dm_vec,
        const bool isforce,
        const bool isstress,
        ModuleBase::matrix* fvl,
        ModuleBase::matrix* svl)
        : nspin_(nspin), vr_eff_(vr_eff), dm_vec_(dm_vec),
          isforce_(isforce), isstress_(isstress), fvl_(fvl), svl_(svl),
          dr3_(gint_info_->get_mgrid_volume()) {}

    void cal_gint();

    private:
    void init_dm_gint_();

    void cal_fvl_svl_();
    
    void transfer_cpu_to_gpu_();
    void transfer_gpu_to_cpu_();
    // input
    const int nspin_;
    std::vector<const double*> vr_eff_;
    std::vector<HContainer<double>*> dm_vec_;
    const bool isforce_;
    const bool isstress_;

    // output
    ModuleBase::matrix* fvl_;
    ModuleBase::matrix* svl_;

    // intermediate variables
    std::vector<HContainer<double>> dm_gint_vec_;

    double dr3_;
    
    // GPU memory
    std::vector<CudaMemWrapper<double>> vr_eff_d_vec_;
    std::vector<CudaMemWrapper<double>> dm_gint_d_vec_;
    CudaMemWrapper<double> fvl_d_;
    CudaMemWrapper<double> svl_d_;
};

}