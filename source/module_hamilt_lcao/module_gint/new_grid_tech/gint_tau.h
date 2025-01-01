#pragma once
#include <memory>
#include <vector>
#include "module_hamilt_lcao/module_hcontainer/hcontainer.h"
#include "gint.h"
#include "gint_info.h"

namespace ModuleGint
{
    
class Gint_tau : public Gint
{
    public:
    gint_tau(
        const std::vector<double*>& DMR_vec,
        const int nspin,
        std::vector<double*>& tau);
    
    void cal_gint() override;
    
    private:
    void init_DMRGint_();

    void cal_tau_();

    // input
    const std::vector<HContainer<double>*> DMR_vec_;
    const int nspin_;

    // output
    double **kin_;

    //========================
    // Intermediate variables
    //========================
    std::vector<std::shared_ptr<HContainer<double>>> DMRGint_vec_;
};

} // namespace ModuleGint
