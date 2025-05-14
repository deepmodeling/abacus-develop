#pragma once

#include <memory>
#include <vector>
#include "module_hamilt_lcao/module_hcontainer/hcontainer.h"
#include "gint.h"
#include "gint_info.h"

namespace ModuleGint
{

class Gint_rho : public Gint
{
    public:
    Gint_rho(
        const std::vector<HContainer<double>*>& dm_vec,
        const int nspin,
        double **rho)
        : dm_vec_(dm_vec), nspin_(nspin), rho_(rho) {};
    
    void cal_gint() override;

    private:
    void init_dm_gint_();

    void cal_rho_();

    // input
    const std::vector<HContainer<double>*> dm_vec_;
    const int nspin_;

    // output
    double **rho_;

    //========================
    // Intermediate variables
    //========================
    std::vector<std::shared_ptr<HContainer<double>>> dm_gint_vec_;
};

}