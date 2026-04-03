#pragma once

#include <memory>
#include <vector>
#include "source_lcao/module_hcontainer/hcontainer.h"
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
        double **rho,
        bool is_dm_symm = true)
        : dm_vec_(dm_vec), nspin_(nspin), is_dm_symm_(is_dm_symm), rho_(rho) {}
    
    void cal_gint();

    private:
    template<typename Real>
    void cal_gint_impl_();

    template<typename Real>
    std::vector<HContainer<Real>> init_dm_gint_() const;

    template<typename Real>
    Real* get_rho_data_(int is, std::vector<std::vector<Real>>& rho_cache) const;

    template<typename Real>
    void cal_rho_(
        const std::vector<HContainer<Real>>& dm_gint_vec,
        const std::vector<Real*>& rho_data) const;

    template<typename Real>
    void transfer_rho_cache_(const std::vector<std::vector<Real>>& rho_cache) const;

    // input
    const std::vector<HContainer<double>*> dm_vec_;
    const int nspin_;
    
    // if true, it means the DMR matrix is symmetric,
    // which leads to faster computations compared to the asymmetric case.
    const bool is_dm_symm_;

    // output
    double ** rho_ = nullptr;
};

}
