#include <algorithm>
#include <type_traits>

#include "source_base/global_function.h"
#include "gint_rho.h"
#include "gint_common.h"
#include "phi_operator.h"

namespace ModuleGint
{

void Gint_rho::cal_gint()
{
    ModuleBase::TITLE("Gint", "cal_gint_rho");
    ModuleBase::timer::start("Gint", "cal_gint_rho");
    switch (gint_info_->get_exec_precision())
    {
    case GintPrecision::fp32:
        cal_gint_impl_<float>();
        break;
    case GintPrecision::fp64:
    default:
        cal_gint_impl_<double>();
        break;
    }
    ModuleBase::timer::end("Gint", "cal_gint_rho");
}

template<typename Real>
void Gint_rho::cal_gint_impl_()
{
    std::vector<HContainer<Real>> dm_gint_vec = init_dm_gint_<Real>();
    std::vector<std::vector<Real>> rho_cache(nspin_);
    std::vector<Real*> rho_data(nspin_);
    for (int is = 0; is < nspin_; ++is)
    {
        rho_data[is] = get_rho_data_<Real>(is, rho_cache);
    }
    transfer_dm_2d_to_gint(*gint_info_, dm_vec_, dm_gint_vec);
    cal_rho_(dm_gint_vec, rho_data);
    transfer_rho_cache_<Real>(rho_cache);
}

template<typename Real>
std::vector<HContainer<Real>> Gint_rho::init_dm_gint_() const
{
    std::vector<HContainer<Real>> dm_gint_vec(nspin_);
    for (int is = 0; is < nspin_; is++)
    {
        dm_gint_vec[is] = gint_info_->get_hr<Real>();
    }
    return dm_gint_vec;
}

// Overloaded helpers (C++11-compatible alternative to if constexpr).
// The double overload is preferred by overload resolution when Real=double.

inline double* get_rho_data(double* const* rho, int is, int /*local_mgrid_num*/,
                            std::vector<std::vector<double>>& /*rho_cache*/)
{
    return rho[is];
}

template<typename Real>
Real* get_rho_data(double* const* rho, int is, int local_mgrid_num,
                   std::vector<std::vector<Real>>& rho_cache)
{
    rho_cache[is].resize(local_mgrid_num);
    std::transform(rho[is], rho[is] + local_mgrid_num, rho_cache[is].begin(), [](const double value) {
        return static_cast<Real>(value);
    });
    return rho_cache[is].data();
}

inline void transfer_rho_back(double* const* /*rho*/, int /*nspin*/, int /*local_mgrid_num*/,
                               const std::vector<std::vector<double>>& /*rho_cache*/)
{
    // Nothing to do: double rho was written directly.
}

template<typename Real>
void transfer_rho_back(double* const* rho, int nspin, int local_mgrid_num,
                       const std::vector<std::vector<Real>>& rho_cache)
{
    for (int is = 0; is < nspin; ++is)
    {
        for (int ir = 0; ir < local_mgrid_num; ++ir)
        {
            rho[is][ir] = static_cast<double>(rho_cache[is][ir]);
        }
    }
}



template<typename Real>
Real* Gint_rho::get_rho_data_(int is, std::vector<std::vector<Real>>& rho_cache) const
{
    return get_rho_data(
        rho_, is, gint_info_->get_local_mgrid_num(), rho_cache);
}

template<typename Real>
void Gint_rho::cal_rho_(
    const std::vector<HContainer<Real>>& dm_gint_vec,
    const std::vector<Real*>& rho_data) const
{
#pragma omp parallel
    {
        PhiOperator phi_op;
        std::vector<Real> phi;
        std::vector<Real> phi_dm;
#pragma omp for schedule(dynamic)
        for (int i = 0; i < gint_info_->get_bgrids_num(); i++)
        {
            const auto& biggrid = gint_info_->get_biggrids()[i];
            if (biggrid->get_atoms().empty())
            {
                continue;
            }
            phi_op.set_bgrid(biggrid);
            const int phi_len = phi_op.get_rows() * phi_op.get_cols();
            phi.resize(phi_len);
            phi_dm.resize(phi_len);
            phi_op.set_phi(phi.data());
            for (int is = 0; is < nspin_; is++)
            {
                phi_op.phi_mul_dm(phi.data(), dm_gint_vec[is], is_dm_symm_, phi_dm.data());
                phi_op.phi_dot_phi(phi.data(), phi_dm.data(), rho_data[is]);
            }
        }
    }
}

template<typename Real>
void Gint_rho::transfer_rho_cache_(const std::vector<std::vector<Real>>& rho_cache) const
{
    transfer_rho_back(
        rho_, nspin_, gint_info_->get_local_mgrid_num(), rho_cache);
}

}
