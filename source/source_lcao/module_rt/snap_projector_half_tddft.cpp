#include "snap_projector_half_tddft.h"

#include "source_base/constants.h"
#include "source_base/global_function.h"
#include "source_base/math_integral.h"
#include "source_base/math_lebedev_laikov.h"
#include "source_base/math_polyint.h"
#include "source_base/timer.h"
#include "source_base/ylm.h"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <complex>
#include <map>
#include <memory>
#include <mutex>
#include <set>
#include <string>

namespace module_rt
{
namespace
{
constexpr int default_radial_grid_num = 140;
constexpr int default_lebedev_grid_points = 110;

/**
 * @brief Cached Gauss-Legendre radial grid for a requested grid size.
 */
struct GaussLegendreGrid
{
    explicit GaussLegendreGrid(const int ngrid) : x(ngrid), w(ngrid)
    {
        ModuleBase::Integral::Gauss_Legendre_grid_and_weight(ngrid, x.data(), w.data());
    }

    std::vector<double> x;
    std::vector<double> w;
};

const GaussLegendreGrid& gauss_legendre_grid(const int ngrid)
{
    // Tests may request non-default radial grids, so cache by grid size.
    static std::map<int, std::shared_ptr<const GaussLegendreGrid>> cache;
    static std::mutex cache_mutex;

    std::lock_guard<std::mutex> lock(cache_mutex);
    std::shared_ptr<const GaussLegendreGrid>& grid = cache[ngrid];
    if (!grid)
    {
        grid.reset(new GaussLegendreGrid(ngrid));
    }
    return *grid;
}

/**
 * @brief Owned Lebedev-Laikov angular grid generated at runtime.
 */
struct AngularGridData
{
    explicit AngularGridData(const int ngrid) : x(ngrid), y(ngrid), z(ngrid), w(ngrid)
    {
        ModuleBase::Lebedev_laikov_grid grid(ngrid);
        grid.generate_grid_points();
        const ModuleBase::Vector3<double>* grid_coor = grid.get_grid_coor();
        const double* weight = grid.get_weight();
        for (int i = 0; i < ngrid; ++i)
        {
            x[i] = grid_coor[i].x;
            y[i] = grid_coor[i].y;
            z[i] = grid_coor[i].z;
            w[i] = weight[i];
        }
    }

    std::vector<double> x;
    std::vector<double> y;
    std::vector<double> z;
    std::vector<double> w;
};

/**
 * @brief Non-owning view used by the integration loops.
 */
struct AngularGridView
{
    int size = 0;
    const double* x = nullptr;
    const double* y = nullptr;
    const double* z = nullptr;
    const double* w = nullptr;
};

bool is_supported_lebedev_grid(const int ngrid)
{
    static const std::set<int> supported_grids
        = {6,   14,  26,  38,   50,   74,   86,   110,  146,  170,  194,  230,  266,  302,  350,  434,
           590, 770, 974, 1202, 1454, 1730, 2030, 2354, 2702, 3074, 3470, 3890, 4334, 4802, 5294, 5810};
    return supported_grids.find(ngrid) != supported_grids.end();
}

AngularGridView angular_grid(const int ngrid)
{
    if (!is_supported_lebedev_grid(ngrid))
    {
        ModuleBase::WARNING_QUIT("snap_projector_half_tddft",
                                 "Unsupported Lebedev-Laikov grid size: " + std::to_string(ngrid));
    }

    if (ngrid == default_lebedev_grid_points)
    {
        // Keep the production path on the historical static 110-point table.
        AngularGridView view;
        view.size = default_lebedev_grid_points;
        view.x = ModuleBase::Integral::Lebedev_Laikov_grid110_x;
        view.y = ModuleBase::Integral::Lebedev_Laikov_grid110_y;
        view.z = ModuleBase::Integral::Lebedev_Laikov_grid110_z;
        view.w = ModuleBase::Integral::Lebedev_Laikov_grid110_w;
        return view;
    }

    // Higher-order grids are generated lazily for tests and future callers.
    static std::map<int, std::shared_ptr<const AngularGridData>> cache;
    static std::mutex cache_mutex;

    std::lock_guard<std::mutex> lock(cache_mutex);
    std::shared_ptr<const AngularGridData>& data = cache[ngrid];
    if (!data)
    {
        data.reset(new AngularGridData(ngrid));
    }

    AngularGridView view;
    view.size = ngrid;
    view.x = data->x.data();
    view.y = data->y.data();
    view.z = data->z.data();
    view.w = data->w.data();
    return view;
}

double radial_factor(const ProjectorChannel& channel, const double r, const double w_radial)
{
    const double projector_val
        = ModuleBase::PolyInt::Polynomial_Interpolation(channel.radial_values, channel.mesh, channel.dk, r);

    return projector_val * r * w_radial;
}
} // namespace

void snap_projector_half_tddft(const LCAO_Orbitals& orb,
                               const std::vector<ProjectorChannel>& projector_channels,
                               std::vector<std::vector<std::complex<double>>>& nlm,
                               const ModuleBase::Vector3<double>& R1,
                               const int& T1,
                               const int& L1,
                               const int& m1,
                               const int& N1,
                               const ModuleBase::Vector3<double>& R0,
                               const ModuleBase::Vector3<double>& A,
                               const bool& calc_r,
                               const char* timer_name)
{
    // Preserve the production default while allowing tests to call the overload.
    SnapIntegrationOptions options;
    options.radial_grid_num = default_radial_grid_num;
    options.lebedev_grid_points = default_lebedev_grid_points;
    snap_projector_half_tddft(orb, projector_channels, nlm, R1, T1, L1, m1, N1, R0, A, calc_r, options, timer_name);
}

void snap_projector_half_tddft(const LCAO_Orbitals& orb,
                               const std::vector<ProjectorChannel>& projector_channels,
                               std::vector<std::vector<std::complex<double>>>& nlm,
                               const ModuleBase::Vector3<double>& R1,
                               const int& T1,
                               const int& L1,
                               const int& m1,
                               const int& N1,
                               const ModuleBase::Vector3<double>& R0,
                               const ModuleBase::Vector3<double>& A,
                               const bool& calc_r,
                               const SnapIntegrationOptions& options,
                               const char* timer_name)
{
    ModuleBase::timer::start("module_rt", timer_name);
    if (options.radial_grid_num <= 0)
    {
        ModuleBase::WARNING_QUIT("snap_projector_half_tddft", "The radial grid size must be positive.");
    }
    const int radial_grid_num = options.radial_grid_num;
    const AngularGridView lebedev = angular_grid(options.lebedev_grid_points);

    const int required_size = calc_r ? 4 : 1;
    if (nlm.size() != required_size)
    {
        nlm.resize(required_size);
    }

    int natomwfc = 0;
    std::vector<bool> active(projector_channels.size(), false);

    const double Rcut1 = orb.Phi[T1].getRcut();
    const ModuleBase::Vector3<double> dRa = R0 - R1;
    const double distance10 = dRa.norm();

    bool any_active = false;
    for (int ich = 0; ich < static_cast<int>(projector_channels.size()); ++ich)
    {
        const ProjectorChannel& channel = projector_channels[ich];
        natomwfc += 2 * channel.l + 1;
        if (distance10 <= Rcut1 + channel.rcut)
        {
            active[ich] = true;
            any_active = true;
        }
    }

    for (auto& x: nlm)
    {
        x.assign(natomwfc, 0.0);
    }

    if (natomwfc == 0 || !any_active)
    {
        ModuleBase::timer::end("module_rt", timer_name);
        return;
    }

    // The LCAO orbital is sampled at r + R0 - R1 around the projector center.
    const auto& phi_ln = orb.Phi[T1].PhiLN(L1, N1);
    const int mesh_r1 = phi_ln.getNr();
    const double* psi_1 = phi_ln.getPsi();
    const double dk_1 = phi_ln.getDk();

    const GaussLegendreGrid& gl = gauss_legendre_grid(radial_grid_num);
    std::vector<double> r_radial(radial_grid_num);
    std::vector<double> w_radial(radial_grid_num);

    std::vector<double> A_dot_lebedev(lebedev.size);
    for (int ian = 0; ian < lebedev.size; ++ian)
    {
        A_dot_lebedev[ian] = A.x * lebedev.x[ian] + A.y * lebedev.y[ian] + A.z * lebedev.z[ian];
    }

    std::vector<std::complex<double>> result_angular;
    std::vector<std::complex<double>> res_ang_x;
    std::vector<std::complex<double>> res_ang_y;
    std::vector<std::complex<double>> res_ang_z;
    std::vector<double> rly1((L1 + 1) * (L1 + 1));
    std::vector<std::vector<double>> rly0_cache(lebedev.size);

    int index_offset = 0;
    for (int ich = 0; ich < static_cast<int>(projector_channels.size()); ++ich)
    {
        const ProjectorChannel& channel = projector_channels[ich];
        const int L0 = channel.l;
        const int num_m0 = 2 * L0 + 1;

        if (!active[ich])
        {
            index_offset += num_m0;
            continue;
        }

        assert(channel.mesh > 0);
        assert(channel.radial_values != nullptr);
        assert(channel.radial_grid != nullptr);

        const double r_min = channel.radial_grid[0];
        const double r_max = channel.radial_grid[channel.mesh - 1];
        const double xl = (r_max - r_min) * 0.5;
        const double xmean = (r_max + r_min) * 0.5;

        for (int i = 0; i < radial_grid_num; ++i)
        {
            r_radial[i] = xmean + xl * gl.x[i];
            w_radial[i] = xl * gl.w[i];
        }

        const double A_phase = A * R0;
        const std::complex<double> exp_iAR0 = std::exp(ModuleBase::IMAG_UNIT * A_phase);

        // Y_lm(projector direction) only depends on the angular grid.
        for (int ian = 0; ian < lebedev.size; ++ian)
        {
            ModuleBase::Ylm::rl_sph_harm(L0, lebedev.x[ian], lebedev.y[ian], lebedev.z[ian], rly0_cache[ian]);
        }

        if (result_angular.size() < static_cast<size_t>(num_m0))
        {
            result_angular.resize(num_m0);
            if (calc_r)
            {
                res_ang_x.resize(num_m0);
                res_ang_y.resize(num_m0);
                res_ang_z.resize(num_m0);
            }
        }

        for (int ir = 0; ir < radial_grid_num; ++ir)
        {
            const double r_val = r_radial[ir];

            std::fill(result_angular.begin(), result_angular.begin() + num_m0, 0.0);
            if (calc_r)
            {
                std::fill(res_ang_x.begin(), res_ang_x.begin() + num_m0, 0.0);
                std::fill(res_ang_y.begin(), res_ang_y.begin() + num_m0, 0.0);
                std::fill(res_ang_z.begin(), res_ang_z.begin() + num_m0, 0.0);
            }

            for (int ian = 0; ian < lebedev.size; ++ian)
            {
                const double x = lebedev.x[ian];
                const double y = lebedev.y[ian];
                const double z = lebedev.z[ian];
                const double w_ang = lebedev.w[ian];

                const double rx = r_val * x;
                const double ry = r_val * y;
                const double rz = r_val * z;

                const double tx = rx + dRa.x;
                const double ty = ry + dRa.y;
                const double tz = rz + dRa.z;
                const double tnorm = std::sqrt(tx * tx + ty * ty + tz * tz);

                if (tnorm > Rcut1)
                {
                    continue;
                }

                if (tnorm > 1e-10)
                {
                    const double inv_tnorm = 1.0 / tnorm;
                    ModuleBase::Ylm::rl_sph_harm(L1, tx * inv_tnorm, ty * inv_tnorm, tz * inv_tnorm, rly1);
                }
                else
                {
                    ModuleBase::Ylm::rl_sph_harm(L1, 0.0, 0.0, 1.0, rly1);
                }

                const double phase = r_val * A_dot_lebedev[ian];
                const std::complex<double> exp_iAr = std::exp(ModuleBase::IMAG_UNIT * phase);
                const double interp_psi = ModuleBase::PolyInt::Polynomial_Interpolation(psi_1, mesh_r1, dk_1, tnorm);
                const double ylm_L1_val = rly1[L1 * L1 + m1];
                const std::complex<double> common_factor = exp_iAr * ylm_L1_val * interp_psi * w_ang;

                // Accumulate all magnetic components of the same projector channel.
                const std::vector<double>& rly0_vec = rly0_cache[ian];
                const int offset_L0 = L0 * L0;
                for (int m0 = 0; m0 < num_m0; ++m0)
                {
                    const std::complex<double> term = common_factor * rly0_vec[offset_L0 + m0];
                    result_angular[m0] += term;

                    if (calc_r)
                    {
                        res_ang_x[m0] += term * (rx + R0.x);
                        res_ang_y[m0] += term * (ry + R0.y);
                        res_ang_z[m0] += term * (rz + R0.z);
                    }
                }
            }

            const double factor = radial_factor(channel, r_val, w_radial[ir]);
            int current_idx = index_offset;
            for (int m0 = 0; m0 < num_m0; ++m0)
            {
                nlm[0][current_idx] += factor * result_angular[m0] * exp_iAR0;
                if (calc_r)
                {
                    nlm[1][current_idx] += factor * res_ang_x[m0] * exp_iAR0;
                    nlm[2][current_idx] += factor * res_ang_y[m0] * exp_iAR0;
                    nlm[3][current_idx] += factor * res_ang_z[m0] * exp_iAR0;
                }
                ++current_idx;
            }
        }

        index_offset += num_m0;
    }

    for (auto& dim: nlm)
    {
        for (auto& x: dim)
        {
            x = std::conj(x);
        }
    }

    assert(index_offset == natomwfc);
    ModuleBase::timer::end("module_rt", timer_name);
}

} // namespace module_rt
