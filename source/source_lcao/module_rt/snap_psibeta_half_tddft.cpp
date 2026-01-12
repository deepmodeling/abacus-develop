#include "snap_psibeta_half_tddft.h"

#include "source_base/constants.h"
#include "source_base/math_integral.h"
#include "source_base/math_polyint.h"
#include "source_base/timer.h"
#include "source_base/ylm.h"

#include <cmath>
#include <complex>
#include <vector>

namespace module_rt
{

/**
 * @brief Initialize Gauss-Legendre grid points and weights.
 *        Thread-safe initialization using static local variable.
 *
 * @param grid_size Number of grid points (140)
 * @param gl_x Output: Grid points in [-1, 1]
 * @param gl_w Output: Weights
 */
static void init_gauss_legendre_grid(int grid_size, std::vector<double>& gl_x, std::vector<double>& gl_w)
{
    static bool init = false;
// Thread-safe initialization
#pragma omp critical(init_gauss_legendre)
    {
        if (!init)
        {
            ModuleBase::Integral::Gauss_Legendre_grid_and_weight(grid_size, gl_x.data(), gl_w.data());
            init = true;
        }
    }
}

/**
 * @brief Main function to calculate overlap integrals <phi|exp^{-iAr}|beta>
 *        and its derivatives (if calc_r is true).
 *
 *        This function integrates the overlap between a local orbital phi (at R1)
 *        and a non-local projector beta (at R0), modulated by a plane-wave-like phase factor
 *        exp^{-iAr}, where A is a vector potential.
 *
 * @param orb LCAO Orbitals information
 * @param infoNL_ Non-local pseudopotential information
 * @param nlm Output:
 *            nlm[0] : <phi|exp^{-iAr}|beta>
 *            nlm[1, 2, 3] : <phi|r_a * exp^{-iAr}|beta>, a = x, y, z (if calc_r=true)
 * @param R1 Position of atom 1 (orbital phi)
 * @param T1 Type of atom 1
 * @param L1 Angular momentum of orbital phi
 * @param m1 Magnetic quantum number of orbital phi
 * @param N1 Radial quantum number of orbital phi
 * @param R0 Position of atom 0 (projector beta)
 * @param T0 Type of atom 0
 * @param A Vector potential A (or related field vector)
 * @param calc_r Whether to calculate position operator matrix elements
 */
void snap_psibeta_half_tddft(const LCAO_Orbitals& orb,
                             const InfoNonlocal& infoNL_,
                             std::vector<std::vector<std::complex<double>>>& nlm,
                             const ModuleBase::Vector3<double>& R1,
                             const int& T1,
                             const int& L1,
                             const int& m1,
                             const int& N1,
                             const ModuleBase::Vector3<double>& R0,
                             const int& T0,
                             const ModuleBase::Vector3<double>& A,
                             const bool& calc_r)
{
    ModuleBase::timer::tick("module_rt", "snap_psibeta_half_tddft");

    // 1. Initialization and Early Exits
    const int nproj = infoNL_.nproj[T0];

    // Resize output vector based on whether position operator matrix elements are needed
    int required_size = calc_r ? 4 : 1;
    if (nlm.size() != required_size)
        nlm.resize(required_size);

    if (nproj == 0)
        return;

    // 2. Determine total number of projectors and identify active ones based on cutoff
    int natomwfc = 0;
    std::vector<bool> calproj(nproj, false);

    const double Rcut1 = orb.Phi[T1].getRcut();
    const ModuleBase::Vector3<double> dRa = R0 - R1;
    const double distance10 = dRa.norm();

    bool any_active = false;
    for (int ip = 0; ip < nproj; ip++)
    {
        const int L0 = infoNL_.Beta[T0].Proj[ip].getL();
        natomwfc += 2 * L0 + 1;

        const double Rcut0 = infoNL_.Beta[T0].Proj[ip].getRcut();
        if (distance10 <= (Rcut1 + Rcut0))
        {
            calproj[ip] = true;
            any_active = true;
        }
    }

    // Initialize output values to zero and resize inner vectors
    for (auto& x: nlm)
    {
        x.assign(natomwfc, 0.0);
    }

    if (!any_active)
    {
        ModuleBase::timer::tick("module_rt", "snap_psibeta_half_tddft");
        return;
    }

    // 3. Prepare Orbital Data (Phi)
    const auto& phi_ln = orb.Phi[T1].PhiLN(L1, N1);
    const int mesh_r1 = phi_ln.getNr();
    const double* psi_1 = phi_ln.getPsi();
    const double dk_1 = phi_ln.getDk();

    // 4. Prepare Integration Grids
    const int radial_grid_num = 140;
    const int angular_grid_num = 110;

    // Cached standard Gauss-Legendre grid
    static std::vector<double> gl_x(radial_grid_num);
    static std::vector<double> gl_w(radial_grid_num);
    init_gauss_legendre_grid(radial_grid_num, gl_x, gl_w);

    // Buffers for mapped radial grid
    std::vector<double> r_radial(radial_grid_num);
    std::vector<double> w_radial(radial_grid_num);

    // Precompute A dot r_angular (A * u_angle) for the Lebedev grid
    std::vector<double> A_dot_lebedev(angular_grid_num);
    for (int ian = 0; ian < angular_grid_num; ++ian)
    {
        A_dot_lebedev[ian] = A.x * ModuleBase::Integral::Lebedev_Laikov_grid110_x[ian]
                             + A.y * ModuleBase::Integral::Lebedev_Laikov_grid110_y[ian]
                             + A.z * ModuleBase::Integral::Lebedev_Laikov_grid110_z[ian];
    }

    // Reuseable buffers for inner loops to avoid allocation
    std::vector<std::complex<double>> result_angular; // Accumulator for angular integration
    // Accumulators for position operator components
    std::vector<std::complex<double>> res_ang_x, res_ang_y, res_ang_z;

    std::vector<double> rly1((L1 + 1) * (L1 + 1));                 // Spherical harmonics buffer for L1
    std::vector<std::vector<double>> rly0_cache(angular_grid_num); // Cache for L0 Ylm

    // 5. Loop over Projectors (Beta)
    int index_offset = 0;
    for (int nb = 0; nb < nproj; nb++)
    {
        const int L0 = infoNL_.Beta[T0].Proj[nb].getL();
        const int num_m0 = 2 * L0 + 1;

        if (!calproj[nb])
        {
            index_offset += num_m0;
            continue;
        }

        const auto& proj = infoNL_.Beta[T0].Proj[nb];
        const int mesh_r0 = proj.getNr();
        const double* beta_r = proj.getBeta_r();
        const double* radial0 = proj.getRadial();
        const double dk_0 = proj.getDk();
        const double Rcut0 = proj.getRcut();

        // 5.1 Map Gauss-Legendre grid to radial interval [r_min, r_max]
        double r_min = radial0[0];
        double r_max = radial0[mesh_r0 - 1];
        double xl = (r_max - r_min) * 0.5;
        double xmean = (r_max + r_min) * 0.5;

        for (int i = 0; i < radial_grid_num; ++i)
        {
            r_radial[i] = xmean + xl * gl_x[i];
            w_radial[i] = xl * gl_w[i];
        }

        const double A_phase = A * R0;
        const std::complex<double> exp_iAR0 = std::exp(ModuleBase::IMAG_UNIT * A_phase);

        // 5.2 Precompute Spherical Harmonics (Ylm) for L0 on angular grid
        // Since L0 changes with projector, we compute this per projector loop.
        for (int ian = 0; ian < angular_grid_num; ++ian)
        {
            ModuleBase::Ylm::rl_sph_harm(L0,
                                         ModuleBase::Integral::Lebedev_Laikov_grid110_x[ian],
                                         ModuleBase::Integral::Lebedev_Laikov_grid110_y[ian],
                                         ModuleBase::Integral::Lebedev_Laikov_grid110_z[ian],
                                         rly0_cache[ian]);
        }

        // Resize accumulators if needed
        if (result_angular.size() < num_m0)
        {
            result_angular.resize(num_m0);
            if (calc_r)
            {
                res_ang_x.resize(num_m0);
                res_ang_y.resize(num_m0);
                res_ang_z.resize(num_m0);
            }
        }

        // 5.3 Radial Integration Loop
        for (int ir = 0; ir < radial_grid_num; ir++)
        {
            const double r_val = r_radial[ir];

            // Reset angular accumulators for this radial shell
            std::fill(result_angular.begin(), result_angular.begin() + num_m0, 0.0);
            if (calc_r)
            {
                std::fill(res_ang_x.begin(), res_ang_x.begin() + num_m0, 0.0);
                std::fill(res_ang_y.begin(), res_ang_y.begin() + num_m0, 0.0);
                std::fill(res_ang_z.begin(), res_ang_z.begin() + num_m0, 0.0);
            }

            // 5.4 Angular Integration Loop (Lebedev Grid)
            for (int ian = 0; ian < angular_grid_num; ian++)
            {
                const double x = ModuleBase::Integral::Lebedev_Laikov_grid110_x[ian];
                const double y = ModuleBase::Integral::Lebedev_Laikov_grid110_y[ian];
                const double z = ModuleBase::Integral::Lebedev_Laikov_grid110_z[ian];
                const double w_ang = ModuleBase::Integral::Lebedev_Laikov_grid110_w[ian];

                // Vector r = r_val * u_angle
                double rx = r_val * x;
                double ry = r_val * y;
                double rz = r_val * z;

                // Vector r' = r + R0 - R1 = r + dRa
                double tx = rx + dRa.x;
                double ty = ry + dRa.y;
                double tz = rz + dRa.z;

                double tnorm = std::sqrt(tx * tx + ty * ty + tz * tz);

                // If r' is outside the cutoff of Phi(r'), skip
                if (tnorm > Rcut1)
                    continue;

                // Compute Ylm for L1 at direction r'
                if (tnorm > 1e-10)
                {
                    double inv_tnorm = 1.0 / tnorm;
                    ModuleBase::Ylm::rl_sph_harm(L1, tx * inv_tnorm, ty * inv_tnorm, tz * inv_tnorm, rly1);
                }
                else
                {
                    // At origin, only Y_00 is non-zero (if using real spherical harmonics convention)
                    ModuleBase::Ylm::rl_sph_harm(L1, 0.0, 0.0, 1.0, rly1);
                }

                // Calculate common phase and weight factor
                // phase = A * r = r_val * (A * u_angle)
                const double phase = r_val * A_dot_lebedev[ian];
                const std::complex<double> exp_iAr = std::exp(ModuleBase::IMAG_UNIT * phase);

                // Interpolate Psi at |r'|
                double interp_psi = ModuleBase::PolyInt::Polynomial_Interpolation(psi_1, mesh_r1, dk_1, tnorm);

                const int offset_L1 = L1 * L1 + m1;
                const double ylm_L1_val = rly1[offset_L1];

                // Combined factor: exp(iAr) * Y_L1m1(r') * Psi(|r'|) * weight_angle
                const std::complex<double> common_factor = exp_iAr * ylm_L1_val * interp_psi * w_ang;

                // Retrieve precomputed Y_L0m0(r)
                const std::vector<double>& rly0_vec = rly0_cache[ian];
                const int offset_L0 = L0 * L0;

                // Accumulate results for all m0 components
                for (int m0 = 0; m0 < num_m0; m0++)
                {
                    std::complex<double> term = common_factor * rly0_vec[offset_L0 + m0];
                    result_angular[m0] += term;

                    if (calc_r)
                    {
                        // Position operator r_op = r + R0
                        // Note: Term involves (r_op)_a * exp(...).
                        double r_op_x = rx + R0.x;
                        double r_op_y = ry + R0.y;
                        double r_op_z = rz + R0.z;

                        res_ang_x[m0] += term * r_op_x;
                        res_ang_y[m0] += term * r_op_y;
                        res_ang_z[m0] += term * r_op_z;
                    }
                }
            } // End Angular Loop

            // 5.5 Combine Radial and Angular parts
            // Interpolate Beta(|r|)
            // Note: The original code implies beta_r stores values that might need scaling or are just the function
            // values. Typically radial integration is \int f(r) r^2 dr. Here we have factor: beta_val * r_radial[ir] *
            // w_radial[ir] w_radial includes the Jacobian for the change of variable from [-1,1] to [r_min, r_max]. The
            // extra r_radial[ir] suggests either beta is stored as r*beta, or we are doing \int ... r dr (2D?), or
            // Jacobian r^2 is split. Assuming original logic is correct.

            double beta_val = ModuleBase::PolyInt::Polynomial_Interpolation(beta_r, mesh_r0, dk_0, r_radial[ir]);

            double radial_factor = beta_val * r_radial[ir] * w_radial[ir];

            int current_idx = index_offset;
            for (int m0 = 0; m0 < num_m0; m0++)
            {
                // Final accumulation into global nlm array
                // Add phase exp(i A * R0)
                nlm[0][current_idx] += radial_factor * result_angular[m0] * exp_iAR0;

                if (calc_r)
                {
                    nlm[1][current_idx] += radial_factor * res_ang_x[m0] * exp_iAR0;
                    nlm[2][current_idx] += radial_factor * res_ang_y[m0] * exp_iAR0;
                    nlm[3][current_idx] += radial_factor * res_ang_z[m0] * exp_iAR0;
                }
                current_idx++;
            }

        } // End Radial Loop

        index_offset += num_m0;
    } // End Projector Loop

    // 6. Final Conjugation
    // Apply conjugation to all elements as per convention <phi|beta> = <beta|phi>*
    for (int dim = 0; dim < nlm.size(); dim++)
    {
        for (auto& x: nlm[dim])
        {
            x = std::conj(x);
        }
    }

    assert(index_offset == natomwfc);
    ModuleBase::timer::tick("module_rt", "snap_psibeta_half_tddft");
}

} // namespace module_rt