#include "snap_psibeta_half_tddft.h"

#include "source_base/array_pool.h"
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
 *        and its derivatives (if calc_r and/or calc_deri are true).
 *
 *        This function integrates the overlap between a local orbital phi (at R1)
 *        and a non-local projector beta (at R0), modulated by a plane-wave-like phase factor
 *        exp^{-iAr}, where A is a vector potential.
 * IMPORTANT MATHEMATICAL IDENTITY FOR DERIVATIVES:
 * ================================================
 * For the integral M_ij = <phi_i | e^{iA·r} | beta_j>, the gradients with respect
 * to atomic centers R_i (phi center) and R_j (beta center) satisfy:
 *
 *     ∇_{R_i} M_ij + ∇_{R_j} M_ij = i·A·M_ij
 *
 * Or equivalently (solving for one gradient in terms of the other):
 *
 *     ∇_{R_i} M_ij = -∇_{R_j} M_ij + i·A·M_ij
 *
 * This means:
 * - You CANNOT simply negate ∇_{R_j} to get ∇_{R_i}
 * - You MUST add a correction term proportional to the vector potential A
 *   and the integral value M_ij itself
 * - The correction term i·A·M_ij arises from the phase factor e^{iA·r}
 *
 * Physical Interpretation:
 * - In the absence of vector potential (A=0), the gradients are simply negatives
 * - With vector potential, the phase factor breaks this simple symmetry
 * - The correction ensures gauge invariance and proper electromagnetic coupling
 *
 * @param orb LCAO Orbitals information
 * @param infoNL_ Non-local pseudopotential information
 * @param nlm Output (size depends on calc_r and calc_deri):
 *            calc_r=false, calc_deri=false (size=1):
 *              nlm[0] : <phi|exp^{-iAr}|beta>
 *            calc_r=true, calc_deri=false (size=4):
 *              nlm[0] : <phi|exp^{-iAr}|beta>
 *              nlm[1-3] : <phi|r_a * exp^{-iAr}|beta>, a = x, y, z
 *            calc_r=false, calc_deri=true (size=4):
 *              nlm[0] : <phi|exp^{-iAr}|beta>
 *              nlm[1-3] : <phi|exp^{-iAr}|∂beta/∂τ_a>, a = x, y, z
 *            calc_r=true, calc_deri=true (size=16):
 *              nlm[0] : <phi|exp^{-iAr}|beta>
 *              nlm[1-3] : <phi|r_a * exp^{-iAr}|beta>, a = x, y, z
 *              nlm[4-6] : <phi|exp^{-iAr}|∂beta/∂τ_a>, a = x, y, z
 *              nlm[7-15] : <phi|r_a * exp^{-iAr}|∂beta/∂τ_b>, a,b = x,y,z (3x3 tensor)
 *                         Row-major ordering: [r_x*∂/∂τ_x, r_x*∂/∂τ_y, r_x*∂/∂τ_z,
 *                                              r_y*∂/∂τ_x, r_y*∂/∂τ_y, r_y*∂/∂τ_z,
 *                                              r_z*∂/∂τ_x, r_z*∂/∂τ_y, r_z*∂/∂τ_z]
 * @param R1 Position of atom 1 (orbital phi)
 * @param T1 Type of atom 1
 * @param L1 Angular momentum of orbital phi
 * @param m1 Magnetic quantum number of orbital phi
 * @param N1 Radial quantum number of orbital phi
 * @param R0 Position of atom 0 (projector beta)
 * @param T0 Type of atom 0
 * @param A Vector potential A (or related field vector)
 * @param calc_r Whether to calculate position operator matrix elements
 * @param calc_deri Whether to calculate derivative terms (∂beta/∂τ)
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
                             const bool& calc_r,
                             const bool& calc_deri)
{
    ModuleBase::timer::tick("module_rt", "snap_psibeta_half_tddft");

    // 1. Initialization and Early Exits
    const int nproj = infoNL_.nproj[T0];

    // Determine required output size based on flags
    int required_size = 1;  // Base: overlap
    if (calc_r && !calc_deri)
        required_size = 4;  // overlap + position operators (x,y,z)
    else if (!calc_r && calc_deri)
        required_size = 4;  // overlap + derivatives (dx,dy,dz)
    else if (calc_r && calc_deri)
        required_size = 16; // overlap(1) + position(3) + derivatives(3) + tensor(9)

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
    // Accumulators for derivative components
    std::vector<std::complex<double>> res_ang_dx, res_ang_dy, res_ang_dz;
    // Accumulators for full 3x3 position*derivative tensor
    // Index mapping: [0]=r_x*∂/∂τ_x, [1]=r_x*∂/∂τ_y, [2]=r_x*∂/∂τ_z,
    //                [3]=r_y*∂/∂τ_x, [4]=r_y*∂/∂τ_y, [5]=r_y*∂/∂τ_z,
    //                [6]=r_z*∂/∂τ_x, [7]=r_z*∂/∂τ_y, [8]=r_z*∂/∂τ_z
    std::vector<std::complex<double>> res_ang_tensor[9];
    // Accumulators for radial contribution to tensor
    std::vector<std::complex<double>> res_ang_tensor_radial[9];

    std::vector<double> rly1((L1 + 1) * (L1 + 1));                 // Spherical harmonics buffer for L1
    std::vector<std::vector<double>> rly0_cache(angular_grid_num); // Cache for L0 Ylm
    ModuleBase::Array_Pool<double> grly0_cache; // Cache for L0 Ylm gradients (will be initialized per projector)

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

        // Precompute spherical harmonics gradients if derivatives are needed
        if (calc_deri)
        {
            // Initialize Array_Pool for current L0
            const int ylm_size = (L0 + 1) * (L0 + 1);
            grly0_cache = ModuleBase::Array_Pool<double>(ylm_size * angular_grid_num, 3);

            std::vector<double> rly0_temp(ylm_size);
            for (int ian = 0; ian < angular_grid_num; ++ian)
            {
                ModuleBase::Ylm::grad_rl_sph_harm(L0,
                                                  ModuleBase::Integral::Lebedev_Laikov_grid110_x[ian],
                                                  ModuleBase::Integral::Lebedev_Laikov_grid110_y[ian],
                                                  ModuleBase::Integral::Lebedev_Laikov_grid110_z[ian],
                                                  rly0_temp.data(),
                                                  grly0_cache.get_ptr_2D() + ian * ylm_size);
            }
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
            if (calc_deri)
            {
                res_ang_dx.resize(num_m0);
                res_ang_dy.resize(num_m0);
                res_ang_dz.resize(num_m0);
                if (calc_r)
                {
                    for (int i = 0; i < 9; i++)
                    {
                        res_ang_tensor[i].resize(num_m0);
                        res_ang_tensor_radial[i].resize(num_m0);
                    }
                }
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
            if (calc_deri)
            {
                std::fill(res_ang_dx.begin(), res_ang_dx.begin() + num_m0, 0.0);
                std::fill(res_ang_dy.begin(), res_ang_dy.begin() + num_m0, 0.0);
                std::fill(res_ang_dz.begin(), res_ang_dz.begin() + num_m0, 0.0);
                if (calc_r)
                {
                    for (int i = 0; i < 9; i++)
                    {
                        std::fill(res_ang_tensor[i].begin(), res_ang_tensor[i].begin() + num_m0, 0.0);
                        std::fill(res_ang_tensor_radial[i].begin(), res_ang_tensor_radial[i].begin() + num_m0, 0.0);
                    }
                }
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

                    if (calc_deri)
                    {
                        // Get gradient of spherical harmonics from cache
                        const int ylm_idx = ian * ((L0 + 1) * (L0 + 1)) + offset_L0 + m0;
                        const double grad_ylm_x = grly0_cache[ylm_idx][0];
                        const double grad_ylm_y = grly0_cache[ylm_idx][1];
                        const double grad_ylm_z = grly0_cache[ylm_idx][2];

                        // Accumulate angular part: beta * (∂Y_lm/∂τ_a)
                        // Note: Radial part (∂beta/∂r)*(r_a/r) will be added in radial loop
                        res_ang_dx[m0] += common_factor * grad_ylm_x;
                        res_ang_dy[m0] += common_factor * grad_ylm_y;
                        res_ang_dz[m0] += common_factor * grad_ylm_z;

                        if (calc_r)
                        {
                            // Compute full 3x3 tensor: r_a * ∂/∂τ_b
                            double r_op_x = rx + R0.x;
                            double r_op_y = ry + R0.y;
                            double r_op_z = rz + R0.z;

                            // Store gradient components in array for cleaner indexing
                            const double grad_ylm[3] = {grad_ylm_x, grad_ylm_y, grad_ylm_z};
                            const double r_op[3] = {r_op_x, r_op_y, r_op_z};

                            // Compute all 9 tensor components: r_op[a] * grad_ylm[b]
                            // This is the angular contribution to the tensor
                            for (int a = 0; a < 3; a++)
                            {
                                for (int b = 0; b < 3; b++)
                                {
                                    res_ang_tensor[a * 3 + b][m0] += common_factor * r_op[a] * grad_ylm[b];
                                }
                            }

                            // Store angular integrals for radial derivative contribution
                            // These will be multiplied by dbeta_dr later
                            // For radial contribution: r_op_a * (r_b/r)
                            const double r_unit[3] = {x, y, z};  // Unit vector

                            for (int a = 0; a < 3; a++)
                            {
                                for (int b = 0; b < 3; b++)
                                {
                                    res_ang_tensor_radial[a * 3 + b][m0] += common_factor * r_op[a] * r_unit[b];
                                }
                            }
                        }
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

            // Compute beta radial derivative if needed
            double dbeta_dr = 0.0;
            if (calc_deri && r_radial[ir] > 1e-10)
            {
                // Adaptive step size: smaller of 1e-6 or 0.1% of r
                const double dr = std::min(1e-6, r_radial[ir] * 1e-3);
                double beta_plus = ModuleBase::PolyInt::Polynomial_Interpolation(beta_r, mesh_r0, dk_0, r_radial[ir] + dr);
                double beta_minus = ModuleBase::PolyInt::Polynomial_Interpolation(beta_r, mesh_r0, dk_0, r_radial[ir] - dr);
                dbeta_dr = (beta_plus - beta_minus) / (2.0 * dr);
            }

            double radial_factor = beta_val * r_radial[ir] * w_radial[ir];

            // Compute radial derivative contribution: (∂beta/∂r) * (r_a/r)
            // This needs to be integrated over angular directions
            std::vector<std::complex<double>> radial_contrib_x(num_m0, 0.0);
            std::vector<std::complex<double>> radial_contrib_y(num_m0, 0.0);
            std::vector<std::complex<double>> radial_contrib_z(num_m0, 0.0);

            if (calc_deri && r_radial[ir] > 1e-10)
            {
                const double inv_r = 1.0 / r_radial[ir];

                // Integrate radial contribution over angular grid
                for (int ian = 0; ian < angular_grid_num; ian++)
                {
                    const double x = ModuleBase::Integral::Lebedev_Laikov_grid110_x[ian];
                    const double y = ModuleBase::Integral::Lebedev_Laikov_grid110_y[ian];
                    const double z = ModuleBase::Integral::Lebedev_Laikov_grid110_z[ian];
                    const double w_ang = ModuleBase::Integral::Lebedev_Laikov_grid110_w[ian];

                    double rx = r_radial[ir] * x;
                    double ry = r_radial[ir] * y;
                    double rz = r_radial[ir] * z;

                    double tx = rx + dRa.x;
                    double ty = ry + dRa.y;
                    double tz = rz + dRa.z;
                    double tnorm = std::sqrt(tx * tx + ty * ty + tz * tz);

                    if (tnorm > Rcut1)
                        continue;

                    // Compute Ylm for L1
                    if (tnorm > 1e-10)
                    {
                        double inv_tnorm = 1.0 / tnorm;
                        ModuleBase::Ylm::rl_sph_harm(L1, tx * inv_tnorm, ty * inv_tnorm, tz * inv_tnorm, rly1);
                    }
                    else
                    {
                        ModuleBase::Ylm::rl_sph_harm(L1, 0.0, 0.0, 1.0, rly1);
                    }

                    const double phase = r_radial[ir] * A_dot_lebedev[ian];
                    const std::complex<double> exp_iAr = std::exp(ModuleBase::IMAG_UNIT * phase);
                    double interp_psi = ModuleBase::PolyInt::Polynomial_Interpolation(psi_1, mesh_r1, dk_1, tnorm);

                    const int offset_L1 = L1 * L1 + m1;
                    const double ylm_L1_val = rly1[offset_L1];
                    const std::complex<double> common_factor = exp_iAr * ylm_L1_val * interp_psi * w_ang;

                    const std::vector<double>& rly0_vec = rly0_cache[ian];
                    const int offset_L0 = L0 * L0;

                    for (int m0 = 0; m0 < num_m0; m0++)
                    {
                        std::complex<double> term = common_factor * rly0_vec[offset_L0 + m0];

                        // Add radial contribution: (∂beta/∂r) * (r_a/r)
                        radial_contrib_x[m0] += term * (rx * inv_r);
                        radial_contrib_y[m0] += term * (ry * inv_r);
                        radial_contrib_z[m0] += term * (rz * inv_r);
                    }
                }
            }

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

                if (calc_deri)
                {
                    int deri_offset = calc_r ? 4 : 1;

                    // Accumulate complete derivative: angular part + radial part
                    // ∂beta/∂τ_a = beta * (∂Y_lm/∂τ_a) + (∂beta/∂r) * (r_a/r)
                    double radial_factor_beta = beta_val * r_radial[ir] * w_radial[ir];
                    double radial_factor_deriv = dbeta_dr * r_radial[ir] * w_radial[ir];

                    nlm[deri_offset][current_idx] += (radial_factor_beta * res_ang_dx[m0]
                                                      + radial_factor_deriv * radial_contrib_x[m0]) * exp_iAR0;
                    nlm[deri_offset + 1][current_idx] += (radial_factor_beta * res_ang_dy[m0]
                                                          + radial_factor_deriv * radial_contrib_y[m0]) * exp_iAR0;
                    nlm[deri_offset + 2][current_idx] += (radial_factor_beta * res_ang_dz[m0]
                                                          + radial_factor_deriv * radial_contrib_z[m0]) * exp_iAR0;

                    if (calc_r)
                    {
                        // Accumulate full 3x3 tensor with both angular and radial contributions
                        // Angular contribution (beta * ∂Y_lm/∂τ_b)
                        for (int i = 0; i < 9; i++)
                        {
                            nlm[7 + i][current_idx] += radial_factor_beta * res_ang_tensor[i][m0] * exp_iAR0;
                        }

                        // Radial contribution ((∂beta/∂r) * (r_b/r))
                        for (int i = 0; i < 9; i++)
                        {
                            nlm[7 + i][current_idx] += radial_factor_deriv * res_ang_tensor_radial[i][m0] * exp_iAR0;
                        }
                    }
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