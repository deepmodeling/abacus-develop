#include "snap_psibeta_half_tddft.h"

#include "source_base/constants.h"
#include "source_base/math_integral.h"
#include "source_base/math_polyint.h"
#include "source_base/timer.h"
#include "source_base/ylm.h"

namespace module_rt
{

// nlm[0] : <phi|exp^{-iAr}|beta>
// nlm[1, 2, 3,] : <phi|r_a * exp^{-iAr}|beta>, which a = x, y, z.
void snap_psibeta_half_tddft(const LCAO_Orbitals& orb,
                             const InfoNonlocal& infoNL_,
                             std::vector<std::vector<std::complex<double>>>& nlm,
                             const ModuleBase::Vector3<double>& R1,
                             const int& T1,
                             const int& L1,
                             const int& m1,
                             const int& N1,
                             const ModuleBase::Vector3<double>& R0, // The projector.
                             const int& T0,
                             const ModuleBase::Vector3<double>& A,
                             const bool& calc_r)
{
    ModuleBase::timer::tick("module_rt", "snap_psibeta_half_tddft");

    // find number of projectors on atom R0
    const int nproj = infoNL_.nproj[T0];
    if (nproj == 0)
    {
        if (calc_r)
        {
            nlm.resize(4);
        }
        else
        {
            nlm.resize(1);
        }
        return;
    }

    std::vector<bool> calproj(nproj);
    std::vector<int> rmesh1(nproj);

    if (calc_r)
    {
        nlm.resize(4);
    }
    else
    {
        nlm.resize(1);
    }

    // Count number of projectors (l,m)
    int natomwfc = 0;
    for (int ip = 0; ip < nproj; ip++)
    {
        //============================
        // Use pseudo-atomic orbitals
        //============================

        const int L0 = infoNL_.Beta[T0].Proj[ip].getL(); // mohan add 2021-05-07
        natomwfc += 2 * L0 + 1;
    }

    for (int dim = 0; dim < nlm.size(); dim++)
    {
        nlm[dim].resize(natomwfc);
        for (auto& x: nlm[dim])
        {
            x = 0.0;
        }
    }

    // rcut of orbtials and projectors
    // in our calculation, we always put orbital phi at the left side of <phi|beta>
    // because <phi|beta> = <beta|phi>
    const double Rcut1 = orb.Phi[T1].getRcut();
    const ModuleBase::Vector3<double> dRa = R0 - R1;

    double distance10 = dRa.norm();

    bool all_out = true;
    for (int ip = 0; ip < nproj; ip++)
    {
        const double Rcut0 = infoNL_.Beta[T0].Proj[ip].getRcut();
        if (distance10 > (Rcut1 + Rcut0))
        {
            calproj[ip] = false;
        }
        else
        {
            all_out = false;
            calproj[ip] = true;
        }
    }

    if (all_out)
    {
        ModuleBase::timer::tick("module_rt", "snap_psibeta_half_tddft");
        return;
    }

    const int mesh_r1 = orb.Phi[T1].PhiLN(L1, N1).getNr();
    const double* psi_1 = orb.Phi[T1].PhiLN(L1, N1).getPsi();
    const double dk_1 = orb.Phi[T1].PhiLN(L1, N1).getDk();
    const double inv_dk_1 = 1.0 / dk_1;

    const int ridial_grid_num = 140;
    const int angular_grid_num = 110;
    std::vector<double> r_ridial(ridial_grid_num);
    std::vector<double> weights_ridial(ridial_grid_num);

    // OPTIMIZATION START: Cache standard Gauss-Legendre grid
    static std::vector<double> gl_x(ridial_grid_num);
    static std::vector<double> gl_w(ridial_grid_num);
    static bool gl_init = false;

    // Thread-safe initialization
    if (!gl_init)
    {
        #pragma omp critical(init_gauss_legendre)
        {
            if (!gl_init)
            {
                ModuleBase::Integral::Gauss_Legendre_grid_and_weight(ridial_grid_num, gl_x.data(), gl_w.data());
                gl_init = true;
            }
        }
    }
    // OPTIMIZATION END

    // Precompute A dot r_angular for Lebedev grid
    std::vector<double> A_dot_lebedev(angular_grid_num);
    for (int ian = 0; ian < angular_grid_num; ++ian)
    {
        A_dot_lebedev[ian] = A.x * ModuleBase::Integral::Lebedev_Laikov_grid110_x[ian] +
                             A.y * ModuleBase::Integral::Lebedev_Laikov_grid110_y[ian] +
                             A.z * ModuleBase::Integral::Lebedev_Laikov_grid110_z[ian];
    }

    // Buffers to reuse
    std::vector<std::complex<double>> result_angular;
    std::vector<std::complex<double>> result_angular_r_commu_x;
    std::vector<std::complex<double>> result_angular_r_commu_y;
    std::vector<std::complex<double>> result_angular_r_commu_z;
    std::vector<double> rly1;
    std::vector<std::vector<double>> rly0_all(angular_grid_num);

    int index = 0;
    for (int nb = 0; nb < nproj; nb++)
    {
        const int L0 = infoNL_.Beta[T0].Proj[nb].getL();
        if (!calproj[nb])
        {
            index += 2 * L0 + 1;
            continue;
        }

        const int mesh_r0 = infoNL_.Beta[T0].Proj[nb].getNr();
        const double* beta_r = infoNL_.Beta[T0].Proj[nb].getBeta_r();
        const double* radial0 = infoNL_.Beta[T0].Proj[nb].getRadial();
        const double dk_0 = infoNL_.Beta[T0].Proj[nb].getDk();

        const double Rcut0 = infoNL_.Beta[T0].Proj[nb].getRcut();

        // OPTIMIZATION: Use precomputed standard Gauss-Legendre grid
        double r_min = radial0[0];
        double r_max = radial0[mesh_r0 - 1];
        double xl = (r_max - r_min) * 0.5;
        double xmean = (r_max + r_min) * 0.5;
        
        for(int i=0; i<ridial_grid_num; ++i)
        {
            r_ridial[i] = xmean + xl * gl_x[i];
            weights_ridial[i] = xl * gl_w[i];
        }

        const double A_phase = A * R0;
        const std::complex<double> exp_iAR0 = std::exp(ModuleBase::IMAG_UNIT * A_phase);

        // Precompute rly0 for all angular points
        for(int ian = 0; ian < angular_grid_num; ++ian) {
            ModuleBase::Ylm::rl_sph_harm(L0, 
                ModuleBase::Integral::Lebedev_Laikov_grid110_x[ian], 
                ModuleBase::Integral::Lebedev_Laikov_grid110_y[ian], 
                ModuleBase::Integral::Lebedev_Laikov_grid110_z[ian], 
                rly0_all[ian]);
        }

        // Resize buffers
        if (result_angular.size() < 2 * L0 + 1)
        {
            result_angular.resize(2 * L0 + 1);
            if (calc_r)
            {
                result_angular_r_commu_x.resize(2 * L0 + 1);
                result_angular_r_commu_y.resize(2 * L0 + 1);
                result_angular_r_commu_z.resize(2 * L0 + 1);
            }
        }

        for (int ir = 0; ir < ridial_grid_num; ir++)
        {
            // Reset result accumulators
            std::fill(result_angular.begin(), result_angular.begin() + (2 * L0 + 1), 0.0);
            if (calc_r)
            {
                std::fill(result_angular_r_commu_x.begin(), result_angular_r_commu_x.begin() + (2 * L0 + 1), 0.0);
                std::fill(result_angular_r_commu_y.begin(), result_angular_r_commu_y.begin() + (2 * L0 + 1), 0.0);
                std::fill(result_angular_r_commu_z.begin(), result_angular_r_commu_z.begin() + (2 * L0 + 1), 0.0);
            }

            for (int ian = 0; ian < angular_grid_num; ian++)
            {
                const double x = ModuleBase::Integral::Lebedev_Laikov_grid110_x[ian];
                const double y = ModuleBase::Integral::Lebedev_Laikov_grid110_y[ian];
                const double z = ModuleBase::Integral::Lebedev_Laikov_grid110_z[ian];
                const double weights_angular = ModuleBase::Integral::Lebedev_Laikov_grid110_w[ian];
                
                const double r_val = r_ridial[ir];
                const ModuleBase::Vector3<double> r_coor(r_val * x, r_val * y, r_val * z);
                
                const ModuleBase::Vector3<double> tmp_r_coor = r_coor + dRa;
                const double tmp_r_coor_norm = tmp_r_coor.norm();
                
                if (tmp_r_coor_norm > Rcut1) 
                {
                    continue;
                }

                ModuleBase::Vector3<double> tmp_r_unit;
                if (tmp_r_coor_norm > 1e-10)
                {
                    tmp_r_unit = tmp_r_coor / tmp_r_coor_norm;
                }

                const std::vector<double>& rly0_vec = rly0_all[ian];

                ModuleBase::Ylm::rl_sph_harm(L1, tmp_r_unit.x, tmp_r_unit.y, tmp_r_unit.z, rly1);

                const double phase = r_val * A_dot_lebedev[ian];
                const std::complex<double> exp_iAr = std::exp(ModuleBase::IMAG_UNIT * phase);

                const ModuleBase::Vector3<double> tmp_r_coor_r_commu = r_coor + R0;
                
                // OPTIMIZATION: Inline Polynomial Interpolation
                double position = tmp_r_coor_norm * inv_dk_1;
                int iq = static_cast<int>(position);
                double interp_v = 0.0;
                
                if (iq <= mesh_r1 - 4)
                {
                    const double x0 = position - static_cast<double>(iq);
                    const double x1 = 1.0 - x0;
                    const double x2 = 2.0 - x0;
                    const double x3 = 3.0 - x0;
                    interp_v = x1*x2*(psi_1[iq]*x3+psi_1[iq+3]*x0)/6.0
                             + x0*x3*(psi_1[iq+1]*x2-psi_1[iq+2]*x1)/2.0;
                }
                
                const double weight_interp = interp_v * weights_angular;
                const int offset_L0 = L0 * L0;
                const int offset_L1 = L1 * L1 + m1;
                const double rly1_val = rly1[offset_L1];
                
                const std::complex<double> common_factor = exp_iAr * rly1_val * weight_interp;

                for (int m0 = 0; m0 < 2 * L0 + 1; m0++)
                {
                    std::complex<double> temp = common_factor * rly0_vec[offset_L0 + m0];
                    result_angular[m0] += temp;

                    if (calc_r)
                    {
                        result_angular_r_commu_x[m0] += temp * tmp_r_coor_r_commu.x;
                        result_angular_r_commu_y[m0] += temp * tmp_r_coor_r_commu.y;
                        result_angular_r_commu_z[m0] += temp * tmp_r_coor_r_commu.z;
                    }
                }
            }

            int index_tmp = index;
            const double temp = ModuleBase::PolyInt::Polynomial_Interpolation(beta_r,
                     mesh_r0, dk_0, r_ridial[ir]) * r_ridial[ir] * weights_ridial[ir];

            if (!calc_r)
            {
                for (int m0 = 0; m0 < 2 * L0 + 1; m0++)
                {
                    nlm[0][index_tmp] += temp * result_angular[m0] * exp_iAR0;
                    index_tmp++;
                }
            }
            else
            {
                for (int m0 = 0; m0 < 2 * L0 + 1; m0++)
                {
                    nlm[0][index_tmp] += temp * result_angular[m0] * exp_iAR0;
                    nlm[1][index_tmp] += temp * result_angular_r_commu_x[m0] * exp_iAR0;
                    nlm[2][index_tmp] += temp * result_angular_r_commu_y[m0] * exp_iAR0;
                    nlm[3][index_tmp] += temp * result_angular_r_commu_z[m0] * exp_iAR0;
                    index_tmp++;
                }
            }
        }

        index += 2 * L0 + 1;
    }

    for(int dim = 0; dim < nlm.size(); dim++)
    {
        for (auto &x : nlm[dim])
        {
            // nlm[0] is <phi|exp^{-iAr}|beta>
            // nlm[1 or 2 or 3] is <phi|r_a * exp^{-iAr}|beta>, a = x, y, z
            x = std::conj(x); 
        }
    }

    assert(index == natomwfc);
    ModuleBase::timer::tick("module_rt", "snap_psibeta_half_tddft");

    return;
}

} // namespace module_rt
