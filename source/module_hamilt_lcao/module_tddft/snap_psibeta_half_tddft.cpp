#include "snap_psibeta_half_tddft.h"

#include "module_base/constants.h"
#include "module_base/math_integral.h"
#include "module_base/math_polyint.h"
#include "module_base/timer.h"
#include "module_base/ylm.h"

namespace module_tddft
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
    ModuleBase::timer::tick("module_tddft", "snap_psibeta_half_tddft");

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
        ModuleBase::timer::tick("module_tddft", "snap_psibeta_half_tddft");
        return;
    }

    const int mesh_r1 = orb.Phi[T1].PhiLN(L1, N1).getNr();
    const double* psi_1 = orb.Phi[T1].PhiLN(L1, N1).getPsi();
    const double dk_1 = orb.Phi[T1].PhiLN(L1, N1).getDk();

    const int ridial_grid_num = 140;
    const int angular_grid_num = 110;
    std::vector<double> r_ridial(ridial_grid_num);
    std::vector<double> weights_ridial(ridial_grid_num);

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
        ModuleBase::Integral::Gauss_Legendre_grid_and_weight(radial0[0],
                                                             radial0[mesh_r0 - 1],
                                                             ridial_grid_num,
                                                             r_ridial.data(),
                                                             weights_ridial.data());

        const double A_phase = A * R0;
        const std::complex<double> exp_iAR0 = std::exp(ModuleBase::IMAG_UNIT * A_phase);

        std::vector<double> rly0(L0);
        std::vector<double> rly1(L1);
        for (int ir = 0; ir < ridial_grid_num; ir++)
        {
            std::vector<std::complex<double>> result_angular(2 * L0 + 1, 0.0);
            std::vector<std::complex<double>> result_angular_r_commu_x;
            std::vector<std::complex<double>> result_angular_r_commu_y;
            std::vector<std::complex<double>> result_angular_r_commu_z;
            if (calc_r)
            {
                result_angular_r_commu_x.resize(2 * L0 + 1, 0.0);
                result_angular_r_commu_y.resize(2 * L0 + 1, 0.0);
                result_angular_r_commu_z.resize(2 * L0 + 1, 0.0);
            }

            for (int ian = 0; ian < angular_grid_num; ian++)
            {
                const double x = ModuleBase::Integral::Lebedev_Laikov_grid110_x[ian];
                const double y = ModuleBase::Integral::Lebedev_Laikov_grid110_y[ian];
                const double z = ModuleBase::Integral::Lebedev_Laikov_grid110_z[ian];
                const double weights_angular = ModuleBase::Integral::Lebedev_Laikov_grid110_w[ian];
                const ModuleBase::Vector3<double> r_angular_tmp(x, y, z);

                const ModuleBase::Vector3<double> r_coor = r_ridial[ir] * r_angular_tmp;
                const ModuleBase::Vector3<double> tmp_r_coor = r_coor + dRa;
                const double tmp_r_coor_norm = tmp_r_coor.norm();
                if (tmp_r_coor_norm > Rcut1) {
                    continue;
                }

                ModuleBase::Vector3<double> tmp_r_unit;
                if (tmp_r_coor_norm > 1e-10)
                {
                    tmp_r_unit = tmp_r_coor / tmp_r_coor_norm;
                }

                ModuleBase::Ylm::rl_sph_harm(L0, x, y, z, rly0);

                ModuleBase::Ylm::rl_sph_harm(L1, tmp_r_unit.x, tmp_r_unit.y, tmp_r_unit.z, rly1);

                const double phase = A * r_coor;
                const std::complex<double> exp_iAr = std::exp(ModuleBase::IMAG_UNIT * phase);

                const ModuleBase::Vector3<double> tmp_r_coor_r_commu = r_coor + R0;
                const double temp_interpolation_value = ModuleBase::PolyInt::Polynomial_Interpolation(psi_1, mesh_r1, dk_1, tmp_r_coor_norm);
                for (int m0 = 0; m0 < 2 * L0 + 1; m0++)
                {
                    std::complex<double> temp = exp_iAr * rly0[L0 * L0 + m0] * rly1[L1 * L1 + m1]
                                                * temp_interpolation_value
                                                * weights_angular;
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
            const double temp = ModuleBase::PolyInt::Polynomial_Interpolation(beta_r, mesh_r0, dk_0, r_ridial[ir]) * r_ridial[ir] * weights_ridial[ir];
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
    ModuleBase::timer::tick("module_tddft", "snap_psibeta_half_tddft");

    return;
}

} // namespace module_tddft
