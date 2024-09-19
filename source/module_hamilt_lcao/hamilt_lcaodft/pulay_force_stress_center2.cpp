#include "pulay_force_stress.h"
namespace PulayForceStress
{
    template<>  // gamma-only, provided xy
    void cal_pulay_fs(
        ModuleBase::matrix& f,
        ModuleBase::matrix& s,
        const elecstate::DensityMatrix<double, double>& dm,
        const UnitCell& ucell,
        const Parallel_Orbitals& pv,
        const double* (&dHSx)[3],
        const double* (&dHSxy)[6],
        const bool& isforce,
        const bool& isstress,
        Record_adj* ra,
        const double& factor_force,
        const double& factor_stress)
    {
        ModuleBase::TITLE("Force_LCAO", "cal_pulay_fs_center2");
        ModuleBase::timer::tick("Force_LCAO", "cal_pulay_fs_center2");

        const int nspin = PARAM.inp.nspin;
        const int nlocal = PARAM.globalv.nlocal;

        for (int i = 0; i < nlocal; ++i)
        {
            const int iat = ucell.iwt2iat[i];
            for (int j = 0; j < nlocal; ++j)
            {
                const int mu = pv.global2local_row(j);
                const int nu = pv.global2local_col(i);

                if (mu >= 0 && nu >= 0)
                {
                    const int index = mu * pv.ncol + nu;
                    double sum = 0.0;
                    for (int is = 0; is < nspin; ++is) { sum += dm.get_DMK(is + 1, 0, nu, mu); }
                    if (isforce)
                    {
                        const double sumf = sum * factor_force;
                        for (int i = 0; i < 3; ++i) { f(iat, i) += sumf * 2.0 * dHSx[i][index]; }
                    }
                    if (isstress)
                    {
                        const double sums = sum * factor_stress;
                        int ij = 0;
                        for (int i = 0; i < 3;++i) { for (int j = i; j < 3; ++j) { s(i, j) += sums * dHSxy[ij++][index]; } }
                    }
                }
            }
        }

        if (isstress) { StressTools::stress_fill(ucell.lat0, ucell.omega, s); }

        ModuleBase::timer::tick("Force_LCAO", "cal_pulay_fs_center2");
    }

    template<>  //multi-k, provided xy
    void cal_pulay_fs(
        ModuleBase::matrix& f,
        ModuleBase::matrix& s,
        const elecstate::DensityMatrix<std::complex<double>, double>& dm,
        const UnitCell& ucell,
        const Parallel_Orbitals& pv,
        const double* (&dHSx)[3],
        const double* (&dHSxy)[6],
        const bool& isforce,
        const bool& isstress,
        Record_adj* ra,
        const double& factor_force,
        const double& factor_stress)
    {
        auto stress_func = [](ModuleBase::matrix& local_s, const double& dm2d1_s, const double** dHSx, const double** dHSxy, const double* dtau, const int& irr)
            {
                int ij = 0;
                for (int i = 0; i < 3; ++i) { for (int j = i; j < 3; ++j) { local_s(i, j) += dm2d1_s * dHSxy[ij++][irr]; } }
            };
        cal_pulay_fs(f, s, dm, ucell, pv, dHSx, dHSxy, nullptr, isforce, isstress, ra, factor_force, factor_stress, stress_func);
    }

    template<>  // multi-k, provided x
    void cal_pulay_fs(
        ModuleBase::matrix& f,
        ModuleBase::matrix& s,
        const elecstate::DensityMatrix<std::complex<double>, double>& dm,
        const UnitCell& ucell,
        const Parallel_Orbitals& pv,
        const double* (&dHSx)[3],
        const double* dtau,
        const bool& isforce,
        const bool& isstress,
        Record_adj* ra,
        const double& factor_force,
        const double& factor_stress)
    {
        auto stress_func = [](ModuleBase::matrix& local_s, const double& dm2d1_s, const double** dHSx, const double** dHSxy, const double* dtau, const int& irr)
            {
                for (int i = 0; i < 3; ++i) { for (int j = i; j < 3; ++j) { local_s(i, j) += dm2d1_s * dHSx[i][irr] * dtau[irr * 3 + j]; } }
            };
        cal_pulay_fs(f, s, dm, ucell, pv, dHSx, nullptr, dtau, isforce, isstress, ra, factor_force, factor_stress, stress_func);
    }

}