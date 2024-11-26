#pragma once
#include "pulay_force_stress.h"
#include "module_base/timer.h"
#include "module_parameter/parameter.h"
namespace PulayForceStress
{
    // common kernel
    template <typename TK, typename TR, typename Tfunc>
    inline void cal_pulay_fs(
        ModuleBase::matrix& f,
        ModuleBase::matrix& s,
        const elecstate::DensityMatrix<TK, TR>& dm,
        const UnitCell& ucell,
        const Parallel_Orbitals& pv,
        const double** dHSx,
        const double** dHSxy,
        const double* dtau,
        const bool& isforce,
        const bool& isstress,
        Record_adj* ra,
        const double& factor_force,
        const double& factor_stress,
        Tfunc& stress_func)
    {
        ModuleBase::TITLE("Force_LCAO", "cal_pulay_fs_center2");
        ModuleBase::timer::tick("Force_LCAO", "cal_pulay_fs_center2");

        const int nspin_DMR = (PARAM.inp.nspin == 2) ? 2 : 1;
        int total_irr = 0;
#ifdef _OPENMP
#pragma omp parallel
        {
            int num_threads = omp_get_num_threads();
            ModuleBase::matrix local_s(3, 3);
            int local_total_irr = 0;
#else
            ModuleBase::matrix& local_s = s;
            int& local_total_irr = total_irr;
#endif

#ifdef _OPENMP
#pragma omp for schedule(dynamic)
#endif
            for (int iat = 0; iat < ucell.nat; iat++)
            {
                const int T1 = ucell.iat2it[iat];
                Atom* atom1 = &ucell.atoms[T1];
                const int I1 = ucell.iat2ia[iat];
                // get iat1
                int iat1 = ucell.itia2iat(T1, I1);
                double* f_iat;
                if (isforce) { f_iat = &f(iat, 0); }
#ifdef _OPENMP
                // using local stack to avoid false sharing in multi-threaded case
                double f_tmp[3] = { 0.0, 0.0, 0.0 };
                if (num_threads > 1) { f_iat = f_tmp; }
#endif
                int irr = pv.nlocstart[iat];
                const int start1 = ucell.itiaiw2iwt(T1, I1, 0);
                for (int cb = 0; cb < ra->na_each[iat]; ++cb)
                {
                    const int T2 = ra->info[iat][cb][3];
                    const int I2 = ra->info[iat][cb][4];
                    const int start2 = ucell.itiaiw2iwt(T2, I2, 0);
                    Atom* atom2 = &ucell.atoms[T2];
                    // get iat2
                    int iat2 = ucell.itia2iat(T2, I2);
                    double Rx = ra->info[iat][cb][0];
                    double Ry = ra->info[iat][cb][1];
                    double Rz = ra->info[iat][cb][2];
                    // get BaseMatrix
                    if (pv.get_row_size(iat1) <= 0 || pv.get_col_size(iat2) <= 0) { continue; }
                    std::vector<hamilt::BaseMatrix<double>*> tmp_matrix;
                    for (int is = 0; is < nspin_DMR; ++is)
                    {
                        tmp_matrix.push_back(dm.get_DMR_pointer(is + 1)->find_matrix(iat1, iat2, Rx, Ry, Rz));
                    }
                    for (int mu = 0; mu < pv.get_row_size(iat1); ++mu)
                    {
                        for (int nu = 0; nu < pv.get_col_size(iat2); ++nu)
                        {
                            // the DMR should not be summed over spin, do the summation here
                            double dm2d1 = 0.0;
                            for (int is = 0; is < nspin_DMR; ++is) { dm2d1 += tmp_matrix[is]->get_value(mu, nu); }
                            double dm2d2 = 2.0 * dm2d1;
                            if (isforce)
                            {
                                const double dm2d2_f = dm2d2 * factor_force;
                                for (int i = 0; i < 3; ++i) { f_iat[i] += dm2d2_f * dHSx[i][irr]; }
                            }
                            if (isstress)
                            {
                                const double dm2d1_s = dm2d1 * factor_stress;
                                stress_func(local_s, dm2d1_s, dHSx, dHSxy, dtau, irr);
                            }
                            ++local_total_irr;
                            ++irr;
                        }
                    }
                }
#ifdef _OPENMP
                if (isforce && num_threads > 1) { for (int i = 0; i < 3; ++i) { f(iat, i) += f_iat[i]; } }
#endif
            } // end iat
#ifdef _OPENMP
#pragma omp critical(cal_foverlap_k_reduce)
            {
                total_irr += local_total_irr;
                if (isstress)
                {
                    for (int i = 0; i < 3; ++i) { for (int j = i; j < 3; ++j) { s(i, j) += local_s(i, j); } }
                }
            }
        }
#endif

        if (total_irr != pv.nnr)
        {
            ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running, "wrong irr", total_irr);
            ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running, "wrong pv.nnr", pv.nnr);
            ModuleBase::WARNING_QUIT("Force_LCAO::cal_pulay_fs_center2", "irr!=pv.nnr");
        }

        if (isstress) { StressTools::stress_fill(ucell.lat0, ucell.omega, s); }

        ModuleBase::timer::tick("Force_LCAO", "cal_pulay_fs_center2");
    }
}
