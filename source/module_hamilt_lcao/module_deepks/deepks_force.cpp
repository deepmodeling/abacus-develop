#include "module_parameter/parameter.h"

#ifdef __DEEPKS

#include "deepks_force.h"
#include "deepks_iterate.h"
#include "module_base/constants.h"
#include "module_base/libm/libm.h"
#include "module_base/timer.h"
#include "module_base/vector3.h"
#include "module_hamilt_lcao/module_hcontainer/atom_pair.h"

template <typename TK>
void DeePKS_domain::cal_f_delta(const std::vector<std::vector<TK>>& dm,
                                const UnitCell& ucell,
                                const LCAO_Orbitals& orb,
                                const Grid_Driver& GridD,
                                const Parallel_Orbitals& pv,
                                const int nks,
                                const std::vector<ModuleBase::Vector3<double>>& kvec_d,
                                std::vector<hamilt::HContainer<double>*> phialpha,
                                double** gedm,
                                ModuleBase::IntArray* inl_index,
                                ModuleBase::matrix& f_delta,
                                const bool isstress,
                                ModuleBase::matrix& svnl_dalpha)
{
    ModuleBase::TITLE("DeePKS_domain", "cal_f_delta");
    ModuleBase::timer::tick("DeePKS_domain", "cal_f_delta");
    f_delta.zero_out();

    const int lmaxd = orb.get_lmax_d();

    DeePKS_domain::iterate_ad2(
        ucell,
        GridD,
        orb,
        false, // no trace_alpha
        [&](const int iat,
            const ModuleBase::Vector3<double>& tau0,
            const int ibt1,
            const ModuleBase::Vector3<double>& tau1,
            const int start1,
            const int nw1_tot,
            ModuleBase::Vector3<int> dR1,
            const int ibt2,
            const ModuleBase::Vector3<double>& tau2,
            const int start2,
            const int nw2_tot,
            ModuleBase::Vector3<int> dR2)
        {
            double r1[3] = {0, 0, 0};
            double r2[3] = {0, 0, 0};
            if (isstress)
            {
                r1[0] = (tau1.x - tau0.x);
                r1[1] = (tau1.y - tau0.y);
                r1[2] = (tau1.z - tau0.z);
                r2[0] = (tau2.x - tau0.x);
                r2[1] = (tau2.y - tau0.y);
                r2[2] = (tau2.z - tau0.z);
            }

            auto row_indexes = pv.get_indexes_row(ibt1);
            auto col_indexes = pv.get_indexes_col(ibt2);

            if (row_indexes.size() * col_indexes.size() == 0)
            {
                return; // to next loop
            }

            int dRx = 0;
            int dRy = 0;
            int dRz = 0;
            if constexpr (std::is_same<TK, std::complex<double>>::value)
            {
                dRx = dR2.x - dR1.x;
                dRy = dR2.y - dR1.y;
                dRz = dR2.z - dR1.z;
            }
            ModuleBase::Vector3<double> dR(dRx, dRy, dRz);

            hamilt::AtomPair<double> dm_pair(ibt1, ibt2, dRx, dRy, dRz, &pv);

            dm_pair.allocate(nullptr, true);

            if constexpr (std::is_same<TK, double>::value) // for gamma-only
            {
                for (int is = 0; is < dm.size(); is++)
                {
                    if (ModuleBase::GlobalFunc::IS_COLUMN_MAJOR_KS_SOLVER(PARAM.inp.ks_solver))
                    {
                        dm_pair.add_from_matrix(dm[is].data(), pv.get_row_size(), 1.0, 1);
                    }
                    else
                    {
                        dm_pair.add_from_matrix(dm[is].data(), pv.get_col_size(), 1.0, 0);
                    }
                }
            }
            else // for multi-k
            {
                for (int ik = 0; ik < nks; ik++)
                {
                    const double arg = -(kvec_d[ik] * dR) * ModuleBase::TWO_PI;
                    double sinp, cosp;
                    ModuleBase::libm::sincos(arg, &sinp, &cosp);
                    const std::complex<double> kphase = std::complex<double>(cosp, sinp);
                    if (ModuleBase::GlobalFunc::IS_COLUMN_MAJOR_KS_SOLVER(PARAM.inp.ks_solver))
                    {
                        dm_pair.add_from_matrix(dm[ik].data(), pv.get_row_size(), kphase, 1);
                    }
                    else
                    {
                        dm_pair.add_from_matrix(dm[ik].data(), pv.get_col_size(), kphase, 0);
                    }
                }
            }

            hamilt::BaseMatrix<double>* overlap_1 = phialpha[0]->find_matrix(iat, ibt1, dR1);
            hamilt::BaseMatrix<double>* overlap_2 = phialpha[0]->find_matrix(iat, ibt2, dR2);
            if (overlap_1 == nullptr || overlap_2 == nullptr)
            {
                return; // to next loop
            }
            std::vector<hamilt::BaseMatrix<double>*> grad_overlap_1(3);
            std::vector<hamilt::BaseMatrix<double>*> grad_overlap_2(3);
            for (int i = 0; i < 3; ++i)
            {
                grad_overlap_1[i] = phialpha[i + 1]->find_matrix(iat, ibt1, dR1);
                grad_overlap_2[i] = phialpha[i + 1]->find_matrix(iat, ibt2, dR2);
            }

            assert(overlap_1->get_col_size() == overlap_2->get_col_size());

            const double* dm_current = dm_pair.get_pointer();

            for (int iw1 = 0; iw1 < row_indexes.size(); ++iw1)
            {
                for (int iw2 = 0; iw2 < col_indexes.size(); ++iw2)
                {
                    double nlm[3] = {0, 0, 0};
                    double nlm_t[3] = {0, 0, 0}; // for stress

                    if (!PARAM.inp.deepks_equiv)
                    {
                        int ib = 0;
                        for (int L0 = 0; L0 <= orb.Alpha[0].getLmax(); ++L0)
                        {
                            for (int N0 = 0; N0 < orb.Alpha[0].getNchi(L0); ++N0)
                            {
                                const int inl = inl_index[ucell.iat2it[iat]](ucell.iat2ia[iat], L0, N0);
                                const int nm = 2 * L0 + 1;
                                for (int m1 = 0; m1 < nm; ++m1)
                                {
                                    for (int m2 = 0; m2 < nm; ++m2)
                                    {
                                        for (int dim = 0; dim < 3; dim++)
                                        {
                                            nlm[dim] += gedm[inl][m1 * nm + m2]
                                                        * overlap_1->get_value(row_indexes[iw1], ib + m1)
                                                        * grad_overlap_2[dim]->get_value(col_indexes[iw2], ib + m2);
                                        }
                                    }
                                }
                                ib += nm;
                            }
                        }
                        assert(ib == overlap_1->get_col_size());
                    }
                    else
                    {
                        int nproj = 0;
                        for (int il = 0; il < lmaxd + 1; il++)
                        {
                            nproj += (2 * il + 1) * orb.Alpha[0].getNchi(il);
                        }
                        for (int iproj = 0; iproj < nproj; iproj++)
                        {
                            for (int jproj = 0; jproj < nproj; jproj++)
                            {
                                for (int dim = 0; dim < 3; dim++)
                                {
                                    nlm[dim] += gedm[iat][iproj * nproj + jproj]
                                                * overlap_1->get_value(row_indexes[iw1], iproj)
                                                * grad_overlap_2[dim]->get_value(col_indexes[iw2], jproj);
                                }
                            }
                        }
                    }

                    // HF term is minus, only one projector for each atom force.
                    f_delta(iat, 0) -= 2.0 * *dm_current * nlm[0];
                    f_delta(iat, 1) -= 2.0 * *dm_current * nlm[1];
                    f_delta(iat, 2) -= 2.0 * *dm_current * nlm[2];

                    // Pulay term is plus, only one projector for each atom force.
                    f_delta(ibt2, 0) += 2.0 * *dm_current * nlm[0];
                    f_delta(ibt2, 1) += 2.0 * *dm_current * nlm[1];
                    f_delta(ibt2, 2) += 2.0 * *dm_current * nlm[2];

                    if (isstress)
                    {
                        if (!PARAM.inp.deepks_equiv)
                        {
                            int ib = 0;
                            for (int L0 = 0; L0 <= orb.Alpha[0].getLmax(); ++L0)
                            {
                                for (int N0 = 0; N0 < orb.Alpha[0].getNchi(L0); ++N0)
                                {
                                    const int inl = inl_index[ucell.iat2it[iat]](ucell.iat2ia[iat], L0, N0);
                                    const int nm = 2 * L0 + 1;

                                    for (int m1 = 0; m1 < nm; ++m1)
                                    {
                                        for (int m2 = 0; m2 < nm; ++m2)
                                        {
                                            for (int dim = 0; dim < 3; ++dim)
                                            {
                                                nlm_t[dim]
                                                    += gedm[inl][m1 * nm + m2]
                                                       * overlap_2->get_value(col_indexes[iw2], ib + m1)
                                                       * grad_overlap_1[dim]->get_value(row_indexes[iw1], ib + m2);
                                            }
                                        }
                                    }
                                    ib += nm;
                                }
                            }
                            assert(ib == overlap_2->get_col_size());
                        }
                        else
                        {
                            int nproj = 0;
                            for (int il = 0; il < lmaxd + 1; il++)
                            {
                                nproj += (2 * il + 1) * orb.Alpha[0].getNchi(il);
                            }
                            for (int iproj = 0; iproj < nproj; iproj++)
                            {
                                for (int jproj = 0; jproj < nproj; jproj++)
                                {
                                    for (int dim = 0; dim < 3; dim++)
                                    {
                                        nlm_t[dim] += gedm[iat][iproj * nproj + jproj]
                                                      * overlap_2->get_value(col_indexes[iw2], iproj)
                                                      * grad_overlap_1[dim]->get_value(row_indexes[iw1], jproj);
                                    }
                                }
                            }
                        }

                        for (int ipol = 0; ipol < 3; ipol++)
                        {
                            for (int jpol = ipol; jpol < 3; jpol++)
                            {
                                svnl_dalpha(ipol, jpol)
                                    += *dm_current * (nlm[ipol] * r2[jpol] + nlm_t[ipol] * r1[jpol]);
                            }
                        }
                    }
                    dm_current++;
                } // iw2
            }     // iw1
        }
    );

    if (isstress)
    {
        assert(ucell.omega > 0.0);
        const double weight = ucell.lat0 / ucell.omega;
        // use upper triangle to make symmetric stress tensor
        for (int i = 0; i < 3; ++i)
        {
            for (int j = 0; j < 3; ++j)
            {
                if (j > i)
                {
                    svnl_dalpha(j, i) = svnl_dalpha(i, j);
                }
                svnl_dalpha(i, j) *= weight;
            }
        }
    }
    ModuleBase::timer::tick("DeePKS_domain", "cal_f_delta");
    return;
}

// prints forces and stress from DeePKS (LCAO)
void DeePKS_domain::check_f_delta(const int nat, ModuleBase::matrix& f_delta, ModuleBase::matrix& svnl_dalpha)
{
    ModuleBase::TITLE("DeePKS_domain", "check_F_delta");

    std::ofstream ofs("F_delta.dat");
    ofs << std::setprecision(10);

    for (int iat = 0; iat < nat; iat++)
    {
        ofs << f_delta(iat, 0) << " " << f_delta(iat, 1) << " " << f_delta(iat, 2) << std::endl;
    }

    std::ofstream ofs1("stress_delta.dat");
    ofs1 << std::setprecision(10);
    for (int ipol = 0; ipol < 3; ipol++)
    {
        for (int jpol = 0; jpol < 3; jpol++)
        {
            ofs1 << svnl_dalpha(ipol, jpol) << " ";
        }
        ofs1 << std::endl;
    }
    return;
}

template void DeePKS_domain::cal_f_delta<double>(const std::vector<std::vector<double>>& dm,
                                                 const UnitCell& ucell,
                                                 const LCAO_Orbitals& orb,
                                                 const Grid_Driver& GridD,
                                                 const Parallel_Orbitals& pv,
                                                 const int nks,
                                                 const std::vector<ModuleBase::Vector3<double>>& kvec_d,
                                                 std::vector<hamilt::HContainer<double>*> phialpha,
                                                 double** gedm,
                                                 ModuleBase::IntArray* inl_index,
                                                 ModuleBase::matrix& f_delta,
                                                 const bool isstress,
                                                 ModuleBase::matrix& svnl_dalpha);

template void DeePKS_domain::cal_f_delta<std::complex<double>>(const std::vector<std::vector<std::complex<double>>>& dm,
                                                               const UnitCell& ucell,
                                                               const LCAO_Orbitals& orb,
                                                               const Grid_Driver& GridD,
                                                               const Parallel_Orbitals& pv,
                                                               const int nks,
                                                               const std::vector<ModuleBase::Vector3<double>>& kvec_d,
                                                               std::vector<hamilt::HContainer<double>*> phialpha,
                                                               double** gedm,
                                                               ModuleBase::IntArray* inl_index,
                                                               ModuleBase::matrix& f_delta,
                                                               const bool isstress,
                                                               ModuleBase::matrix& svnl_dalpha);

#endif
