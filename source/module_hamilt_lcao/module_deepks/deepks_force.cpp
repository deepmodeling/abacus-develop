#include "module_parameter/parameter.h"

#ifdef __DEEPKS

#include "deepks_force.h"
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
                                const int lmaxd,
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

    f_delta.zero_out();

    const double Rcut_Alpha = orb.Alpha[0].getRcut();

    const int nrow = pv.nrow;

    for (int T0 = 0; T0 < ucell.ntype; T0++)
    {
        Atom* atom0 = &ucell.atoms[T0];
        for (int I0 = 0; I0 < atom0->na; I0++)
        {
            const int iat = ucell.itia2iat(T0, I0);
            const ModuleBase::Vector3<double> tau0 = atom0->tau[I0];
            GridD.Find_atom(ucell, atom0->tau[I0], T0, I0);

            for (int ad1 = 0; ad1 < GridD.getAdjacentNum() + 1; ++ad1)
            {
                const int T1 = GridD.getType(ad1);
                const int I1 = GridD.getNatom(ad1);
                const int ibt1 = ucell.itia2iat(T1, I1);
                const ModuleBase::Vector3<double> tau1 = GridD.getAdjacentTau(ad1);
                const Atom* atom1 = &ucell.atoms[T1];
                const int nw1_tot = atom1->nw * PARAM.globalv.npol;
                const double Rcut_AO1 = orb.Phi[T1].getRcut();

                ModuleBase::Vector3<int> dR1(GridD.getBox(ad1).x, GridD.getBox(ad1).y, GridD.getBox(ad1).z);

                for (int ad2 = 0; ad2 < GridD.getAdjacentNum() + 1; ad2++)
                {
                    const int T2 = GridD.getType(ad2);
                    const int I2 = GridD.getNatom(ad2);
                    const int ibt2 = ucell.itia2iat(T2, I2);
                    const ModuleBase::Vector3<double> tau2 = GridD.getAdjacentTau(ad2);
                    const Atom* atom2 = &ucell.atoms[T2];
                    const int nw2_tot = atom2->nw * PARAM.globalv.npol;
                    ModuleBase::Vector3<int> dR2(GridD.getBox(ad2).x, GridD.getBox(ad2).y, GridD.getBox(ad2).z);

                    const double Rcut_AO2 = orb.Phi[T2].getRcut();
                    const double dist1 = (tau1 - tau0).norm() * ucell.lat0;
                    const double dist2 = (tau2 - tau0).norm() * ucell.lat0;

                    if (dist1 > Rcut_Alpha + Rcut_AO1 || dist2 > Rcut_Alpha + Rcut_AO2)
                    {
                        continue;
                    }

                    double r0[3];
                    double r1[3];
                    if (isstress)
                    {
                        r1[0] = (tau1.x - tau0.x);
                        r1[1] = (tau1.y - tau0.y);
                        r1[2] = (tau1.z - tau0.z);
                        r0[0] = (tau2.x - tau0.x);
                        r0[1] = (tau2.y - tau0.y);
                        r0[2] = (tau2.z - tau0.z);
                    }

                    auto row_indexes = pv.get_indexes_row(ibt1);
                    auto col_indexes = pv.get_indexes_col(ibt2);

                    if (row_indexes.size() * col_indexes.size() == 0)
                    {
                        continue;
                    }

                    int dRx;
                    int dRy;
                    int dRz;
                    if constexpr (std::is_same<TK, double>::value) // for gamma-only
                    {
                        dRx = 0;
                        dRy = 0;
                        dRz = 0;
                    }
                    else // for multi-k
                    {
                        dRx = dR2.x - dR1.x;
                        dRy = dR2.y - dR1.y;
                        dRz = dR2.z - dR1.z;
                    }
                    ModuleBase::Vector3<double> dR(dRx, dRy, dRz);

                    hamilt::AtomPair<double> dm_pair(ibt1, ibt2, dRx, dRy, dRz, &pv);

                    dm_pair.allocate(nullptr, 1);

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

                    const double* dm_current = dm_pair.get_pointer();

                    for (int iw1 = 0; iw1 < row_indexes.size(); ++iw1)
                    {
                        for (int iw2 = 0; iw2 < col_indexes.size(); ++iw2)
                        {
                            double nlm[3] = {0, 0, 0};
                            double nlm_t[3] = {0, 0, 0}; // for stress

                            hamilt::BaseMatrix<double>* overlap_1 = phialpha[0]->find_matrix(iat, ibt1, dR1);
                            hamilt::BaseMatrix<double>* overlap_2 = phialpha[0]->find_matrix(iat, ibt2, dR2);
                            std::vector<hamilt::BaseMatrix<double>*> grad_overlap_1(3);
                            std::vector<hamilt::BaseMatrix<double>*> grad_overlap_2(3);
                            for (int i = 0; i < 3; ++i)
                            {
                                grad_overlap_1[i] = phialpha[i + 1]->find_matrix(iat, ibt1, dR1);
                                grad_overlap_2[i] = phialpha[i + 1]->find_matrix(iat, ibt2, dR2);
                            }

                            assert(overlap_1->get_col_size() == overlap_2->get_col_size());

                            if (!PARAM.inp.deepks_equiv)
                            {
                                int ib = 0;
                                for (int L0 = 0; L0 <= orb.Alpha[0].getLmax(); ++L0)
                                {
                                    for (int N0 = 0; N0 < orb.Alpha[0].getNchi(L0); ++N0)
                                    {
                                        const int inl = inl_index[T0](I0, L0, N0);
                                        const int nm = 2 * L0 + 1;
                                        for (int m1 = 0; m1 < nm; ++m1)
                                        {
                                            for (int m2 = 0; m2 < nm; ++m2)
                                            {
                                                for (int dim = 0; dim < 3; dim++)
                                                {
                                                    nlm[dim]
                                                        += gedm[inl][m1 * nm + m2]
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
                                            const int inl = inl_index[T0](I0, L0, N0);
                                            const int nm = 2 * L0 + 1;

                                            for (int m1 = 0; m1 < nm; ++m1)
                                            {
                                                for (int m2 = 0; m2 < nm; ++m2)
                                                {
                                                    for (int dim = 0; dim < 3; ++dim)
                                                    {
                                                        nlm_t[dim] += gedm[inl][m1 * nm + m2]
                                                                      * overlap_2->get_value(col_indexes[iw2], ib + m1)
                                                                      * grad_overlap_1[dim]->get_value(row_indexes[iw1],
                                                                                                       ib + m2);
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
                                            += *dm_current * (nlm[ipol] * r0[jpol] + nlm_t[ipol] * r1[jpol]);
                                    }
                                }
                            }
                            dm_current++;
                        } // iw2
                    }     // iw1
                }         // ad2
            }             // ad1
        }                 // end I0
    }                     // end T0

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

    return;
}

// prints forces and stress from DeePKS (LCAO)
void DeePKS_domain::check_f_delta(const int nat, ModuleBase::matrix& f_delta, ModuleBase::matrix& svnl_dalpha)
{
    ModuleBase::TITLE("LCAO_Deepks", "check_F_delta");

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
                                                 const int lmaxd,
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
                                                               const int lmaxd,
                                                               const int nks,
                                                               const std::vector<ModuleBase::Vector3<double>>& kvec_d,
                                                               std::vector<hamilt::HContainer<double>*> phialpha,
                                                               double** gedm,
                                                               ModuleBase::IntArray* inl_index,
                                                               ModuleBase::matrix& f_delta,
                                                               const bool isstress,
                                                               ModuleBase::matrix& svnl_dalpha);

#endif
