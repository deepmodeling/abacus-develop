#include "module_parameter/parameter.h"

#ifdef __MLALGO

#include "deepks_force.h"
#include "deepks_iterate.h"
#include "source_base/constants.h"
#include "source_base/libm/libm.h"
#include "source_base/timer.h"
#include "source_base/vector3.h"
#include "source_lcao/module_hcontainer/atom_pair.h"

template <typename TK>
void DeePKS_domain::cal_f_delta(const hamilt::HContainer<double>* dmr,
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
    const double Rcut_Alpha = orb.Alpha[0].getRcut();

#pragma omp parallel
{
    ModuleBase::matrix f_delta_local(f_delta.nr, f_delta.nc);
    ModuleBase::matrix svnl_dalpha_local(svnl_dalpha.nr, svnl_dalpha.nc);
#pragma omp for schedule(dynamic)
    for (int iat = 0; iat < ucell.nat; iat++)
    {
        const int T0 = ucell.iat2it[iat];
        const int I0 = ucell.iat2ia[iat];
        Atom* atom0 = &ucell.atoms[T0];
        const ModuleBase::Vector3<double> tau0 = atom0->tau[I0];
        AdjacentAtomInfo adjs;
        GridD.Find_atom(ucell, atom0->tau[I0], T0, I0, &adjs);

        for (int ad1 = 0; ad1 < adjs.adj_num + 1; ++ad1)
        {
            const int T1 = adjs.ntype[ad1];
            const int I1 = adjs.natom[ad1];
            const int ibt1 = ucell.itia2iat(T1, I1);
            const ModuleBase::Vector3<double> tau1 = adjs.adjacent_tau[ad1];
            const Atom* atom1 = &ucell.atoms[T1];
            const int nw1_tot = atom1->nw * PARAM.globalv.npol;
            const double Rcut_AO1 = orb.Phi[T1].getRcut();

            ModuleBase::Vector3<int> dR1(adjs.box[ad1].x, adjs.box[ad1].y, adjs.box[ad1].z);

            for (int ad2 = 0; ad2 < adjs.adj_num + 1; ad2++)
            {
                const int T2 = adjs.ntype[ad2];
                const int I2 = adjs.natom[ad2];
                const int ibt2 = ucell.itia2iat(T2, I2);
                const ModuleBase::Vector3<double> tau2 = adjs.adjacent_tau[ad2];
                const Atom* atom2 = &ucell.atoms[T2];
                const int nw2_tot = atom2->nw * PARAM.globalv.npol;
                ModuleBase::Vector3<int> dR2(adjs.box[ad2].x, adjs.box[ad2].y, adjs.box[ad2].z);

                const double Rcut_AO2 = orb.Phi[T2].getRcut();
                const double dist1 = (tau1 - tau0).norm() * ucell.lat0;
                const double dist2 = (tau2 - tau0).norm() * ucell.lat0;

                if (dist1 > Rcut_Alpha + Rcut_AO1 || dist2 > Rcut_Alpha + Rcut_AO2)
                {
                    continue;
                }

                double r1[3]{};
                double r2[3]{};
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
                    continue;
                }

                int dRx = 0;
                int dRy = 0;
                int dRz = 0;
                if (std::is_same<TK, std::complex<double>>::value) // for multi-k
                {
                    dRx = dR1.x - dR2.x;
                    dRy = dR1.y - dR2.y;
                    dRz = dR1.z - dR2.z;
                }
                ModuleBase::Vector3<double> dR(dRx, dRy, dRz);
                const double* dm_current = dmr->find_matrix(ibt1, ibt2, dR.x, dR.y, dR.z)->get_pointer();

                hamilt::BaseMatrix<double>* overlap_1 = phialpha[0]->find_matrix(iat, ibt1, dR1);
                hamilt::BaseMatrix<double>* overlap_2 = phialpha[0]->find_matrix(iat, ibt2, dR2);
                if (overlap_1 == nullptr || overlap_2 == nullptr)
                {
                    continue;
                }
                std::vector<hamilt::BaseMatrix<double>*> grad_overlap_1(3);
                std::vector<hamilt::BaseMatrix<double>*> grad_overlap_2(3);
                for (int i = 0; i < 3; ++i)
                {
                    grad_overlap_1[i] = phialpha[i + 1]->find_matrix(iat, ibt1, dR1);
                    grad_overlap_2[i] = phialpha[i + 1]->find_matrix(iat, ibt2, dR2);
                }

                assert(overlap_1->get_col_size() == overlap_2->get_col_size());

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
                                                if (isstress)
                                                {
                                                    nlm_t[dim] += gedm[inl][m1 * nm + m2]
                                                                    * overlap_2->get_value(col_indexes[iw2], ib + m1)
                                                                    * grad_overlap_1[dim]->get_value(row_indexes[iw1],
                                                                                                    ib + m2);
                                                }
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
                                        if (isstress)
                                        {
                                            nlm_t[dim] += gedm[iat][iproj * nproj + jproj]
                                                            * overlap_2->get_value(col_indexes[iw2], iproj)
                                                            * grad_overlap_1[dim]->get_value(row_indexes[iw1], jproj);
                                        }
                                    }
                                }
                            }
                        }

                        // HF term is minus, only one projector for each atom force.
                        f_delta_local(iat, 0) -= 2.0 * *dm_current * nlm[0];
                        f_delta_local(iat, 1) -= 2.0 * *dm_current * nlm[1];
                        f_delta_local(iat, 2) -= 2.0 * *dm_current * nlm[2];

                        // Pulay term is plus, only one projector for each atom force.
                        f_delta_local(ibt2, 0) += 2.0 * *dm_current * nlm[0];
                        f_delta_local(ibt2, 1) += 2.0 * *dm_current * nlm[1];
                        f_delta_local(ibt2, 2) += 2.0 * *dm_current * nlm[2];

                        if (isstress)
                        {
                            for (int ipol = 0; ipol < 3; ipol++)
                            {
                                for (int jpol = ipol; jpol < 3; jpol++)
                                {
                                    svnl_dalpha_local(ipol, jpol)
                                        += *dm_current * (nlm[ipol] * r2[jpol] + nlm_t[ipol] * r1[jpol]);
                                }
                            }
                        }
                        dm_current++;
                    } // iw2
                }     // iw1
            }         // ad2
        }             // ad1
    }                 // iat
    if (isstress)
    {
        for (int ipol = 0; ipol < 3; ipol++)
        {
            for (int jpol = ipol; jpol < 3; jpol++)
            {
                #pragma omp atomic
                svnl_dalpha(ipol, jpol) += svnl_dalpha_local(ipol, jpol);
            }
        }
    }
    for (int iat = 0; iat < ucell.nat; iat++)
    {
        #pragma omp atomic
        f_delta(iat, 0) += f_delta_local(iat, 0);
        #pragma omp atomic
        f_delta(iat, 1) += f_delta_local(iat, 1);
        #pragma omp atomic
        f_delta(iat, 2) += f_delta_local(iat, 2);
    }
}

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

template void DeePKS_domain::cal_f_delta<double>(const hamilt::HContainer<double>* dmr,
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

template void DeePKS_domain::cal_f_delta<std::complex<double>>(const hamilt::HContainer<double>* dmr,
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
