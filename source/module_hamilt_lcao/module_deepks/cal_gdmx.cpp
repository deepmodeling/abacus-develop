#ifdef __DEEPKS

#include "LCAO_deepks.h"
#include "module_base/constants.h"
#include "module_base/libm/libm.h"
#include "module_base/parallel_reduce.h"
#include "module_base/timer.h"
#include "module_base/vector3.h"
#include "module_hamilt_lcao/module_hcontainer/atom_pair.h"
#include "module_parameter/parameter.h"

/// this subroutine calculates the gradient of projected density matrices
/// gdmx_m,m = d/dX sum_{mu,nu} rho_{mu,nu} <chi_mu|alpha_m><alpha_m'|chi_nu>

// There are 2 subroutines in this file:
// 1. cal_gdmx, calculating gdmx
// 2. check_gdmx, which prints gdmx to a series of .dat files

template <typename TK>
void LCAO_Deepks::cal_gdmx(const std::vector<std::vector<TK>>& dm,
                           const UnitCell& ucell,
                           const LCAO_Orbitals& orb,
                           const Grid_Driver& GridD,
                           const int nks,
                           const std::vector<ModuleBase::Vector3<double>>& kvec_d,
                           std::vector<hamilt::HContainer<double>*> phialpha,
                           torch::Tensor& gdmx)
{
    ModuleBase::TITLE("LCAO_Deepks", "cal_gdmx");
    ModuleBase::timer::tick("LCAO_Deepks", "cal_gdmx");
    // get DS_alpha_mu and S_nu_beta

    int nrow = this->pv->nrow;
    const int nm = 2 * lmaxd + 1;
    // gdmx: dD/dX
    // \sum_{mu,nu} 2*c_mu*c_nu * <dphi_mu/dx|alpha_m><alpha_m'|phi_nu>
    // size: [3][natom][tot_Inl][2l+1][2l+1]
    gdmx = torch::zeros({3, ucell.nat, inlmax, nm, nm}, torch::dtype(torch::kFloat64));
    auto accessor = gdmx.accessor<double, 5>();

    const double Rcut_Alpha = orb.Alpha[0].getRcut();

    for (int T0 = 0; T0 < ucell.ntype; T0++)
    {
        Atom* atom0 = &ucell.atoms[T0];
        for (int I0 = 0; I0 < atom0->na; I0++)
        {
            const int iat = ucell.itia2iat(T0, I0); // on which alpha is located
            const ModuleBase::Vector3<double> tau0 = atom0->tau[I0];
            GridD.Find_atom(ucell, atom0->tau[I0], T0, I0);

            for (int ad1 = 0; ad1 < GridD.getAdjacentNum() + 1; ++ad1)
            {
                const int T1 = GridD.getType(ad1);
                const int I1 = GridD.getNatom(ad1);
                const int ibt1 = ucell.itia2iat(T1, I1); // on which chi_mu is located
                const int start1 = ucell.itiaiw2iwt(T1, I1, 0);

                const ModuleBase::Vector3<double> tau1 = GridD.getAdjacentTau(ad1);
                const Atom* atom1 = &ucell.atoms[T1];
                const int nw1_tot = atom1->nw * PARAM.globalv.npol;
                const double Rcut_AO1 = orb.Phi[T1].getRcut();

                ModuleBase::Vector3<int> dR1(GridD.getBox(ad1).x, GridD.getBox(ad1).y, GridD.getBox(ad1).z);

                for (int ad2 = 0; ad2 < GridD.getAdjacentNum() + 1; ad2++)
                {
                    const int T2 = GridD.getType(ad2);
                    const int I2 = GridD.getNatom(ad2);
                    const int start2 = ucell.itiaiw2iwt(T2, I2, 0);
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

                    auto row_indexes = pv->get_indexes_row(ibt1);
                    auto col_indexes = pv->get_indexes_col(ibt2);
                    if (row_indexes.size() * col_indexes.size() == 0)
                    {
                        continue;
                    }

                    double* dm_current;
                    int dRx;
                    int dRy;
                    int dRz;
                    if constexpr (std::is_same<TK, double>::value)
                    {
                        dRx = 0;
                        dRy = 0;
                        dRz = 0;
                    }
                    else
                    {
                        dRx = (dR2 - dR1).x;
                        dRy = (dR2 - dR1).y;
                        dRz = (dR2 - dR1).z;
                    }
                    ModuleBase::Vector3<double> dR(dRx, dRy, dRz);

                    hamilt::AtomPair<double> dm_pair(ibt1, ibt2, dRx, dRy, dRz, pv);
                    dm_pair.allocate(nullptr, 1);
                    for (int ik = 0; ik < nks; ik++)
                    {
                        TK kphase;
                        if constexpr (std::is_same<TK, double>::value)
                        {
                            kphase = 1.0;
                        }
                        else
                        {
                            const double arg = -(kvec_d[ik] * dR) * ModuleBase::TWO_PI;
                            double sinp, cosp;
                            ModuleBase::libm::sincos(arg, &sinp, &cosp);
                            kphase = TK(cosp, sinp);
                        }
                        if (ModuleBase::GlobalFunc::IS_COLUMN_MAJOR_KS_SOLVER(PARAM.inp.ks_solver))
                        {
                            dm_pair.add_from_matrix(dm[ik].data(), pv->get_row_size(), kphase, 1);
                        }
                        else
                        {
                            dm_pair.add_from_matrix(dm[ik].data(), pv->get_col_size(), kphase, 0);
                        }
                    }

                    dm_current = dm_pair.get_pointer();

                    for (int iw1 = 0; iw1 < row_indexes.size(); ++iw1)
                    {
                        for (int iw2 = 0; iw2 < col_indexes.size(); ++iw2)
                        {
                            hamilt::BaseMatrix<double>* overlap_1 = phialpha[0]->find_matrix(iat, ibt1, dR1);
                            hamilt::BaseMatrix<double>* overlap_2 = phialpha[0]->find_matrix(iat, ibt2, dR2);
                            std::vector<hamilt::BaseMatrix<double>*> grad_overlap_1(3);
                            std::vector<hamilt::BaseMatrix<double>*> grad_overlap_2(3);

                            assert(overlap_1->get_col_size() == overlap_2->get_col_size());

                            for (int i = 0; i < 3; ++i)
                            {
                                grad_overlap_1[i] = phialpha[i + 1]->find_matrix(iat, ibt1, dR1);
                                grad_overlap_2[i] = phialpha[i + 1]->find_matrix(iat, ibt2, dR2);
                            }

                            int ib = 0;
                            for (int L0 = 0; L0 <= orb.Alpha[0].getLmax(); ++L0)
                            {
                                for (int N0 = 0; N0 < orb.Alpha[0].getNchi(L0); ++N0)
                                {
                                    const int inl = this->inl_index[T0](I0, L0, N0);
                                    const int nm = 2 * L0 + 1;
                                    for (int m1 = 0; m1 < nm; ++m1)
                                    {
                                        for (int m2 = 0; m2 < nm; ++m2)
                                        {
                                            for (int i = 0; i < 3; i++)
                                            {
                                                //(<d/dX chi_mu|alpha_m>)<chi_nu|alpha_m'>
                                                accessor[i][iat][inl][m1][m2]
                                                    += grad_overlap_2[i]->get_value(col_indexes[iw2], ib + m2)
                                                       * overlap_1->get_value(row_indexes[iw1], ib + m1) * *dm_current;

                                                //(<d/dX chi_nu|alpha_m'>)<chi_mu|alpha_m>
                                                accessor[i][iat][inl][m2][m1]
                                                    += grad_overlap_2[i]->get_value(col_indexes[iw2], ib + m2)
                                                       * overlap_1->get_value(row_indexes[iw1], ib + m1) * *dm_current;

                                                // (<chi_mu|d/dX alpha_m>)<chi_nu|alpha_m'> = -(<d/dX
                                                // chi_mu|alpha_m>)<chi_nu|alpha_m'>
                                                accessor[i][ibt2][inl][m1][m2]
                                                    -= grad_overlap_2[i]->get_value(col_indexes[iw2], ib + m2)
                                                       * overlap_1->get_value(row_indexes[iw1], ib + m1) * *dm_current;

                                                //(<chi_nu|d/dX alpha_m'>)<chi_mu|alpha_m> = -(<d/dX
                                                //chi_nu|alpha_m'>)<chi_mu|alpha_m>
                                                accessor[i][ibt2][inl][m2][m1]
                                                    -= grad_overlap_2[i]->get_value(col_indexes[iw2], ib + m2)
                                                       * overlap_1->get_value(row_indexes[iw1], ib + m1) * *dm_current;
                                            }
                                        }
                                    }
                                    ib += nm;
                                }
                            }
                            assert(ib == overlap_1->get_col_size());
                            dm_current++;
                        } // iw2
                    }     // iw1
                }         // ad2
            }             // ad1
        }                 // I0
    }                     // T0

#ifdef __MPI
    Parallel_Reduce::reduce_all(gdmx.data_ptr<double>(), 3 * ucell.nat * inlmax * nm * nm);
#endif
    ModuleBase::timer::tick("LCAO_Deepks", "cal_gdmx");
    return;
}

void LCAO_Deepks::check_gdmx(const int nat, const torch::Tensor& gdmx)
{
    std::stringstream ss;
    std::ofstream ofs_x;
    std::ofstream ofs_y;
    std::ofstream ofs_z;

    ofs_x << std::setprecision(10);
    ofs_y << std::setprecision(10);
    ofs_z << std::setprecision(10);

    const int nm = 2 * this->lmaxd + 1;
    auto accessor = gdmx.accessor<double, 5>();
    for (int ia = 0; ia < nat; ia++)
    {
        ss.str("");
        ss << "gdmx_" << ia << ".dat";
        ofs_x.open(ss.str().c_str());
        ss.str("");
        ss << "gdmy_" << ia << ".dat";
        ofs_y.open(ss.str().c_str());
        ss.str("");
        ss << "gdmz_" << ia << ".dat";
        ofs_z.open(ss.str().c_str());

        for (int inl = 0; inl < inlmax; inl++)
        {
            for (int m1 = 0; m1 < nm; m1++)
            {
                for (int m2 = 0; m2 < nm; m2++)
                {
                    ofs_x << accessor[0][ia][inl][m1][m2] << " ";
                    ofs_y << accessor[1][ia][inl][m1][m2] << " ";
                    ofs_z << accessor[2][ia][inl][m1][m2] << " ";
                }
            }
            ofs_x << std::endl;
            ofs_y << std::endl;
            ofs_z << std::endl;
        }
        ofs_x.close();
        ofs_y.close();
        ofs_z.close();
    }
}

template void LCAO_Deepks::cal_gdmx<double>(const std::vector<std::vector<double>>& dm,
                                            const UnitCell& ucell,
                                            const LCAO_Orbitals& orb,
                                            const Grid_Driver& GridD,
                                            const int nks,
                                            const std::vector<ModuleBase::Vector3<double>>& kvec_d,
                                            std::vector<hamilt::HContainer<double>*> phialpha,
                                            torch::Tensor& gdmx);

template void LCAO_Deepks::cal_gdmx<std::complex<double>>(const std::vector<std::vector<std::complex<double>>>& dm,
                                                          const UnitCell& ucell,
                                                          const LCAO_Orbitals& orb,
                                                          const Grid_Driver& GridD,
                                                          const int nks,
                                                          const std::vector<ModuleBase::Vector3<double>>& kvec_d,
                                                          std::vector<hamilt::HContainer<double>*> phialpha,
                                                          torch::Tensor& gdmx);

#endif
