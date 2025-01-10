// wenfei 2022-1-11
// This file contains 2 subroutines:
// 1. build_phialpha, which calculates the overlap
// between atomic basis and projector alpha : <phi_mu|alpha>
// which will be used in calculating pdm, gdmx, H_V_delta, F_delta;
// 2. check_phialpha, which prints the results into .dat files
// for checking

#ifdef __DEEPKS

#include "deepks_phialpha.h"

#include "module_base/timer.h"
#include "module_base/vector3.h"
#include "module_parameter/parameter.h"

void DeePKS_domain::allocate_phialpha(const bool& cal_deri,
                                      const UnitCell& ucell,
                                      const LCAO_Orbitals& orb,
                                      const Grid_Driver& GridD,
                                      const Parallel_Orbitals* pv,
                                      std::vector<hamilt::HContainer<double>*>& phialpha)
{
    ModuleBase::TITLE("DeePKS_domain", "allocate_phialpha");

    phialpha.resize(cal_deri ? 4 : 1);

    phialpha[0] = new hamilt::HContainer<double>(pv); // phialpha is always real
    // Do not use fix_gamma, since it may find wrong matrix for gamma-only case in DeePKS

    // cutoff for alpha is same for all types of atoms
    const double Rcut_Alpha = orb.Alpha[0].getRcut();

    // Total number of atomic basis of projected orbitals
    int nw_alpha = 0;
    for (int l = 0; l <= orb.Alpha[0].getLmax(); l++)
    {
        nw_alpha += orb.Alpha[0].getNchi(l) * (2 * l + 1);
    }

    for (int T0 = 0; T0 < ucell.ntype; T0++)
    {
        Atom* atom0 = &ucell.atoms[T0];
        for (int I0 = 0; I0 < atom0->na; I0++)
        {
            const int iat = ucell.itia2iat(T0, I0);
            auto tau_a = atom0->tau[I0];
            GridD.Find_atom(ucell, tau_a, T0, I0);
            for (int ad = 0; ad < GridD.getAdjacentNum() + 1; ++ad)
            {
                const int T1 = GridD.getType(ad);
                const int I1 = GridD.getNatom(ad);
                const int ibt = ucell.itia2iat(T1, I1);
                auto tau_b = GridD.getAdjacentTau(ad);
                const int iw1_start = ucell.itiaiw2iwt(T1, I1, 0);
                const Atom* atom1 = &ucell.atoms[T1];
                auto R_index = GridD.getBox(ad);

                const double Rcut_AO1 = orb.Phi[T1].getRcut();
                const double dist = (tau_b - tau_a).norm() * ucell.lat0;
                if (dist > Rcut_Alpha + Rcut_AO1)
                {
                    continue;
                }

                hamilt::AtomPair<double> pair(iat, ibt, R_index, pv);
                // Notice: in AtomPair, the usage is set_size(ncol, nrow)
                pair.set_size(nw_alpha, atom1->nw * PARAM.globalv.npol);
                phialpha[0]->insert_pair(pair);
            }
        }
    }

    phialpha[0]->allocate(nullptr, true);
    // whether to calculate the derivative of phialpha
    if (cal_deri)
    {
        for (int i = 1; i < 4; ++i)
        {
            phialpha[i] = new hamilt::HContainer<double>(*phialpha[0], nullptr); // copy constructor
        }
    }
    return;
}

void DeePKS_domain::build_phialpha(const bool& cal_deri,
                                   const UnitCell& ucell,
                                   const LCAO_Orbitals& orb,
                                   const Grid_Driver& GridD,
                                   const Parallel_Orbitals* pv,
                                   const TwoCenterIntegrator& overlap_orb_alpha,
                                   std::vector<hamilt::HContainer<double>*>& phialpha)
{
    ModuleBase::TITLE("DeePKS_domain", "build_phialpha");
    ModuleBase::timer::tick("DeePKS_domain", "build_phialpha");

    // cutoff for alpha is same for all types of atoms
    const double Rcut_Alpha = orb.Alpha[0].getRcut();

    // Total number of atomic basis of projected orbitals
    // nw_alpha will be used frequently, better to add a function in Numerical Orbital to get it
    int nw_alpha = 0;
    for (int l = 0; l <= orb.Alpha[0].getLmax(); l++)
    {
        nw_alpha += orb.Alpha[0].getNchi(l) * (2 * l + 1);
    }

    for (int T0 = 0; T0 < ucell.ntype; T0++)
    {
        Atom* atom0 = &ucell.atoms[T0];
        for (int I0 = 0; I0 < atom0->na; I0++)
        {
            const int iat = ucell.itia2iat(T0, I0);
            auto tau_a = atom0->tau[I0];
            GridD.Find_atom(ucell, tau_a, T0, I0);
            for (int ad = 0; ad < GridD.getAdjacentNum() + 1; ++ad)
            {
                const int T1 = GridD.getType(ad);
                const int I1 = GridD.getNatom(ad);
                const int ibt = ucell.itia2iat(T1, I1);
                auto tau_b = GridD.getAdjacentTau(ad);
                const int iw1_start = ucell.itiaiw2iwt(T1, I1, 0);
                const Atom* atom1 = &ucell.atoms[T1];
                auto R_index = GridD.getBox(ad);
                int R[3] = {R_index.x, R_index.y, R_index.z};

                const double Rcut_AO1 = orb.Phi[T1].getRcut();
                const double dist = (tau_b - tau_a).norm() * ucell.lat0;
                if (dist > Rcut_Alpha + Rcut_AO1)
                {
                    continue;
                }

                double* data_pointer = phialpha[0]->data(iat, ibt, R);
                std::vector<double*> grad_pointer(3);
                if (cal_deri)
                {
                    for (int i = 0; i < 3; ++i)
                    {
                        grad_pointer[i] = phialpha[i + 1]->data(iat, ibt, R);
                    }
                }

                // Calculate phialpha
                // Get all indexes of atomic basis on the neighbour atom in this core
                // Notice that atom pair (a,b) and (b,a) are different when the former is changed to projected orbitals
                // So we need both row and col indexes for complete calculation
                auto all_indexes = pv->get_indexes_row(ibt);
                auto col_indexes = pv->get_indexes_col(ibt);
                all_indexes.insert(all_indexes.end(), col_indexes.begin(), col_indexes.end());
                std::sort(all_indexes.begin(), all_indexes.end());
                all_indexes.erase(std::unique(all_indexes.begin(), all_indexes.end()), all_indexes.end()); // for unique

                // inner loop : all atomic basis on the adjacent atom ad
                for (int iw1l = 0; iw1l < all_indexes.size(); iw1l += PARAM.globalv.npol)
                {
                    const int iw1 = all_indexes[iw1l] / PARAM.globalv.npol;

                    std::vector<std::vector<double>> nlm;
                    // 2D, dim 0 contains the overlap <phi_{iw1}|alpha_{all}>
                    // dim 1-3 contains the gradient of overlap

                    const int L1 = atom1->iw2l[iw1];
                    const int N1 = atom1->iw2n[iw1];
                    const int m1 = atom1->iw2m[iw1];

                    // convert m (0,1,...2l) to M (-l, -l+1, ..., l-1, l)
                    const int M1 = (m1 % 2 == 0) ? -m1 / 2 : (m1 + 1) / 2;

                    const ModuleBase::Vector3<double> dtau = tau_a - tau_b;
                    const int T0_fixed = 0; // All the projected orbitals are the same, use 0 here
                    overlap_orb_alpha.snap(T1, L1, N1, M1, T0_fixed, dtau * ucell.lat0, cal_deri, nlm);

                    const int index_begin = all_indexes[iw1l] * nw_alpha * PARAM.globalv.npol;
                    for (int iw0 = 0; iw0 < nw_alpha; iw0++)
                    {
                        data_pointer[index_begin + iw0] = nlm[0][iw0];
                        if (cal_deri)
                        {
                            grad_pointer[0][index_begin + iw0] = nlm[1][iw0];
                            grad_pointer[1][index_begin + iw0] = nlm[2][iw0];
                            grad_pointer[2][index_begin + iw0] = nlm[3][iw0];
                        }
                        if (PARAM.globalv.npol == 2)
                        {
                            data_pointer[index_begin + iw0 + nw_alpha] = nlm[0][iw0];
                            if (cal_deri)
                            {
                                grad_pointer[0][index_begin + iw0 + nw_alpha] = nlm[1][iw0];
                                grad_pointer[1][index_begin + iw0 + nw_alpha] = nlm[2][iw0];
                                grad_pointer[2][index_begin + iw0 + nw_alpha] = nlm[3][iw0];
                            }
                        }
                    }
                } // end iw
            }
        }
    }

    ModuleBase::timer::tick("DeePKS_domain", "build_phialpha");
    return;
}

void DeePKS_domain::check_phialpha(const bool& cal_deri,
                                   const UnitCell& ucell,
                                   const LCAO_Orbitals& orb,
                                   const Grid_Driver& GridD,
                                   const Parallel_Orbitals* pv,
                                   std::vector<hamilt::HContainer<double>*>& phialpha)
{
    ModuleBase::TITLE("DeePKS_domain", "check_phialpha");
    ModuleBase::timer::tick("DeePKS_domain", "check_phialpha");

    const double Rcut_Alpha = orb.Alpha[0].getRcut();
    // same for all types of atoms
    int nw_alpha = 0;
    for (int l = 0; l <= orb.Alpha[0].getLmax(); l++)
    {
        nw_alpha += orb.Alpha[0].getNchi(l) * (2 * l + 1);
    }

    std::ofstream ofs("phialpha.dat");
    std::ofstream ofs_x("dphialpha_x.dat");
    std::ofstream ofs_y("dphialpha_y.dat");
    std::ofstream ofs_z("dphialpha_z.dat");

    ofs << std::setprecision(10);
    ofs_x << std::setprecision(10);
    ofs_y << std::setprecision(10);
    ofs_z << std::setprecision(10);

    for (int T0 = 0; T0 < ucell.ntype; T0++)
    {
        Atom* atom0 = &ucell.atoms[T0];
        for (int I0 = 0; I0 < atom0->na; I0++)
        {
            const int iat = ucell.itia2iat(T0, I0);
            //=======================================================
            // Step 1 :
            // saves <alpha|phi>, where alpha runs over all projectors
            // and phi runs over atomic basis sets on the current core
            //=======================================================

            const ModuleBase::Vector3<double> tau0 = atom0->tau[I0];
            GridD.Find_atom(ucell, atom0->tau[I0], T0, I0);

            ofs << "iat : " << iat << std::endl;
            ofs_x << "iat : " << iat << std::endl;
            ofs_y << "iat : " << iat << std::endl;
            ofs_z << "iat : " << iat << std::endl;

            for (int ad = 0; ad < GridD.getAdjacentNum() + 1; ++ad)
            {
                const int T1 = GridD.getType(ad);
                const int I1 = GridD.getNatom(ad);
                const int start1 = ucell.itiaiw2iwt(T1, I1, 0);
                const double Rcut_AO1 = orb.Phi[T1].getRcut();

                const ModuleBase::Vector3<double> tau1 = GridD.getAdjacentTau(ad);
                const Atom* atom1 = &ucell.atoms[T1];
                const int nw1_tot = atom1->nw * PARAM.globalv.npol;

                const double dist1 = (tau1 - tau0).norm() * ucell.lat0;

                if (dist1 > Rcut_Alpha + Rcut_AO1)
                {
                    continue;
                }

                int ibt = ucell.itia2iat(T1, I1);
                int R[3];

                ofs << "ad : " << ad << " " << dist1 << std::endl;
                ofs_x << "ad : " << ad << " " << dist1 << std::endl;
                ofs_y << "ad : " << ad << " " << dist1 << std::endl;
                ofs_z << "ad : " << ad << " " << dist1 << std::endl;

                R[0] = GridD.getBox(ad).x;
                R[1] = GridD.getBox(ad).y;
                R[2] = GridD.getBox(ad).z;

                if (!PARAM.globalv.gamma_only_local)
                {
                    ofs << "R : " << R[0] << " " << R[1] << " " << R[2] << std::endl;
                    ofs_x << "R : " << R[0] << " " << R[1] << " " << R[2] << std::endl;
                    ofs_y << "R : " << R[0] << " " << R[1] << " " << R[2] << std::endl;
                    ofs_z << "R : " << R[0] << " " << R[1] << " " << R[2] << std::endl;
                }

                const double* data_pointer = phialpha[0]->data(iat, ibt, R);
                std::vector<double*> grad_pointer(3, nullptr);
                if (cal_deri)
                {
                    grad_pointer[0] = phialpha[1]->data(iat, ibt, R);
                    grad_pointer[1] = phialpha[2]->data(iat, ibt, R);
                    grad_pointer[2] = phialpha[3]->data(iat, ibt, R);
                }

                for (int iw1 = 0; iw1 < nw1_tot; ++iw1)
                {
                    const int iw1_all = start1 + iw1;
                    ofs << "iw : " << iw1_all << std::endl;
                    ofs_x << "iw : " << iw1_all << std::endl;
                    ofs_y << "iw : " << iw1_all << std::endl;
                    ofs_z << "iw : " << iw1_all << std::endl;
                    const int iw1_local = pv->global2local_row(iw1_all);
                    const int iw2_local = pv->global2local_col(iw1_all);
                    if (iw1_local < 0 && iw2_local < 0)
                        continue;

                    for (int ind = 0; ind < nw_alpha; ind++)
                    {
                        ofs << data_pointer[iw1 * nw_alpha + ind] << " ";
                        if (cal_deri)
                        {
                            ofs_x << grad_pointer[0][iw1 * nw_alpha + ind] << " ";
                            ofs_y << grad_pointer[1][iw1 * nw_alpha + ind] << " ";
                            ofs_z << grad_pointer[2][iw1 * nw_alpha + ind] << " ";
                        }
                        // 6 numbers per line
                        if (ind % 6 == 5)
                        {
                            ofs << "\n";
                            if (cal_deri)
                            {
                                ofs_x << "\n";
                                ofs_y << "\n";
                                ofs_z << "\n";
                            }
                        }
                    }
                    ofs << std::endl;
                    if (cal_deri)
                    {
                        ofs_x << std::endl;
                        ofs_y << std::endl;
                        ofs_z << std::endl;
                    }
                } // end iw
            }     // end ad
        }         // end I0
    }             // end T0

    ModuleBase::timer::tick("DeePKS_domain", "check_phialpha");
    return;
}

#endif
