#include <functional>

#ifdef __DEEPKS
#include "deepks_iterate.h"

void DeePKS_domain::iterate_ad1(const UnitCell& ucell,
                                const Grid_Driver& GridD,
                                const LCAO_Orbitals& orb,
                                const bool with_trace,
                                std::function<void(const int /*iat*/,
                                                   const ModuleBase::Vector3<double>& /*tau0*/,
                                                   const int /*ibt*/,
                                                   const ModuleBase::Vector3<double>& /*tau1*/,
                                                   const int /*start*/,
                                                   const int /*nw1_tot*/,
                                                   ModuleBase::Vector3<int> /*dR*/)> callback)
{
    const double Rcut_Alpha = orb.Alpha[0].getRcut();
    for (int iat = 0; iat < ucell.nat; iat++)
    {
        const int T0 = ucell.iat2it[iat];
        const int I0 = ucell.iat2ia[iat];
        Atom* atom0 = &ucell.atoms[T0];
        const ModuleBase::Vector3<double> tau0 = atom0->tau[I0];
        GridD.Find_atom(ucell, tau0, T0, I0);
        for (int ad = 0; ad < GridD.getAdjacentNum() + 1; ++ad)
        {
            const int T1 = GridD.getType(ad);
            const int I1 = GridD.getNatom(ad);
            const int ibt = ucell.itia2iat(T1, I1);
            const int start = ucell.itiaiw2iwt(T1, I1, 0);

            const ModuleBase::Vector3<double> tau1 = GridD.getAdjacentTau(ad);
            const Atom* atom1 = &ucell.atoms[T1];
            const int nw1_tot = atom1->nw * PARAM.globalv.npol;
            const double Rcut_AO1 = orb.Phi[T1].getRcut();
            const double dist1 = (tau1 - tau0).norm() * ucell.lat0;

            if (dist1 > Rcut_Alpha + Rcut_AO1)
            {
                continue;
            }

            ModuleBase::Vector3<int> dR(GridD.getBox(ad).x, GridD.getBox(ad).y, GridD.getBox(ad).z);

            callback(iat, tau0, ibt, tau1, start, nw1_tot, dR);
        }
    }
}

// void DeePKS_domain::iterate_ad1(const UnitCell& ucell,
//     const Grid_Driver& GridD,
//     const LCAO_Orbitals& orb,
//     const bool with_trace,
//     std::function<void(const int /*iat*/,
//                        const ModuleBase::Vector3<double>& /*tau0*/,
//                        const int /*ibt*/,
//                        const ModuleBase::Vector3<double>& /*tau1*/,
//                        const int /*start*/,
//                        const int /*nw1_tot*/,
//                        ModuleBase::Vector3<int> /*dR*/,
//                        std::vector<int> /*trace_alpha_row*/,
//                        std::vector<int> /*trace_alpha_col*/)> callback)
// {
//     const double Rcut_Alpha = orb.Alpha[0].getRcut();
//     for (int T0 = 0; T0 < ucell.ntype; T0++)
//     {
//         Atom* atom0 = &ucell.atoms[T0];
//         for (int I0 = 0; I0 < atom0->na; I0++)
//         {
//             const int iat = ucell.itia2iat(T0, I0);
//             const ModuleBase::Vector3<double> tau0 = atom0->tau[I0];
//             GridD.Find_atom(ucell, tau0, T0, I0);

//             // trace alpha orbital
//             std::vector<int> trace_alpha_row;
//             std::vector<int> trace_alpha_col;
//             int ib = 0;
//             for (int L0 = 0; L0 <= orb.Alpha[0].getLmax(); ++L0)
//             {
//                 for (int N0 = 0; N0 < orb.Alpha[0].getNchi(L0); ++N0)
//                 {
//                     const int inl = inl_index[T0](I0, L0, N0);
//                     const int nm = 2 * L0 + 1;

//                     for (int m1 = 0; m1 < nm; ++m1) // m1 = 1 for s, 3 for p, 5 for d
//                     {
//                         for (int m2 = 0; m2 < nm; ++m2) // m1 = 1 for s, 3 for p, 5 for d
//                         {
//                             trace_alpha_row.push_back(ib + m1);
//                             trace_alpha_col.push_back(ib + m2);
//                         }
//                     }
//                     ib += nm;
//                 }
//             }

//             for (int ad = 0; ad < GridD.getAdjacentNum() + 1; ++ad)
//             {
//                 const int T1 = GridD.getType(ad);
//                 const int I1 = GridD.getNatom(ad);
//                 const int ibt = ucell.itia2iat(T1, I1); // on which chi_mu is located
//                 const int start = ucell.itiaiw2iwt(T1, I1, 0);

//                 const ModuleBase::Vector3<double> tau1 = GridD.getAdjacentTau(ad);
//                 const Atom* atom1 = &ucell.atoms[T1];
//                 const int nw1_tot = atom1->nw * PARAM.globalv.npol;
//                 const double Rcut_AO1 = orb.Phi[T1].getRcut();
//                 const double dist1 = (tau1 - tau0).norm() * ucell.lat0;

//                 if (dist1 > Rcut_Alpha + Rcut_AO1)
//                 {
//                     continue;
//                 }

//                 ModuleBase::Vector3<int> dR(GridD.getBox(ad).x, GridD.getBox(ad).y, GridD.getBox(ad).z);

//                 callback(iat, tau0, ibt, tau1, start, nw1_tot, dR, trace_alpha_row, trace_alpha_col);
//             }
//         }
//     }
// }

void DeePKS_domain::iterate_ad2(const UnitCell& ucell,
                                const Grid_Driver& GridD,
                                const LCAO_Orbitals& orb,
                                const bool with_trace,
                                std::function<void(const int /*iat*/,
                                                   const ModuleBase::Vector3<double>& /*tau0*/,
                                                   const int /*ibt1*/,
                                                   const ModuleBase::Vector3<double>& /*tau1*/,
                                                   const int /*start1*/,
                                                   const int /*nw1_tot*/,
                                                   ModuleBase::Vector3<int> /*dR1*/,
                                                   const int /*ibt2*/,
                                                   const ModuleBase::Vector3<double>& /*tau2*/,
                                                   const int /*start2*/,
                                                   const int /*nw2_tot*/,
                                                   ModuleBase::Vector3<int> /*dR2*/)> callback)
{
    const double Rcut_Alpha = orb.Alpha[0].getRcut();
    DeePKS_domain::iterate_ad1(
        ucell,
        GridD,
        orb,
        with_trace,
        [&](const int iat,
            const ModuleBase::Vector3<double>& tau0,
            const int ibt1,
            const ModuleBase::Vector3<double>& tau1,
            const int start1,
            const int nw1_tot,
            ModuleBase::Vector3<int> dR1) {
            for (int ad = 0; ad < GridD.getAdjacentNum() + 1; ++ad)
            {
                const int T2 = GridD.getType(ad);
                const int I2 = GridD.getNatom(ad);
                const int start2 = ucell.itiaiw2iwt(T2, I2, 0);
                const int ibt2 = ucell.itia2iat(T2, I2);
                const ModuleBase::Vector3<double> tau2 = GridD.getAdjacentTau(ad);
                const Atom* atom2 = &ucell.atoms[T2];
                const int nw2_tot = atom2->nw * PARAM.globalv.npol;
                ModuleBase::Vector3<int> dR2(GridD.getBox(ad).x, GridD.getBox(ad).y, GridD.getBox(ad).z);

                const double Rcut_AO2 = orb.Phi[T2].getRcut();
                const double dist2 = (tau2 - tau0).norm() * ucell.lat0;

                if (dist2 >= Rcut_Alpha + Rcut_AO2)
                {
                    continue;
                }

                callback(iat, tau0, ibt1, tau1, start1, nw1_tot, dR1, ibt2, tau2, start2, nw2_tot, dR2);
            }
        }
    );
}

#endif