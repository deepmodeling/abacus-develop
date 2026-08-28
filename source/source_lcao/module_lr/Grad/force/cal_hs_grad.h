#pragma once
#include <cassert>
#include "source_cell/module_neighbor/sltk_grid_driver.h"
#include "source_cell/unitcell.h"
#include "source_basis/module_ao/parallel_orbitals.h"
#include "source_basis/module_nao/two_center_bundle.h"
#include "source_hamilt/module_hcontainer/hcontainer.h"

inline void filter_adjs_by_rcut(const UnitCell& ucell,
    const int iat0, 
    AdjacentAtomInfo& adjs)
{
    std::vector<bool> is_adj(adjs.adj_num + 1, false);
    for (int ad = 0;ad < adjs.adj_num + 1;++ad)
    {
        const int it0 = ucell.iat2it[iat0];
        const int it1 = adjs.ntype[ad];
        const int ia1 = adjs.natom[ad];
        const int iat1 = ucell.itia2iat(it1,ia1);
        const ModuleBase::Vector3<int>& R_index1 = adjs.box[ad];
        if(ucell.cal_dtau(iat0, iat1, R_index1).norm() * ucell.lat0 < ucell.atoms[it0].Rcut + ucell.atoms[it1].Rcut - 1e-15)
        {
            is_adj[ad] = true;
        }
    }
    filter_adjs(is_adj, adjs);
}

/// @brief Return a local Hamiltonian operator in type HContainer. 
/// It can replace OverlapNew::initialize_SR() and EkineticNew::initialize_HR().
template<typename TR>
hamilt::HContainer<TR> build_hcontainer_local_op(const UnitCell& ucell, const Grid_Driver& gd, const Parallel_Orbitals& pv)
{
    hamilt::HContainer<TR> hcontainer(&pv);
    for (int iat0 = 0;iat0 < ucell.nat; ++iat0)
    {
        const int it0 = ucell.iat2it[iat0];
        const int ia0 = ucell.iat2ia[iat0];
        AdjacentAtomInfo adjs;
        gd.Find_atom(ucell, ucell.get_tau(iat0), it0, ia0, &adjs);
        filter_adjs_by_rcut(ucell, iat0, adjs);

        for (int ad = 0;ad < adjs.adj_num + 1;++ad)
        {
            const int it1 = adjs.ntype[ad];
            const int ia1 = adjs.natom[ad];
            const int iat1 = ucell.itia2iat(it1, ia1);
            if (pv.get_nrow_atom(iat0) * pv.get_ncol_atom(iat1))
            {
                hamilt::AtomPair<TR> ap(iat0, iat1, adjs.box[ad], &pv);
                hcontainer.insert_pair(ap);
            }
        }
    }
    hcontainer.allocate(nullptr, true);
    return hcontainer;
}

/// @brief  Calculate <dPhi|phi> or <dPhi|T|phi> by 2-center integration
inline std::vector<hamilt::HContainer<double>> cal_hs_grad(const char job,
    const UnitCell& ucell,
    const Parallel_Orbitals& pv,
    const Grid_Driver& gd,
    const TwoCenterBundle& two_center_bundle)
{
    if(job != 'S' && job != 'T')
    {
        throw std::invalid_argument("job must be 'S' or 'T'");
    }
    // allocate dHS (better to use HContainer<Vector3<double>> for access continuity)
    // hamilt::HContainer<std::array<double, 3>> dHS(&pv);
    std::vector<hamilt::HContainer<double>> dHS(3, build_hcontainer_local_op<double>(ucell, gd, pv));
    // std::vector<hamilt::HContainer<double>> dHS(3, build_hcontainer_local_op<double>(ucell, gd, pv));
    std::vector<double> tmp_deriv(3);
    const int npol = ucell.get_npol();

    for (int ixyz = 0;ixyz < 3;++ixyz)
    {
        // traverse ijR to calculate dHS
        std::vector<int> ijr_info = dHS.at(ixyz).get_ijr_info();
        // std::vector<int> ijr_info = dHS.at(0).get_ijr_info();
        std::vector<int>::iterator it = ijr_info.begin();
        const int npairs = *it++;
        int npairs_count = 0;
        while (it != ijr_info.end())
        {
            ++npairs_count;
            const int iat0 = *it++;
            const int it0 = ucell.iat2it[iat0];
            const Atom& atom0 = ucell.atoms[it0];
            const ModuleBase::Vector3<double> tau0 = ucell.get_tau(iat0);
            auto row_indexes = pv.get_indexes_row(iat0);

            const int iat1 = *it++;
            const int it1 = ucell.iat2it[iat1];
            const Atom& atom1 = ucell.atoms[it1];
            const ModuleBase::Vector3<double> tau1 = ucell.get_tau(iat1);
            auto col_indexes = pv.get_indexes_col(iat1);

            const int nR = *it++;

            for (int iR = 0;iR < nR;++iR)
            {
                // Read the three components in separate statements. The order in which function
                // arguments are evaluated is UNSPECIFIED in C++, so `Vector3<double>(*it++, *it++, *it++)`
                // may store the R triple permuted and silently corrupt periodic systems.
                const int Rx = *it++;
                const int Ry = *it++;
                const int Rz = *it++;
                const ModuleBase::Vector3<double> R(Rx, Ry, Rz);    // int to double
                ModuleBase::Vector3<double> relative_position = (tau1 - tau0 + R * ucell.latvec) * ucell.lat0;
                hamilt::BaseMatrix<double>* dHS_block = dHS[ixyz].find_matrix(iat0, iat1, R.x, R.y, R.z);
                // `ijr_info` came from this very container, so a miss means the indices are wrong.
                assert(dHS_block != nullptr);

                // OMP can be used here
                for (int lw0 = 0;lw0 < row_indexes.size();lw0 += npol)    // spin 1-3 of dHS is not needed at nspin=4
                {
                    const int gw0 = row_indexes[lw0] / npol;
                    const int l0 = atom0.iw2l[gw0];
                    const int n0 = atom0.iw2n[gw0];
                    const int m0 = atom0.iw2m[gw0];
                    const int M0 = (m0 % 2 == 0) ? -m0 / 2 : (m0 + 1) / 2;      // convert m (0,1,...2l) to M (-l, -l+1, ..., l-1, l)
                    for (int lw1 = 0;lw1 < col_indexes.size();lw1 += npol)
                    {
                        const int& gw1 = col_indexes[lw1] / npol;
                        const int l1 = atom1.iw2l[gw1];
                        const int n1 = atom1.iw2n[gw1];
                        const int m1 = atom1.iw2m[gw1];
                        const int M1 = (m1 % 2 == 0) ? -m1 / 2 : (m1 + 1) / 2;
                        switch (job)
                        {
                        case 'S':
                            two_center_bundle.overlap_orb->calculate(it0, l0, n0, M0, it1, l1, n1, M1,
                                relative_position, nullptr, tmp_deriv.data());
                            break;
                        case 'T':
                            two_center_bundle.kinetic_orb->calculate(it0, l0, n0, M0, it1, l1, n1, M1,
                                relative_position, nullptr, tmp_deriv.data());
                            break;
                        }
                        dHS_block->get_value(lw0, lw1) = tmp_deriv[ixyz];
                    }
                }
            }
        }
        assert(npairs == npairs_count);
    }
    return dHS;
}