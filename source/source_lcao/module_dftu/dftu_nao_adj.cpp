#include "dftu_nao_adj.h"

#include "source_base/timer.h"
#include "source_base/tool_title.h"
#include "source_basis/module_ao/parallel_orbitals.h"
#include "source_basis/module_nao/two_center_integrator.h"
#include "source_cell/module_neighbor/sltk_grid_driver.h"
#include "source_cell/unitcell.h"
#include "source_pw/module_pwdft/dftu_base.h"

#include <algorithm>

namespace DFTU_LCAO
{

std::vector<AdjacentAtomInfo> build_adjacent_atoms(const UnitCell* ucell,
                                                   Plus_U_Base* dftu,
                                                   const Grid_Driver* gridD,
                                                   const std::vector<double>& orb_cutoff,
                                                   const double onsite_radius)
{
    std::vector<AdjacentAtomInfo> adjs_all;
    adjs_all.reserve(ucell->nat);
    for (int iat0 = 0; iat0 < ucell->nat; iat0++)
    {
        const ModuleBase::Vector3<double> tau0 = ucell->get_tau(iat0);
        int T0 = 0;
        int I0 = 0;
        ucell->iat2iait(iat0, &I0, &T0);
        if (!dftu->has_l_channel(T0))
        {
            continue;
        }

        AdjacentAtomInfo adjs;
        gridD->Find_atom(*ucell, tau0, T0, I0, &adjs);
        std::vector<bool> is_adj(adjs.adj_num + 1, false);
        for (int ad1 = 0; ad1 < adjs.adj_num + 1; ++ad1)
        {
            const int T1 = adjs.ntype[ad1];
            const int I1 = adjs.natom[ad1];
            const int iat1 = ucell->itia2iat(T1, I1);
            const ModuleBase::Vector3<int>& R_index1 = adjs.box[ad1];
            // choose the real adjacent atoms
            // Note: the distance of atoms should less than the cutoff radius,
            // When equal, the theoretical value of matrix element is zero,
            // but the calculated value is not zero due to the numerical error, which would lead to result changes.
            if (ucell->cal_dtau(iat0, iat1, R_index1).norm() * ucell->lat0
                < orb_cutoff[T1] + onsite_radius)
            {
                is_adj[ad1] = true;
            }
        }
        filter_adjs(is_adj, adjs);
        adjs_all.push_back(adjs);
    }
    return adjs_all;
}

NlmTot cal_nlm_all(const UnitCell& ucell,
                   const Plus_U_Base& dftu,
                   const TwoCenterIntegrator& intor,
                   const std::vector<AdjacentAtomInfo>& adjs_all,
                   const Parallel_Orbitals& pv)
{
    ModuleBase::TITLE("DFTU", "cal_nlm_all");
    ModuleBase::timer::start("DFTU", "cal_nlm_all");

    NlmTot nlm_tot(ucell.nat);
    const int npol = ucell.get_npol();
    int atom_index = 0;
    for (int iat0 = 0; iat0 < ucell.nat; iat0++)
    {
        const ModuleBase::Vector3<double> tau0 = ucell.get_tau(iat0);
        int T0 = 0;
        int I0 = 0;
        ucell.iat2iait(iat0, &I0, &T0);
        if (!dftu.has_l_channel(T0))
        {
            continue;
        }
        const int target_L = dftu.get_l_channel(T0);
        const int tlp1 = 2 * target_L + 1;
        const AdjacentAtomInfo& adjs = adjs_all[atom_index++];

        // calculate and save the table of two-center integrals
        nlm_tot[iat0].resize(adjs.adj_num + 1);

        for (int ad = 0; ad < adjs.adj_num + 1; ++ad)
        {
            const int T1 = adjs.ntype[ad];
            const int I1 = adjs.natom[ad];
            const int iat1 = ucell.itia2iat(T1, I1);
            const ModuleBase::Vector3<double>& tau1 = adjs.adjacent_tau[ad];
            const Atom* atom1 = &ucell.atoms[T1];

            std::vector<int> all_indexes = pv.get_indexes_row(iat1);
            std::vector<int> col_indexes = pv.get_indexes_col(iat1);
            // insert col_indexes into all_indexes to get universal set with no repeat elements
            all_indexes.insert(all_indexes.end(), col_indexes.begin(), col_indexes.end());
            std::sort(all_indexes.begin(), all_indexes.end());
            all_indexes.erase(std::unique(all_indexes.begin(), all_indexes.end()), all_indexes.end());
            for (int iw1l = 0; iw1l < all_indexes.size(); iw1l += npol)
            {
                const int iw1 = all_indexes[iw1l] / npol;
                // only first zeta orbitals in target L of atom iat0 are needed
                std::vector<double> nlm_target(tlp1);
                const int L1 = atom1->iw2l[iw1];
                const int N1 = atom1->iw2n[iw1];
                const int m1 = atom1->iw2m[iw1];
                std::vector<std::vector<double>> nlm;
                // nlm is a vector of vectors, but size of outer vector is only 1 here
                // If we are calculating force, we need also to store the gradient
                // and size of outer vector is then 4
                // inner loop : all projectors (L0,M0)

                // convert m (0,1,...2l) to M (-l, -l+1, ..., l-1, l)
                const int M1 = (m1 % 2 == 0) ? -m1 / 2 : (m1 + 1) / 2;

                ModuleBase::Vector3<double> dtau = tau0 - tau1;
                intor.snap(T1, L1, N1, M1, T0, dtau * ucell.lat0, false /*cal_deri*/, nlm);
                // select the elements of nlm with target_L
                for (int iw = 0; iw < ucell.atoms[T0].nw; iw++)
                {
                    const int L0 = ucell.atoms[T0].iw2l[iw];
                    if (L0 == target_L)
                    {
                        for (int m = 0; m < 2 * L0 + 1; m++)
                        {
                            nlm_target[m] = nlm[0][iw + m];
                        }
                        break;
                    }
                }
                nlm_tot[iat0][ad].insert({all_indexes[iw1l], nlm_target});
            }
        }
    }
    ModuleBase::timer::end("DFTU", "cal_nlm_all");
    return nlm_tot;
}

} // namespace DFTU_LCAO
