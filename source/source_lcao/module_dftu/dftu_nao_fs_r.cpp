/// @file dftu_nao_fs_r.cpp
/// @brief DFT+U force and stress unified entry in real space (r-space) - implementation
///
/// See dftu_nao_fs_r.h for the mathematical formula and detailed documentation.

#include "dftu_nao_fs_r.h"
#include "dftu_nao_for_r.h"
#include "dftu_nao_str_r.h"
#include "dftu_nao_op.h"
#include "dftu_nao_pots.h"
#include "source_basis/module_nao/two_center_integrator.h"
#include "source_pw/module_pwdft/dftu_base.h"
#include "source_base/parallel_reduce.h"
#include "source_base/timer.h"
#include "source_cell/unitcell.h"
#include "source_hamilt/module_hcontainer/hcontainer.h"
#include "source_cell/module_neighbor/sltk_atom_arrange.h"

namespace DFTU_LCAO
{

void cal_fs_nao_r_impl(const UnitCell* ucell,
                       Plus_U_Base* dftu,
                       const TwoCenterIntegrator* intor,
                       int nspin,
                       const std::vector<AdjacentAtomInfo>& adjs_all,
                       const std::vector<const hamilt::HContainer<double>*>& dmR,
                       bool cal_force,
                       bool cal_stress,
                       ModuleBase::matrix& force,
                       ModuleBase::matrix& stress)
{
    ModuleBase::timer::start("DFTU", "cal_fs_nao_r");

    const Parallel_Orbitals* pv = dmR[0]->get_paraV();
    const int npol = ucell->get_npol();
    std::vector<double> stress_tmp;
    if (cal_stress)
    {
        stress_tmp.resize(6, 0.0);
    }
    if (cal_force)
    {
        force.zero_out();
    }
    // calculate atom_index for adjs_all, induced by omp parallel
    int atom_index = 0;
    std::vector<int> atom_index_all(ucell->nat, -1);
    for (int iat0 = 0; iat0 < ucell->nat; iat0++)
    {
        int T0 = 0;
        int I0 = 0;
        ucell->iat2iait(iat0, &I0, &T0);
        if (!dftu->has_l_channel(T0))
        {
            continue;
        }
        atom_index_all[iat0] = atom_index;
        atom_index++;
    }

    // 1. calculate <psi|beta> for each pair of atoms
    // loop over all on-site atoms
#pragma omp parallel
    {
        std::vector<double> stress_local;
        if (cal_stress)
        {
            stress_local.resize(6, 0.0);
        }
        ModuleBase::matrix force_local;
        if (cal_force)
        {
            force_local.create(force.nr, force.nc);
        }
#pragma omp for schedule(dynamic)
        for (int iat0 = 0; iat0 < ucell->nat; iat0++)
        {
            // skip the atoms without plus-U
            if (atom_index_all[iat0] < 0)
            {
                continue;
            }
            const ModuleBase::Vector3<double> tau0 = ucell->get_tau(iat0);
            int T0 = 0;
            int I0 = 0;
            ucell->iat2iait(iat0, &I0, &T0);
            const int target_L = dftu->get_l_channel(T0);
            const int tlp1 = 2 * target_L + 1;
            const AdjacentAtomInfo& adjs = adjs_all[atom_index_all[iat0]];

            std::vector<std::unordered_map<int, std::vector<double>>> nlm_tot;
            nlm_tot.resize(adjs.adj_num + 1);

            for (int ad = 0; ad < adjs.adj_num + 1; ++ad)
            {
                const int T1 = adjs.ntype[ad];
                const int I1 = adjs.natom[ad];
                const int iat1 = ucell->itia2iat(T1, I1);
                const ModuleBase::Vector3<double>& tau1 = adjs.adjacent_tau[ad];
                const Atom* atom1 = &ucell->atoms[T1];

                std::vector<int> all_indexes = pv->get_indexes_row(iat1);
                std::vector<int> col_indexes = pv->get_indexes_col(iat1);
                // insert col_indexes into all_indexes to get universal set with no repeat elements
                all_indexes.insert(all_indexes.end(), col_indexes.begin(), col_indexes.end());
                std::sort(all_indexes.begin(), all_indexes.end());
                all_indexes.erase(std::unique(all_indexes.begin(), all_indexes.end()), all_indexes.end());
                for (int iw1l = 0; iw1l < all_indexes.size(); iw1l += npol)
                {
                    const int iw1 = all_indexes[iw1l] / npol;
                    std::vector<std::vector<double>> nlm;
                    // nlm is a vector of vectors, but size of outer vector is only 1 here
                    // If we are calculating force, we need also to store the gradient
                    // and size of outer vector is then 4
                    // inner loop : all projectors (L0,M0)
                    int L1 = atom1->iw2l[iw1];
                    int N1 = atom1->iw2n[iw1];
                    int m1 = atom1->iw2m[iw1];

                    // convert m (0,1,...2l) to M (-l, -l+1, ..., l-1, l)
                    int M1 = (m1 % 2 == 0) ? -m1 / 2 : (m1 + 1) / 2;

                    ModuleBase::Vector3<double> dtau = tau0 - tau1;
                    intor->snap(T1, L1, N1, M1, T0, dtau * ucell->lat0,
                                1 /*cal_deri*/, nlm);

                    // select the elements of nlm with target_L
                    std::vector<double> nlm_target(tlp1 * 4);
                    const Atom* atom0 = &ucell->atoms[T0];
                    for (int iw = 0; iw < atom0->nw; iw++)
                    {
                        if (atom0->iw2l[iw] == target_L)
                        {
                            for (int n = 0; n < 4; n++) // value, deri_x, deri_y, deri_z
                            {
                                std::copy(nlm[n].begin() + iw, nlm[n].begin() + iw + tlp1,
                                          nlm_target.begin() + n * tlp1);
                            }
                            break;
                        }
                    }
                    nlm_tot[ad].insert({all_indexes[iw1l], nlm_target});
                }
            }
            // first iteration to calculate occupation matrix
            std::vector<double> occ(tlp1 * tlp1 * nspin, 0);
            dftu->occmat().get_flat(iat0, target_L, occ);

            // calculate pot_onsite
            const double u_value = dftu->get_u_current(T0);
            std::vector<double> pot_onsite(occ.size());
            double eu_tmp = 0;
            cal_pot_onsite(occ, tlp1, u_value, &pot_onsite[0], eu_tmp);

            // second iteration to calculate force and stress
            // calculate Force for atom J
            //     DMR_{I,J,R'-R} * <phi_{I,R}|chi_m> U*(1/2*delta(m, m')-occ(m, m'))
            //     d<chi_m'|phi_{J,R'}>/d tau_J for each pair of <IJR> atoms
            // calculate Stress for strain tensor epsilon_{alpha,beta}
            //     -1/Omega * DMR_{I,J,R'-R} * [ d<phi_{I,R}|chi_m>/d tau_{J,alpha} * tau_{J,beta}
            //     U*(1/2*delta(m, m')-occ(m, m'))<chi_m'|phi_{J,R'}>
            //   + <phi_{I,R}|chi_m> U*(1/2*delta(m, m')-occ(m, m'))
            //     d<chi_m'|phi_{J,R'}>/d tau_{J,alpha} * tau_{J,beta} ] for each pair of <IJR> atoms
            for (int ad1 = 0; ad1 < adjs.adj_num + 1; ++ad1)
            {
                const int T1 = adjs.ntype[ad1];
                const int I1 = adjs.natom[ad1];
                const int iat1 = ucell->itia2iat(T1, I1);
                double* force_tmp1 = (cal_force) ? &force_local(iat1, 0) : nullptr;
                double* force_tmp2 = (cal_force) ? &force_local(iat0, 0) : nullptr;
                const ModuleBase::Vector3<int>& R_index1 = adjs.box[ad1];
                ModuleBase::Vector3<double> dis1 = adjs.adjacent_tau[ad1] - tau0;
                for (int ad2 = 0; ad2 < adjs.adj_num + 1; ++ad2)
                {
                    const int T2 = adjs.ntype[ad2];
                    const int I2 = adjs.natom[ad2];
                    const int iat2 = ucell->itia2iat(T2, I2);
                    const ModuleBase::Vector3<int>& R_index2 = adjs.box[ad2];
                    ModuleBase::Vector3<double> dis2 = adjs.adjacent_tau[ad2] - tau0;
                    ModuleBase::Vector3<int> R_vector(R_index2[0] - R_index1[0],
                                                      R_index2[1] - R_index1[1],
                                                      R_index2[2] - R_index1[2]);
                    std::vector<const hamilt::BaseMatrix<double>*> tmp(nspin, nullptr);
                    tmp[0] = dmR[0]->find_matrix(iat1, iat2, R_vector[0], R_vector[1], R_vector[2]);
                    if (nspin == 2)
                    {
                        tmp[1] = dmR[1]->find_matrix(iat1, iat2, R_vector[0], R_vector[1], R_vector[2]);
                    }
                    // if not found , skip this pair of atoms
                    if (tmp[0] != nullptr)
                    {
                        // calculate force
                        if (cal_force)
                        {
                            cal_for_IJR_nao_r(iat1, iat2, pv,
                                            nlm_tot[ad1], nlm_tot[ad2],
                                            pot_onsite, tmp.data(), nspin,
                                            force_tmp1, force_tmp2);
                        }

                        // calculate stress
                        if (cal_stress)
                        {
                            cal_str_IJR_nao_r(iat1, iat2, pv,
                                             nlm_tot[ad1], nlm_tot[ad2],
                                             pot_onsite, tmp.data(), nspin,
                                             dis1, dis2, stress_local.data());
                        }
                    }
                }
            }
        }
#pragma omp critical
        {
            if (cal_force)
            {
                force += force_local;
            }
            if (cal_stress)
            {
                for (int i = 0; i < 6; i++)
                {
                    stress_tmp[i] += stress_local[i];
                }
            }
        }
    }

    if (cal_force)
    {
        Parallel_Reduce::reduce_all(force.c, force.nr * force.nc);
        if (nspin != 4)
        {
            for (int i = 0; i < force.nr * force.nc; i++)
            {
                force.c[i] *= 2.0;
            }
        }
    }

    // stress renormalization
    if (cal_stress)
    {
        Parallel_Reduce::reduce_all(stress_tmp.data(), 6);
        const double weight = ucell->lat0 / ucell->omega;
        for (int i = 0; i < 6; i++)
        {
            stress.c[i] = stress_tmp[i] * weight;
        }
        stress.c[8] = stress.c[5]; // stress(2,2)
        stress.c[7] = stress.c[4]; // stress(2,1)
        stress.c[6] = stress.c[2]; // stress(2,0)
        stress.c[5] = stress.c[4]; // stress(1,2)
        stress.c[4] = stress.c[3]; // stress(1,1)
        stress.c[3] = stress.c[1]; // stress(1,0)
    }

    ModuleBase::timer::end("DFTU", "cal_fs_nao_r");
}

void cal_fs_nao_r(const UnitCell* ucell,
                  Plus_U_Base* dftu,
                  const TwoCenterIntegrator* intor,
                  int nspin,
                  const std::vector<AdjacentAtomInfo>& adjs_all,
                  const std::vector<const hamilt::HContainer<double>*>& dmR,
                  bool cal_force,
                  bool cal_stress,
                  ModuleBase::matrix& force,
                  ModuleBase::matrix& stress)
{
    ModuleBase::TITLE("DFTU", "cal_fs_nao_r");

    if (dmR[0] == nullptr)
    {
        ModuleBase::WARNING_QUIT("DFTU", "dmr is not set");
    }
    if (dmR[0]->size_atom_pairs() == 0)
    {
        return;
    }

    cal_fs_nao_r_impl(ucell, dftu, intor, nspin, adjs_all, dmR,
                      cal_force, cal_stress, force, stress);
}

} // namespace DFTU_LCAO
