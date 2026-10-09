#include "source_base/timer.h"
#include "source_lcao/lcao_domain.h"
#include "source_io/module_parameter/parameter.h"

namespace LCAO_domain
{

typedef std::tuple<int, int, int, int> key_tuple;

#include "source_cell/record_adj.h" //mohan add 2012-07-06

// Per-neighbour context for Step 2 of build_Nonlocal_mu_new: everything that
// is fixed once the projector atom (T0, slot iat) and the two basis centres
// are known, before the (j, k) orbital loop runs. key1/key2 address the two
// <psi|beta> blocks inside nlm_tot / nlm_tot1.
struct NL_pair
{
    const Atom* atom1;
    const Atom* atom2;
    const int start1;
    const int start2;
    const int t0;
    const int iat;
    const key_tuple key1;
    const key_tuple key2;
};

// <psi1|beta><beta|psi2> contribution of one matrix element to the nonlocal
// energy (no derivatives). Defined statically here so the compiler can inline
// it at the hot inner-loop call site of build_Nonlocal_mu_new.
static void accum_nlm_energy(const NL_env& env,
                             const NL_elem& e,
                             const std::unordered_map<int, std::vector<double>>& nlm_cur1,
                             const std::unordered_map<int, std::vector<double>>& nlm_cur2,
                             double* NLloc)
{
    const std::unordered_map<int, std::vector<double>>::const_iterator it1 = nlm_cur1.find(e.iw1_all);
    const std::unordered_map<int, std::vector<double>>::const_iterator it2 = nlm_cur2.find(e.iw2_all);
    if (it1 == nlm_cur1.end() || it2 == nlm_cur2.end())
    {
        return;
    }
    const std::vector<double>& nlm_1 = it1->second;
    const std::vector<double>& nlm_2 = it2->second;
    double nlm_tmp = 0.0;
    const double* tmp_d = nullptr;
    for (int no = 0; no < env.ucell.atoms[e.t0].ncpp.non_zero_count_soc[0]; no++)
    {
        const int p1 = env.ucell.atoms[e.t0].ncpp.index1_soc[0][no];
        const int p2 = env.ucell.atoms[e.t0].ncpp.index2_soc[0][no];
        env.ucell.atoms[e.t0].ncpp.get_d(0, p1, p2, tmp_d);
        nlm_tmp += nlm_2[p2] * nlm_1[p1] * (*tmp_d);
    }

    if (nlm_tmp != 0.0)
    {
        if (env.gamma_only_local)
        {
            // mohan add 2010-12-20
            LCAO_domain::set_mat2d(e.iw1_all, e.iw2_all, nlm_tmp, env.pv, NLloc); // N stands for nonlocal.
        }
        else
        {
            NLloc[e.nnr] += nlm_tmp;
        }
    }
}

// force contribution of one matrix element for nspin==4 (SOC): the four
// spin blocks is0 fold the d_so matrices into DHloc_fixedR_{x,y,z}.
static void accum_nlm_force_soc(const NL_env& env,
                                const NL_elem& e,
                                const std::unordered_map<int, std::vector<std::vector<double>>>& nlm_cur1,
                                const std::unordered_map<int, std::vector<std::vector<double>>>& nlm_cur2,
                                const int is0,
                                ForceStressArrays& fsr)
{
    const std::unordered_map<int, std::vector<std::vector<double>>>::const_iterator it2 = nlm_cur2.find(e.iw2_all);
    const std::unordered_map<int, std::vector<std::vector<double>>>::const_iterator it1 = nlm_cur1.find(e.iw1_all);
    if (it2 == nlm_cur2.end() || it1 == nlm_cur1.end())
    {
        return;
    }
    const std::vector<double>& nlm_1 = it2->second[0];
    std::vector<std::vector<double>> nlm_2;
    nlm_2.resize(3);
    for (int i = 0; i < 3; i++)
    {
        nlm_2[i] = it1->second[i + 1];
    }
    for (int no = 0; no < env.ucell.atoms[e.t0].ncpp.non_zero_count_soc[is0]; no++)
    {
        const int p1 = env.ucell.atoms[e.t0].ncpp.index1_soc[is0][no];
        const int p2 = env.ucell.atoms[e.t0].ncpp.index2_soc[is0][no];
        double coef = 0.0;
        if (is0 == 0)
        {
            coef = (env.ucell.atoms[e.t0].ncpp.d_so(0, p2, p1).real()
                    + env.ucell.atoms[e.t0].ncpp.d_so(3, p2, p1).real())
                   * 0.5;
        }
        else if (is0 == 1)
        {
            coef = (env.ucell.atoms[e.t0].ncpp.d_so(1, p2, p1).real()
                    + env.ucell.atoms[e.t0].ncpp.d_so(2, p2, p1).real())
                   * 0.5;
        }
        else if (is0 == 2)
        {
            coef = (-env.ucell.atoms[e.t0].ncpp.d_so(1, p2, p1).imag()
                    + env.ucell.atoms[e.t0].ncpp.d_so(2, p2, p1).imag())
                   * 0.5;
        }
        else if (is0 == 3)
        {
            coef = (env.ucell.atoms[e.t0].ncpp.d_so(0, p2, p1).real()
                    - env.ucell.atoms[e.t0].ncpp.d_so(3, p2, p1).real())
                   * 0.5;
        }
        fsr.DHloc_fixedR_x[e.nnr] += nlm_2[0][p1] * nlm_1[p2] * coef;
        fsr.DHloc_fixedR_y[e.nnr] += nlm_2[1][p1] * nlm_1[p2] * coef;
        fsr.DHloc_fixedR_z[e.nnr] += nlm_2[2][p1] * nlm_1[p2] * coef;
    }
}

// force contribution of one matrix element for nspin==1/2. gamma_only writes
// the cartesian force arrays, the multi-k branch accumulates DHloc_fixedR.
static void accum_nlm_force(const NL_env& env,
                            const NL_elem& e,
                            const std::unordered_map<int, std::vector<std::vector<double>>>& nlm_cur1,
                            const std::unordered_map<int, std::vector<std::vector<double>>>& nlm_cur2,
                            ForceStressArrays& fsr)
{
    double nlm[3] = {0, 0, 0};
    const std::unordered_map<int, std::vector<std::vector<double>>>::const_iterator it1 = nlm_cur1.find(e.iw1_all);
    const std::unordered_map<int, std::vector<std::vector<double>>>::const_iterator it2 = nlm_cur2.find(e.iw2_all);
    if (it1 == nlm_cur1.end() || it2 == nlm_cur2.end())
    {
        return;
    }
    const std::vector<double>* nlm_1 = nullptr;
    std::vector<std::vector<double>> nlm_2;
    nlm_2.resize(3);
    if (env.gamma_only_local)
    {
        nlm_1 = &it1->second[0];
        for (int i = 0; i < 3; i++)
        {
            nlm_2[i] = it2->second[i + 1];
        }
    }
    else
    {
        // mohan change the order on 2011-06-17
        // origin: < psi1 | beta > < beta | dpsi2/dtau >
        // now: < psi1/dtau | beta > < beta | psi2 >
        nlm_1 = &it2->second[0];
        for (int i = 0; i < 3; i++)
        {
            nlm_2[i] = it1->second[i + 1];
        }
    }

    assert(nlm_1->size() == nlm_2[0].size());

    const double* tmp_d = nullptr;
    for (int no = 0; no < env.ucell.atoms[e.t0].ncpp.non_zero_count_soc[0]; no++)
    {
        const int p1 = env.ucell.atoms[e.t0].ncpp.index1_soc[0][no];
        const int p2 = env.ucell.atoms[e.t0].ncpp.index2_soc[0][no];
        env.ucell.atoms[e.t0].ncpp.get_d(0, p1, p2, tmp_d);
        for (int ir = 0; ir < 3; ir++)
        {
            nlm[ir] += nlm_2[ir][p2] * (*nlm_1)[p1] * (*tmp_d);
        }
    }

    if (env.gamma_only_local)
    {
        LCAO_domain::set_force(env.pv,
                               e.iw1_all,
                               e.iw2_all,
                               nlm[0],
                               nlm[1],
                               nlm[2],
                               'N',
                               fsr.DSloc_x.data(),
                               fsr.DSloc_y.data(),
                               fsr.DSloc_z.data(),
                               fsr.DHloc_fixed_x.data(),
                               fsr.DHloc_fixed_y.data(),
                               fsr.DHloc_fixed_z.data());
    }
    else
    {
        fsr.DHloc_fixedR_x[e.nnr] += nlm[0];
        fsr.DHloc_fixedR_y[e.nnr] += nlm[1];
        fsr.DHloc_fixedR_z[e.nnr] += nlm[2];
    }
}

// Step 1 of build_Nonlocal_mu_new: generate <psi|beta> (and, when derivatives
// are requested, <d psi|beta>) for every atom and its adjacent atoms. Fills
// nlm_tot (energy) or nlm_tot1 (force); only one of the two is populated,
// selected by env.calc_deri. The iat loop is OpenMP-parallel; each thread owns
// its own iat slot of the output containers.
static void build_psi_beta(
    const NL_env& env,
    const LCAO_Orbitals& orb,
    const TwoCenterIntegrator& intor_orb_beta,
    const Grid_Driver* GridD,
    std::vector<std::map<key_tuple, std::unordered_map<int, std::vector<double>>>>& nlm_tot,
    std::vector<std::map<key_tuple, std::unordered_map<int, std::vector<std::vector<double>>>>>& nlm_tot1)
{
    const UnitCell& ucell = env.ucell;
    const Parallel_Orbitals& pv = env.pv;
    const int npol = env.npol;
    const bool calc_deri = env.calc_deri;

#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic)
#endif
    for (int iat = 0; iat < ucell.nat; iat++)
    {
        const int it = ucell.iat2it[iat];
        const int ia = ucell.iat2ia[iat];

        const double Rcut_Beta = ucell.infoNL->get_rcut_max(it);
        const ModuleBase::Vector3<double> tau = ucell.atoms[it].tau[ia];
        AdjacentAtomInfo adjs;
        GridD->Find_atom(ucell, tau, it, ia, &adjs);

        if (!calc_deri)
        {
            nlm_tot[iat].clear();
        }
        else
        {
            nlm_tot1[iat].clear();
        }

        for (int ad = 0; ad < adjs.adj_num + 1; ++ad)
        {
            const int T1 = adjs.ntype[ad];
            const int I1 = adjs.natom[ad];
            const int start1 = ucell.itiaiw2iwt(T1, I1, 0);
            const double Rcut_AO1 = orb.Phi[T1].getRcut();

            const ModuleBase::Vector3<double>& tau1 = adjs.adjacent_tau[ad];
            const Atom* atom1 = &ucell.atoms[T1];
            const int nw1_tot = atom1->nw * npol;

            const ModuleBase::Vector3<double> dtau = tau1 - tau;
            const double dist1 = dtau.norm2() * pow(ucell.lat0, 2);

            if (dist1 > pow(Rcut_Beta + Rcut_AO1, 2))
            {
                continue;
            }
            std::unordered_map<int, std::vector<double>> nlm_cur;
            std::unordered_map<int, std::vector<std::vector<double>>> nlm_cur1;

            for (int iw1 = 0; iw1 < nw1_tot; ++iw1)
            {
                const int iw1_all = start1 + iw1;
                const int iw1_local = pv.global2local_row(iw1_all);
                const int iw2_local = pv.global2local_col(iw1_all);
                if (iw1_local < 0 && iw2_local < 0)
                {
                    continue;
                }
                const int iw1_0 = iw1 / npol;
                std::vector<std::vector<double>> nlm;
                // nlm is a vector of vectors, but size of outer vector is only 1 here
                // If we are calculating force, we need also to store the gradient
                // and size of outer vector is then 4
                // inner loop : all projectors (L0,M0)

                int L1 = atom1->iw2l[iw1_0];
                int N1 = atom1->iw2n[iw1_0];
                int m1 = atom1->iw2m[iw1_0];

                // convert m (0,1,...2l) to M (-l, -l+1, ..., l-1, l)
                int M1 = (m1 % 2 == 0) ? -m1 / 2 : (m1 + 1) / 2;

                ModuleBase::Vector3<double> dtau = tau - tau1;
                intor_orb_beta.snap(T1, L1, N1, M1, it, dtau * ucell.lat0, calc_deri, nlm);

                if (!calc_deri)
                {
                    nlm_cur.insert({iw1_all, nlm[0]});
                }
                else
                {
                    nlm_cur1.insert({iw1_all, nlm});
                }
            } // end iw

            const int iat1 = ucell.itia2iat(T1, I1);
            const int rx1 = adjs.box[ad].x;
            const int ry1 = adjs.box[ad].y;
            const int rz1 = adjs.box[ad].z;
            key_tuple key_1(iat1, rx1, ry1, rz1);

            if (!calc_deri)
            {
                nlm_tot[iat][key_1] = nlm_cur;
            }
            else
            {
                nlm_tot1[iat][key_1] = nlm_cur1;
            }
        } // end ad
    }
}

// Step 2 inner block: run the (j, k) orbital loop for one projector atom and
// accumulate every <psi1|beta><beta|psi2> matrix element into NLloc (energy)
// or fsr (force). Returns the number of elements written, i.e. the nnr_inner
// increment for this neighbour.
static int accum_nlm_block(
    const NL_env& env,
    const NL_pair& pair,
    const int nnr,
    std::vector<std::map<key_tuple, std::unordered_map<int, std::vector<double>>>>& nlm_tot,
    std::vector<std::map<key_tuple, std::unordered_map<int, std::vector<std::vector<double>>>>>& nlm_tot1,
    double* NLloc,
    ForceStressArrays& fsr)
{
    const Parallel_Orbitals& pv = env.pv;
    const int npol = env.npol;
    const int nspin = env.nspin;
    const bool calc_deri = env.calc_deri;

    int nnr_inner = 0;
    for (int j = 0; j < pair.atom1->nw * npol; j++)
    {
        const int j0 = j / npol; // added by zhengdy-soc
        const int iw1_all = pair.start1 + j;
        const int mu = pv.global2local_row(iw1_all);
        if (mu < 0)
        {
            continue;
        }

        // fix a serious bug: atom2[T2] -> atom2
        // mohan 2010-12-20
        for (int k = 0; k < pair.atom2->nw * npol; k++)
        {
            const int k0 = k / npol;
            const int iw2_all = pair.start2 + k;
            const int nu = pv.global2local_col(iw2_all);
            if (nu < 0)
            {
                continue;
            }

            const NL_elem elem{iw1_all, iw2_all, pair.t0, nnr + nnr_inner};
            if (!calc_deri)
            {
                accum_nlm_energy(env, elem, nlm_tot[pair.iat][pair.key1], nlm_tot[pair.iat][pair.key2], NLloc);
            }     // calc_deri
            else  // calculate the derivative
            {
                if (nspin == 4)
                {
                    const int is0 = (j - j0 * npol) + (k - k0 * npol) * 2;
                    accum_nlm_force_soc(env,
                                        elem,
                                        nlm_tot1[pair.iat][pair.key1],
                                        nlm_tot1[pair.iat][pair.key2],
                                        is0,
                                        fsr);
                }
                else if (nspin == 1 || nspin == 2)
                {
                    accum_nlm_force(env,
                                    elem,
                                    nlm_tot1[pair.iat][pair.key1],
                                    nlm_tot1[pair.iat][pair.key2],
                                    fsr);
                }
                else
                {
                    ModuleBase::WARNING_QUIT("LCAO_domain::build_Nonlocal_mu_new", "nspin must be 1, 2 or 4");
                }
            } //! calc_deri
            nnr_inner++;
        } // k
    }     // j
    return nnr_inner;
}

void build_Nonlocal_mu_new(const Parallel_Orbitals& pv,
                           ForceStressArrays& fsr,
                           double* NLloc,
                           const bool& calc_deri,
                           const UnitCell& ucell,
                           const LCAO_Orbitals& orb,
                           const TwoCenterIntegrator& intor_orb_beta,
                           const Grid_Driver* GridD)
{
    ModuleBase::TITLE("LCAO_domain", "vnl_mu_new");
    ModuleBase::timer::start("LCAO_domain", "vnl_mu_new");

    const int nspin = PARAM.inp.nspin;
    const int npol = PARAM.globalv.npol;
    const bool gamma_only_local = PARAM.globalv.gamma_only_local;

    const NL_env env{pv, ucell, nspin, npol, gamma_only_local, calc_deri};

    // < phi1 | beta > < beta | phi2 >
    // phi1 is within the unitcell.
    // while beta is in the supercell.
    // while phi2 is in the supercell.

    // Step 1 : generate <psi|beta>

    // This is the data structure for storing <psi|beta>
    // It is a 4 layer data structure
    // The outmost layer is std::vector with size being number of atoms in unit cell
    // The second layer is a map, the key being a combination of 4 number (iat, dRx, dRy, dRz)
    // which identifies a unique adjacent atom of the first atom
    // The third layer is an unordered map, with key being the index of atomic basis |psi>
    // The inner layer is a vector, each element representing a projector |beta>
    // It then either stores the number <psi|beta> (nlm_tot)
    // or a vector of 4, storing additionally <d/dx_i psi|beta> (nlm_tot1) x_i=x,y,z
    std::vector<std::map<key_tuple, std::unordered_map<int, std::vector<double>>>> nlm_tot;
    std::vector<std::map<key_tuple, std::unordered_map<int, std::vector<std::vector<double>>>>> nlm_tot1;

    if (!calc_deri)
    {
        nlm_tot.resize(ucell.nat);
    }
    else
    {
        nlm_tot1.resize(ucell.nat);
    }

    build_psi_beta(env, orb, intor_orb_beta, GridD, nlm_tot, nlm_tot1);

    //=======================================================
    // Step2:
    // calculate sum_(L0,M0) beta<psi_i|beta><beta|psi_j>
    // and accumulate the value to Hloc_fixedR(i,j)
    //=======================================================
    int total_nnr = 0;
#ifdef _OPENMP
#pragma omp parallel reduction(+ : total_nnr)
    {
#endif
        ModuleBase::Vector3<double> tau1, tau2, dtau;
        ModuleBase::Vector3<double> dtau1, dtau2, tau0;
        double distance = 0.0;
        double rcut = 0.0;
        double rcut1 = 0.0;
        double rcut2 = 0.0;

        //	Record_adj RA;
        //	RA.for_2d();

        // psi1
#ifdef _OPENMP
// use schedule(dynamic) for load balancing because adj_num is various
#pragma omp for schedule(dynamic)
#endif
        for (int iat1 = 0; iat1 < ucell.nat; iat1++)
        {
            const int T1 = ucell.iat2it[iat1];
            const Atom* atom1 = &ucell.atoms[T1];
            const int I1 = ucell.iat2ia[iat1];
            {
                // GridD->Find_atom( atom1->tau[I1] );
                AdjacentAtomInfo adjs;
                GridD->Find_atom(ucell, atom1->tau[I1], T1, I1, &adjs);
                const int start1 = ucell.itiaiw2iwt(T1, I1, 0);
                // Record_adj.for_2d() may not called in some case
                int nnr = 0;
                if (!pv.nlocstart.empty())
                {
                    nnr = pv.nlocstart[iat1];
                }
                tau1 = atom1->tau[I1];

                // psi2
                for (int ad2 = 0; ad2 < adjs.adj_num + 1; ++ad2)
                {
                    const int T2 = adjs.ntype[ad2];
                    const Atom* atom2 = &ucell.atoms[T2];

                    const int I2 = adjs.natom[ad2];
                    const int iat2 = ucell.itia2iat(T2, I2);
                    const int start2 = ucell.itiaiw2iwt(T2, I2, 0);
                    tau2 = adjs.adjacent_tau[ad2];

                    bool is_adj = false;

                    const int rx2 = adjs.box[ad2].x;
                    const int ry2 = adjs.box[ad2].y;
                    const int rz2 = adjs.box[ad2].z;

                    dtau = tau2 - tau1;
                    distance = dtau.norm2() * pow(ucell.lat0, 2);
                    // this rcut is in order to make nnr consistent
                    // with other matrix.
                    rcut = pow(orb.Phi[T1].getRcut() + orb.Phi[T2].getRcut(), 2);
                    if (distance < rcut)
                    {
                        is_adj = true;
                    }
                    else if (distance >= rcut)
                    {
                        for (int ad0 = 0; ad0 < adjs.adj_num + 1; ++ad0)
                        {
                            const int T0 = adjs.ntype[ad0];

                            tau0 = adjs.adjacent_tau[ad0];
                            dtau1 = tau0 - tau1;
                            dtau2 = tau0 - tau2;

                            const double distance1 = dtau1.norm2() * pow(ucell.lat0, 2);
                            const double distance2 = dtau2.norm2() * pow(ucell.lat0, 2);

                            rcut1 = pow(orb.Phi[T1].getRcut() + ucell.infoNL->get_rcut_max(T0), 2);
                            rcut2 = pow(orb.Phi[T2].getRcut() + ucell.infoNL->get_rcut_max(T0), 2);

                            if (distance1 < rcut1 && distance2 < rcut2)
                            {
                                is_adj = true;
                                break;
                            }
                        }
                    }

                    if (is_adj)
                    {
                        // < psi1 | all projectors | psi2 >
                        // ----------------------------- enter the nnr increaing zone -------------------------
                        for (int ad0 = 0; ad0 < adjs.adj_num + 1; ++ad0)
                        {
                            const int T0 = adjs.ntype[ad0];
                            const int I0 = adjs.natom[ad0];
                            const int iat = ucell.itia2iat(T0, I0);

                            // mohan add 2010-12-19
                            if (ucell.infoNL->get_nproj(T0) == 0)
                            {
                                continue;
                            }

                            tau0 = adjs.adjacent_tau[ad0];

                            dtau1 = tau0 - tau1;
                            dtau2 = tau0 - tau2;
                            const double distance1 = dtau1.norm2() * pow(ucell.lat0, 2);
                            const double distance2 = dtau2.norm2() * pow(ucell.lat0, 2);

                            // seems a bug here!! mohan 2011-06-17
                            rcut1 = pow(orb.Phi[T1].getRcut() + ucell.infoNL->get_rcut_max(T0), 2);
                            rcut2 = pow(orb.Phi[T2].getRcut() + ucell.infoNL->get_rcut_max(T0), 2);

                            if (distance1 >= rcut1 || distance2 >= rcut2)
                            {
                                continue;
                            }
                            // const Atom* atom0 = &ucell.atoms[T0];
                            const int rx0 = adjs.box[ad0].x;
                            const int ry0 = adjs.box[ad0].y;
                            const int rz0 = adjs.box[ad0].z;
                            key_tuple key1(iat1, -rx0, -ry0, -rz0);
                            key_tuple key2(iat2, rx2 - rx0, ry2 - ry0, rz2 - rz0);

                            const NL_pair pair{atom1, atom2, start1, start2, T0, iat, key1, key2};
                            accum_nlm_block(env, pair, nnr, nlm_tot, nlm_tot1, NLloc, fsr);
                        }         // ad0

                        // outer circle : accumulate nnr
                        for (int j = 0; j < atom1->nw * npol; j++)
                        {
                            const int j0 = j / npol; // added by zhengdy-soc
                            const int iw1_all = start1 + j;
                            const int mu = pv.global2local_row(iw1_all);
                            if (mu < 0)
                            {
                                continue;
                            }

                            // fix a serious bug: atom2[T2] -> atom2
                            // mohan 2010-12-20
                            for (int k = 0; k < atom2->nw * npol; k++)
                            {
                                const int k0 = k / npol;
                                const int iw2_all = start2 + k;
                                const int nu = pv.global2local_col(iw2_all);
                                if (nu < 0)
                                {
                                    continue;
                                }
                                total_nnr++;
                                nnr++;
                            }
                        }
                    } // end is_adj
                }     // ad2
            }         // I1
        }             // T1
#ifdef _OPENMP
    }
#endif
    if (!gamma_only_local)
    {
        if (total_nnr != pv.nnr)
        {
            GlobalV::ofs_running << " nr=" << total_nnr << std::endl;
            GlobalV::ofs_running << " pv->nnr=" << pv.nnr << std::endl;
            ModuleBase::WARNING_QUIT("LCAO_domain::build_Nonlocal_mu_new", "nnr!=LNNR.nnr");
        }
    }

    ModuleBase::timer::end("LCAO_domain", "vnl_mu_new");
    return;
}

} // namespace LCAO_domain
