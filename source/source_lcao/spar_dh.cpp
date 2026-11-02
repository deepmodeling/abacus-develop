#include "spar_dh.h"

#include "source_base/parallel_reduce.h"
#include "source_base/timer.h"
#include "source_io/module_parameter/parameter.h"
#include "source_lcao/lcao_domain.h"
#include "source_hamilt/module_gint/gint_dvlocal.h"

#include <cmath>
#include <map>
#include <vector>

namespace
{

using SparseMatrixR = std::map<Abfs::Vector3_Order<int>, std::map<size_t, std::map<size_t, double>>>;

// Adds one Cartesian component of the grid-integrated <phi|V|dphi> to dH{x,y,z}R_sparse,
// completing the lower triangle from -R and reducing over the pool.
void distribute_pvdpR_sparse(const int cspin,
                             const int dim,
                             const int nlocal,
                             const double sparse_thr,
                             const SparseMatrixR& pvdpR_sparse,
                             const Parallel_Orbitals& pv,
                             LCAO_HS_Arrays& hs_arrays)
{
    int total_R_num = hs_arrays.all_R_coor.size();
    std::vector<int> nonzero_num(total_R_num);
    std::vector<int> minus_nonzero_num(total_R_num);
    int count = 0;
    for (const auto& R_coor: hs_arrays.all_R_coor)
    {
        auto iter = pvdpR_sparse.find(R_coor);
        if (iter != pvdpR_sparse.end())
        {
            for (auto& row_loop: iter->second)
            {
                nonzero_num[count] += row_loop.second.size();
            }
        }

        auto minus_R_coor = -1 * R_coor;

        iter = pvdpR_sparse.find(minus_R_coor);
        if (iter != pvdpR_sparse.end())
        {
            for (auto& row_loop: iter->second)
            {
                minus_nonzero_num[count] += row_loop.second.size();
            }
        }
        count++;
    }

    Parallel_Reduce::reduce_all(nonzero_num.data(), total_R_num);
    Parallel_Reduce::reduce_all(minus_nonzero_num.data(), total_R_num);

    std::vector<double> tmp(nlocal);
    count = 0;

    const std::vector<int>& trace_lo = ModuleGint::Gint::get_gint_info().get_trace_lo();
    for (const auto& R_coor: hs_arrays.all_R_coor)
    {
        if (nonzero_num[count] != 0 || minus_nonzero_num[count] != 0)
        {
            auto minus_R_coor = -1 * R_coor;

            for (int row = 0; row < nlocal; ++row)
            {
                tmp.assign(tmp.size(), 0);

                auto iter = pvdpR_sparse.find(R_coor);
                if (iter != pvdpR_sparse.end())
                {

                    if (trace_lo[row] >= 0)
                    {
                        auto row_iter = iter->second.find(row);
                        if (row_iter != iter->second.end())
                        {
                            for (auto& value: row_iter->second)
                            {
                                tmp[value.first] = value.second;
                            }
                        }
                    }
                }

                auto minus_R_iter = pvdpR_sparse.find(minus_R_coor);
                if (minus_R_iter != pvdpR_sparse.end())
                {
                    for (int col = 0; col < row; ++col)
                    {
                        if (trace_lo[col] >= 0)
                        {
                            auto row_iter = minus_R_iter->second.find(col);
                            if (row_iter != minus_R_iter->second.end())
                            {
                                auto col_iter = row_iter->second.find(row);
                                if (col_iter != row_iter->second.end())
                                {
                                    tmp[col] = col_iter->second;
                                }
                            }
                        }
                    }
                }

                Parallel_Reduce::reduce_pool(tmp.data(), nlocal);

                if (pv.global2local_row(row) >= 0)
                {
                    for (int col = 0; col < nlocal; ++col)
                    {
                        if (pv.global2local_col(col) >= 0)
                        {
                            if (std::abs(tmp[col]) > sparse_thr)
                            {
                                if (dim == 0)
                                {
                                    double& value = hs_arrays.dHRx_sparse[cspin][R_coor][row][col];
                                    value += tmp[col];
                                    if (std::abs(value) <= sparse_thr)
                                    {
                                        hs_arrays.dHRx_sparse[cspin][R_coor][row].erase(col);
                                    }
                                }
                                if (dim == 1)
                                {
                                    double& value = hs_arrays.dHRy_sparse[cspin][R_coor][row][col];
                                    value += tmp[col];
                                    if (std::abs(value) <= sparse_thr)
                                    {
                                        hs_arrays.dHRy_sparse[cspin][R_coor][row].erase(col);
                                    }
                                }
                                if (dim == 2)
                                {
                                    double& value = hs_arrays.dHRz_sparse[cspin][R_coor][row][col];
                                    value += tmp[col];
                                    if (std::abs(value) <= sparse_thr)
                                    {
                                        hs_arrays.dHRz_sparse[cspin][R_coor][row].erase(col);
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
        count++;
    }
}

// Converts the grid-integrated local-potential derivative <phi|V|dphi>, held per atom pair in
// pvdpR{x,y,z}, to the global-index sparse format of LCAO_HS_Arrays.
void cal_dvlocal_R_sparse(const hamilt::HContainer<double>& pvdpRx,
                          const hamilt::HContainer<double>& pvdpRy,
                          const hamilt::HContainer<double>& pvdpRz,
                          const int npol,
                          const int cspin,
                          const int nlocal,
                          const double sparse_thr,
                          const Parallel_Orbitals& pv,
                          const UnitCell& ucell,
                          LCAO_HS_Arrays& hs_arrays)
{
    ModuleBase::TITLE("Gint", "cal_dvlocal_R_sparse");
    ModuleBase::timer::start("Gint", "cal_dvlocal_R_sparse");
    SparseMatrixR pvdpRx_sparse;
    SparseMatrixR pvdpRy_sparse;
    SparseMatrixR pvdpRz_sparse;

    for (int iap = 0; iap < pvdpRx.size_atom_pairs(); iap++)
    {
        const auto& ap = pvdpRx.get_atom_pair(iap);
        const int iat1 = ap.get_atom_i();
        const int iat2 = ap.get_atom_j();
        const int it1 = ucell.iat2it[iat1];
        const int it2 = ucell.iat2it[iat2];
        const Atom* atom1 = &ucell.atoms[it1];
        const Atom* atom2 = &ucell.atoms[it2];
        const int start1 = ucell.itiaiw2iwt(it1, ucell.iat2ia[iat1], 0);
        const int start2 = ucell.itiaiw2iwt(it2, ucell.iat2ia[iat2], 0);

        for (int ir = 0; ir < ap.get_R_size(); ir++)
        {
            const ModuleBase::Vector3<int> R = ap.get_R_index(ir);
            Abfs::Vector3_Order<int> dR(R.x, R.y, R.z);
            double* p_pvdpRx = pvdpRx.get_atom_pair(iap).get_pointer(ir);
            double* p_pvdpRy = pvdpRy.get_atom_pair(iap).get_pointer(ir);
            double* p_pvdpRz = pvdpRz.get_atom_pair(iap).get_pointer(ir);

            for (int iw = 0; iw < atom1->nw * npol; iw++)
            {
                for (int iw2 = 0; iw2 < atom2->nw * npol; iw2++)
                {
                    const int nw = atom2->nw;
                    const int mug0 = iw / npol;
                    const int nug0 = iw2 / npol;
                    const int iw_nowg = mug0 * nw + nug0;

                    double temp_value = p_pvdpRx[iw_nowg];
                    if (std::abs(temp_value) > sparse_thr)
                    {
                        pvdpRx_sparse[dR][start1 + iw][start2 + iw2] = temp_value;
                    }
                    temp_value = p_pvdpRy[iw_nowg];
                    if (std::abs(temp_value) > sparse_thr)
                    {
                        pvdpRy_sparse[dR][start1 + iw][start2 + iw2] = temp_value;
                    }
                    temp_value = p_pvdpRz[iw_nowg];
                    if (std::abs(temp_value) > sparse_thr)
                    {
                        pvdpRz_sparse[dR][start1 + iw][start2 + iw2] = temp_value;
                    }
                }
            }
        }
    }
    distribute_pvdpR_sparse(cspin, 0, nlocal, sparse_thr, pvdpRx_sparse, pv, hs_arrays);
    distribute_pvdpR_sparse(cspin, 1, nlocal, sparse_thr, pvdpRy_sparse, pv, hs_arrays);
    distribute_pvdpR_sparse(cspin, 2, nlocal, sparse_thr, pvdpRz_sparse, pv, hs_arrays);
    ModuleBase::timer::end("Gint", "cal_dvlocal_R_sparse");
}

} // namespace

void sparse_format::cal_dS(const UnitCell& ucell,
    const Parallel_Orbitals& pv,
    LCAO_HS_Arrays& HS_Arrays,
    const Grid_Driver& grid,
    const TwoCenterBundle& two_center_bundle,
    const LCAO_Orbitals& orb,
    const double& sparse_thr,
    const bool gamma_only_local,
    const int nspin,
    const int npol)
{
ModuleBase::TITLE("sparse_format", "cal_dS");

sparse_format::set_R_range(HS_Arrays.all_R_coor, grid);
const int nnr = gamma_only_local ? pv.nloc : pv.nnr;

ForceStressArrays fsr_dh;
fsr_dh.DHloc_fixedR_x.resize(nnr, 0.0);
fsr_dh.DHloc_fixedR_y.resize(nnr, 0.0);
fsr_dh.DHloc_fixedR_z.resize(nnr, 0.0);
// cal dS=<phi|dphi> in LCAO
const bool cal_deri = true;
const bool cal_stress = false;
LCAO_domain::build_ST_new(fsr_dh,
       'S',
       cal_deri,
       cal_stress,
       ucell,
       orb,
       pv,
       two_center_bundle,
       &grid,
       nullptr,
       false); // delete unused parameter lm.Hloc_fixedR

sparse_format::cal_dSTN_R(ucell, pv, HS_Arrays, fsr_dh, grid, orb.cutoffs(), 0, sparse_thr, nspin, npol);
return;
}
void sparse_format::cal_dH(const UnitCell& ucell,
                           const Parallel_Orbitals& pv,
                           LCAO_HS_Arrays& HS_Arrays,
                           const Grid_Driver& grid,
                           const TwoCenterBundle& two_center_bundle,
                           const LCAO_Orbitals& orb,
                           const int& current_spin,
                           const double& sparse_thr,
                           const ModuleBase::matrix& v_eff,
                           const bool gamma_only_local,
                           const int nspin,
                           const int npol)
{
    ModuleBase::TITLE("sparse_format", "cal_dH");

    sparse_format::set_R_range(HS_Arrays.all_R_coor, grid);

    const int nnr = gamma_only_local ? pv.nloc : pv.nnr;

    ForceStressArrays fsr_dh;

    fsr_dh.DHloc_fixedR_x.resize(nnr, 0.0);
    fsr_dh.DHloc_fixedR_y.resize(nnr, 0.0);
    fsr_dh.DHloc_fixedR_z.resize(nnr, 0.0);
    // cal dT=<phi|kin|dphi> in LCAO
    // cal T + VNL(P1) in LCAO basis
    const bool cal_deri = true;
    const bool cal_stress = false;
    LCAO_domain::build_ST_new(fsr_dh,
                              'T',
                              cal_deri,
                              cal_stress,
                              ucell,
                              orb,
                              pv,
                              two_center_bundle,
                              &grid,
                              nullptr,
                              false); // delete unused parameter lm.Hloc_fixedR

    LCAO_domain::build_Nonlocal_mu_new(pv,
                                       fsr_dh,
                                       nullptr,
                                       true,
                                       ucell,
                                       orb,
                                       *(two_center_bundle.overlap_orb_beta),
                                       &grid);

    sparse_format::cal_dSTN_R(ucell, pv, HS_Arrays, fsr_dh, grid, orb.cutoffs(), current_spin, sparse_thr, nspin, npol);

    if(nspin==2)
    {
        const double* vr_eff1
            = v_eff.nc * v_eff.nr > 0 ? &(v_eff(current_spin, 0)) : nullptr;
        if (!gamma_only_local)
        {
            ModuleGint::Gint_dvlocal gint_dvlocal(vr_eff1, nspin, npol);
            gint_dvlocal.cal_dvlocal();
            cal_dvlocal_R_sparse(*gint_dvlocal.get_pvdpRx(),
                                 *gint_dvlocal.get_pvdpRy(),
                                 *gint_dvlocal.get_pvdpRz(),
                                 npol,
                                 current_spin,
                                 PARAM.globalv.nlocal,
                                 sparse_thr,
                                 pv,
                                 ucell,
                                 HS_Arrays);
        }
    }
    return;
}

void sparse_format::set_R_range(std::set<Abfs::Vector3_Order<int>>& all_R_coor, const Grid_Driver& grid)
{
    int RminX = int(-grid.getGlayerX_minus());
    int RminY = int(-grid.getGlayerY_minus());
    int RminZ = int(-grid.getGlayerZ_minus());

    int Rx = grid.getGlayerX() + grid.getGlayerX_minus();
    int Ry = grid.getGlayerY() + grid.getGlayerY_minus();
    int Rz = grid.getGlayerZ() + grid.getGlayerZ_minus();

    for (int ix = 0; ix < Rx; ix++)
    {
        for (int iy = 0; iy < Ry; iy++)
        {
            for (int iz = 0; iz < Rz; iz++)
            {
                Abfs::Vector3_Order<int> temp_R(ix + RminX, iy + RminY, iz + RminZ);
                all_R_coor.insert(temp_R);
            }
        }
    }

    return;
}

void sparse_format::cal_dSTN_R(const UnitCell& ucell,
                               const Parallel_Orbitals& pv,
                               LCAO_HS_Arrays& HS_Arrays,
                               ForceStressArrays& fsr,
                               const Grid_Driver& grid,
                               const std::vector<double>& orb_cutoff,
                               const int& current_spin,
                               const double& sparse_thr,
                               const int nspin,
                               const int npol)
{
    ModuleBase::TITLE("sparse_format", "cal_dSTN_R");

    int index = 0;
    ModuleBase::Vector3<double> dtau, tau1, tau2;
    ModuleBase::Vector3<double> dtau1, dtau2, tau0;

    double temp_value_double = 0.0;
    std::complex<double> temp_value_complex = 0.0;

    for (int T1 = 0; T1 < ucell.ntype; ++T1)
    {
        Atom* atom1 = &ucell.atoms[T1];
        for (int I1 = 0; I1 < atom1->na; ++I1)
        {
            tau1 = atom1->tau[I1];
            grid.Find_atom(ucell, tau1, T1, I1);
            Atom* atom1 = &ucell.atoms[T1];
            const int start = ucell.itiaiw2iwt(T1, I1, 0);

            for (int ad = 0; ad < grid.getAdjacentNum() + 1; ++ad)
            {
                const int T2 = grid.getType(ad);
                const int I2 = grid.getNatom(ad);
                Atom* atom2 = &ucell.atoms[T2];

                tau2 = grid.getAdjacentTau(ad);
                dtau = tau2 - tau1;
                double distance = dtau.norm() * ucell.lat0;
                double rcut = orb_cutoff[T1] + orb_cutoff[T2];

                bool adj = false;

                if (distance < rcut)
                {
                    adj = true;
                }
                else if (distance >= rcut)
                {
                    for (int ad0 = 0; ad0 < grid.getAdjacentNum() + 1; ++ad0)
                    {
                        const int T0 = grid.getType(ad0);

                        tau0 = grid.getAdjacentTau(ad0);
                        dtau1 = tau0 - tau1;
                        dtau2 = tau0 - tau2;

                        double distance1 = dtau1.norm() * ucell.lat0;
                        double distance2 = dtau2.norm() * ucell.lat0;

                        double rcut1 = orb_cutoff[T1] + ucell.infoNL->get_rcut_max(T0);
                        double rcut2 = orb_cutoff[T2] + ucell.infoNL->get_rcut_max(T0);

                        if (distance1 < rcut1 && distance2 < rcut2)
                        {
                            adj = true;
                            break;
                        }
                    }
                }

                if (adj)
                {
                    const int start2 = ucell.itiaiw2iwt(T2, I2, 0);

                    Abfs::Vector3_Order<int> dR(grid.getBox(ad).x, grid.getBox(ad).y, grid.getBox(ad).z);

                    for (int ii = 0; ii < atom1->nw * npol; ii++)
                    {
                        const int iw1_all = start + ii;
                        const int mu = pv.global2local_row(iw1_all);

                        if (mu < 0)
                        {
                            continue;
                        }

                        for (int jj = 0; jj < atom2->nw * npol; jj++)
                        {
                            int iw2_all = start2 + jj;
                            const int nu = pv.global2local_col(iw2_all);

                            if (nu < 0)
                            {
                                continue;
                            }

                            if (nspin != 4)
                            {
                                temp_value_double = fsr.DHloc_fixedR_x[index];
                                if (std::abs(temp_value_double) > sparse_thr)
                                {
                                    HS_Arrays.dHRx_sparse[current_spin][dR][iw1_all][iw2_all] = temp_value_double;
                                }
                                temp_value_double = fsr.DHloc_fixedR_y[index];
                                if (std::abs(temp_value_double) > sparse_thr)
                                {
                                    HS_Arrays.dHRy_sparse[current_spin][dR][iw1_all][iw2_all] = temp_value_double;
                                }
                                temp_value_double = fsr.DHloc_fixedR_z[index];
                                if (std::abs(temp_value_double) > sparse_thr)
                                {
                                    HS_Arrays.dHRz_sparse[current_spin][dR][iw1_all][iw2_all] = temp_value_double;
                                }
                            }
                            else
                            {
                                ModuleBase::WARNING_QUIT("cal_dSTN_R", "nspin=4 with SOC is not supported yet.");
                            }
                            ++index;
                        }
                    }
                }
            }
        }
    }

    return;
}

void sparse_format::destroy_dH_R_sparse(LCAO_HS_Arrays& HS_Arrays, const int nspin)
{
    ModuleBase::TITLE("LCAO_domain", "destroy_dH_R_sparse");

    if (nspin != 4)
    {
        std::map<Abfs::Vector3_Order<int>, std::map<size_t, std::map<size_t, double>>> empty_dHRx_sparse_up;
        std::map<Abfs::Vector3_Order<int>, std::map<size_t, std::map<size_t, double>>> empty_dHRx_sparse_down;
        std::map<Abfs::Vector3_Order<int>, std::map<size_t, std::map<size_t, double>>> empty_dHRy_sparse_up;
        std::map<Abfs::Vector3_Order<int>, std::map<size_t, std::map<size_t, double>>> empty_dHRy_sparse_down;
        std::map<Abfs::Vector3_Order<int>, std::map<size_t, std::map<size_t, double>>> empty_dHRz_sparse_up;
        std::map<Abfs::Vector3_Order<int>, std::map<size_t, std::map<size_t, double>>> empty_dHRz_sparse_down;

        HS_Arrays.dHRx_sparse[0].swap(empty_dHRx_sparse_up);
        HS_Arrays.dHRx_sparse[1].swap(empty_dHRx_sparse_down);
        HS_Arrays.dHRy_sparse[0].swap(empty_dHRy_sparse_up);
        HS_Arrays.dHRy_sparse[1].swap(empty_dHRy_sparse_down);
        HS_Arrays.dHRz_sparse[0].swap(empty_dHRz_sparse_up);
        HS_Arrays.dHRz_sparse[1].swap(empty_dHRz_sparse_down);
    }
    else
    {
        std::map<Abfs::Vector3_Order<int>, std::map<size_t, std::map<size_t, std::complex<double>>>>
            empty_dHRx_soc_sparse;
        std::map<Abfs::Vector3_Order<int>, std::map<size_t, std::map<size_t, std::complex<double>>>>
            empty_dHRy_soc_sparse;
        std::map<Abfs::Vector3_Order<int>, std::map<size_t, std::map<size_t, std::complex<double>>>>
            empty_dHRz_soc_sparse;

        HS_Arrays.dHRx_soc_sparse.swap(empty_dHRx_soc_sparse);
        HS_Arrays.dHRy_soc_sparse.swap(empty_dHRy_soc_sparse);
        HS_Arrays.dHRz_soc_sparse.swap(empty_dHRz_soc_sparse);
    }

    return;
}
