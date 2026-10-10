#ifndef ABACUS_SOURCE_LCAO_MODULE_LR_PULAY_HC_H
#define ABACUS_SOURCE_LCAO_MODULE_LR_PULAY_HC_H
#include "source_basis/module_nao/two_center_bundle.h"
#include "source_estate/module_dm/density_matrix.h"
#include "source_cell/unitcell.h"
#include "source_lcao/module_lr/utils/lr_util.h"
#include "source_lcao/module_lr/potentials/pot_lr_base.h"
#include "source_hamilt/module_gint/gint_interface.h"
#include "source_base/parallel_reduce.h"

namespace PulayForceStress
{
// add stess later
/// for 2-center-integration terms, provided HS derivatives
    template<typename TK>
ModuleBase::matrix cal_pulay_fs(
    const module_dm::DensityMatrix<TK, double>& dm,  ///< [in] density matrix or energy density matrix
    const UnitCell& ucell,  ///< [in] unit cell
    const std::vector<hamilt::HContainer<double>>& dHS,  ///< [in] dHS x, y, z, for force
    const double& factor_force = 1.0)
{
    ModuleBase::matrix f(ucell.nat, 3);
    const Parallel_Orbitals& pv = *dHS[0].get_paraV();
    const int& npol = ucell.get_npol();
    const int nspin_dmr = dm.get_dmr_vec().size();
    for (int ixyz = 0;ixyz < 3;++ixyz)
    {
        for (int iat0 = 0;iat0 < ucell.nat;++iat0)
        {
            for (int iat1 = 0;iat1 < ucell.nat;++iat1)
            {
                hamilt::AtomPair<double>* ap = dHS[ixyz].find_pair(iat0, iat1);
                if (ap)
                {
                    for (int iR = 0;iR < ap->get_R_size();++iR)
                    {
                        const ModuleBase::Vector3<int>& R = ap->get_R_index(iR);
                        hamilt::BaseMatrix<double>& mat_dhs = ap->get_HR_values(R.x, R.y, R.z);
                        std::vector<const hamilt::BaseMatrix<double>*> mat_dmr;
                        for (int is = 0; is < nspin_dmr; ++is)
                        {
                            mat_dmr.push_back(dm.get_dmr_ptr(is + 1)->find_matrix(iat0, iat1, R.x, R.y, R.z));
                        }

                        for (int mu = 0; mu < pv.get_nrow_atom(iat0); mu += npol)
                        {
                            for (int nu = 0; nu < pv.get_ncol_atom(iat1); nu += npol)
                            {
                                double dm2d = 0.0;
                                for (int is = 0; is < nspin_dmr; ++is) { dm2d += mat_dmr[is]->get_value(mu, nu); }
                                f(iat0, ixyz) += dm2d * factor_force * 2.0 * mat_dhs.get_value(mu, nu);
                            }
                        }
                    }
                }
            }
        }
    }
    Parallel_Reduce::reduce_all(f.c, f.nr * f.nc);
    return f;
}

/// for grid-integration terms
template<typename TK>
ModuleBase::matrix cal_pulay_fs(
    const module_dm::DensityMatrix<TK, double>& dm,  ///< [in] density matrix or energy density matrix
    const UnitCell& ucell,  ///< [in] unit cell
    const LR::PotLRBase* pot ///< [in] potential on grid
)
{
    ModuleBase::matrix force(ucell.nat, 3);
    ModuleBase::matrix stress_tmp(3, 3);

    // The LR kernel yields a single potential channel, so the grid integrals run with
    // nspin=1 on the first density-matrix channel -- the same thing the old
    // `Gint_inout(is=0, ...)` calls did. `dm` may still carry two channels (the
    // `test_force` path feeds in the ground-state DM and doubles the result afterwards);
    // ModuleGint reads only the leading `nspin` entries of `dm_vec`, but `vr_eff` is
    // indexed up to `nspin`, so passing the DM's spin count here would read past
    // `p_vr_hxc` and segfault.
    constexpr int nspin_gint = 1;

    // 1. dm->rho
    double** rho;
    const int& nrxx = pot->nrxx;
    LR_Util::_allocate_2order_nested_ptr(rho, nspin_gint, nrxx);
    ModuleBase::GlobalFunc::ZEROS(rho[0], nrxx);
    ModuleGint::cal_gint_rho(dm.get_dmr_vec(), nspin_gint, rho, false);

    // 2. v_hxc = f_hxc * rho
    ModuleBase::matrix vr_hxc(1, nrxx);   //grid
    pot->cal_v_eff(rho, ucell, vr_hxc);
    LR_Util::_deallocate_2order_nested_ptr(rho, nspin_gint);

    // 3. v(r) -> force
    // An empty local FFT slab has no (0, 0) element. Pass the storage pointer
    // directly; the grid integrator has no local points to read on that rank.
    const std::vector<const double*> p_vr_hxc(nspin_gint, vr_hxc.c);
    ModuleGint::cal_gint_fvl(nspin_gint, p_vr_hxc, dm.get_dmr_vec(), /*isforce=*/true, /*isstress=*/false, &force, &stress_tmp);
    // `cal_gint_fvl` only sums the grid points (and their atom pairs) this rank's share of the
    // real-space FFT box touches; core ABACUS always follows it with this same reduction (see
    // `force_stress_lcao.cpp`'s `Parallel_Reduce::reduce_pool` right after its own grid-based
    // `cal_pulay_fs` call) -- omitted here, so this was silently wrong for nprocs>1.
    Parallel_Reduce::reduce_pool(force.c, force.nr * force.nc);
    return force;
}

/// @brief Open-shell (spin-unrestricted) counterpart of the grid `cal_pulay_fs` above.
///
/// The LR kernel mixes the two channels,
///     $v_\sigma=\sum_{\sigma'}f^{\sigma\sigma'}\rho^1_{\sigma'}$ (+ the full Hartree),
/// so the density has to be built per channel, the potential accumulated over the inner
/// spin, and only then contracted with the density matrix of the outer spin. `PotHxcLR`
/// with `SpinType::S2_updown` selects the $(\sigma,\sigma')$ component via `ispin_op`.
template<typename TK>
ModuleBase::matrix cal_pulay_fs_openshell(
    const module_dm::DensityMatrix<TK, double>& dm,  ///< [in] 2-channel density matrix
    const UnitCell& ucell,
    const LR::PotLRBase* pot)
{
    ModuleBase::matrix force(ucell.nat, 3);
    ModuleBase::matrix stress_tmp(3, 3);
    constexpr int nspin_dm = 2;
    assert(dm.get_dmr_vec().size() == nspin_dm);

    // 1. dm -> rho, one channel each
    double** rho;
    const int& nrxx = pot->nrxx;
    LR_Util::_allocate_2order_nested_ptr(rho, nspin_dm, nrxx);
    for (int is = 0; is < nspin_dm; ++is) { ModuleBase::GlobalFunc::ZEROS(rho[is], nrxx); }
    ModuleGint::cal_gint_rho(dm.get_dmr_vec(), nspin_dm, rho, false);

    // 2. $v_\sigma=\sum_{\sigma'}f^{\sigma\sigma'}\rho_{\sigma'}$
    std::vector<ModuleBase::matrix> vr_hxc(nspin_dm, ModuleBase::matrix(1, nrxx));
    for (int sl = 0; sl < nspin_dm; ++sl)
    {
        for (int sr = 0; sr < nspin_dm; ++sr)
        {
            double* rho_in[1] = { rho[sr] };
            pot->cal_v_eff(rho_in, ucell, vr_hxc[sl], { sl, sr });
        }
    }
    LR_Util::_deallocate_2order_nested_ptr(rho, nspin_dm);

    // 3. v(r) -> force, summed over the outer spin by `cal_gint_fvl`
    std::vector<const double*> p_vr_hxc(nspin_dm);
    // Keep empty-grid ranks in the integration and subsequent pool reduction.
    for (int is = 0; is < nspin_dm; ++is) { p_vr_hxc[is] = vr_hxc[is].c; }
    ModuleGint::cal_gint_fvl(nspin_dm, p_vr_hxc, dm.get_dmr_vec(), /*isforce=*/true, false, &force, &stress_tmp);
    Parallel_Reduce::reduce_pool(force.c, force.nr * force.nc);   // see the closed-shell overload above
    return force;
}
}

#endif // ABACUS_SOURCE_LCAO_MODULE_LR_PULAY_HC_H
