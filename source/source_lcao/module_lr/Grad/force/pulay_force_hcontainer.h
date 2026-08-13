#pragma once 
#include "source_basis/module_nao/two_center_bundle.h"
#include "source_estate/module_dm/density_matrix.h"
#include "source_cell/unitcell.h"
#include "source_lcao/module_lr/utils/lr_util.h"
#include "source_lcao/module_lr/potentials/pot_lr_base.h"
#include "source_hamilt/module_gint/gint_interface.h"

namespace PulayForceStress
{
// add stess later
/// for 2-center-integration terms, provided HS derivatives
    template<typename TK>
ModuleBase::matrix cal_pulay_fs(
    const elecstate::DensityMatrix<TK, double>& dm,  ///< [in] density matrix or energy density matrix
    const UnitCell& ucell,  ///< [in] unit cell
    const std::vector<hamilt::HContainer<double>>& dHS,  ///< [in] dHS x, y, z, for force
    const double& factor_force = 1.0)
{
    ModuleBase::matrix f(ucell.nat, 3);
    const Parallel_Orbitals& pv = *dHS[0].get_paraV();
    const int& npol = ucell.get_npol();
    const int nspin_dmr = dm.get_DMR_vector().size();
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
                            mat_dmr.push_back(dm.get_DMR_pointer(is + 1)->find_matrix(iat0, iat1, R.x, R.y, R.z));
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
    return f;
}

/// for grid-integration terms
template<typename TK>
ModuleBase::matrix cal_pulay_fs(
    const elecstate::DensityMatrix<TK, double>& dm,  ///< [in] density matrix or energy density matrix
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
    ModuleGint::cal_gint_rho(dm.get_DMR_vector(), nspin_gint, rho, false);

    // 2. v_hxc = f_hxc * rho
    ModuleBase::matrix vr_hxc(1, nrxx);   //grid
    pot->cal_v_eff(rho, ucell, vr_hxc);
    LR_Util::_deallocate_2order_nested_ptr(rho, nspin_gint);

    // 3. v(r) -> force
    const std::vector<const double*> p_vr_hxc(nspin_gint, &vr_hxc(0, 0));
    ModuleGint::cal_gint_fvl(nspin_gint, p_vr_hxc, dm.get_DMR_vector(), /*isforce=*/true, /*isstress=*/false, &force, &stress_tmp);
    return force;
}
}
 