#include "force_stress_terms.h"

#include <iostream>

#include "source_base/parallel_reduce.h"
#include "source_base/timer.h"
#include "source_base/tool_quit.h"
#include "source_basis/module_nao/two_center_bundle.h"
#include "source_estate/module_pot/efield.h"
#include "source_estate/module_pot/gatefield.h"
#include "source_estate/module_pot/h_tddft_pw.h"
#include "source_hamilt/module_surchem/surchem.h"
#include "source_hamilt/module_vdw/vdw.h"
#include "source_io/module_parameter/parameter.h"
#include "source_lcao/module_dftu/dftu_nao_adj.h"
#include "source_lcao/module_dftu/dftu_nao_fs_k.h"
#include "source_lcao/module_dftu/dftu_nao_fs_r.h"
#include "source_pw/module_pwdft/dftu_base.h"
#ifdef __MLALGO
#include "source_lcao/module_deepks/deepks_force.h"
#endif

namespace LCAO_domain
{

namespace
{
// Copy the cached vdW force / stress into the parts containers.
void copy_vdw_terms(const vdw::VdwResult* vdw_result,
                    const UnitCell& ucell,
                    const bool isforce,
                    const bool isstress,
                    LCAOForceParts& parts,
                    LCAOStressParts& sparts)
{
    if (vdw_result == nullptr)
    {
        return;
    }
    if (isforce)
    {
        if (!vdw_result->has_force || vdw_result->force.size() != static_cast<std::size_t>(ucell.nat))
        {
            ModuleBase::WARNING_QUIT("Force_Stress_LCAO::getForceStress",
                                     "The cached vdW force is unavailable or has an invalid size.");
        }
        parts.force_vdw.create(ucell.nat, 3);
        for (int iat = 0; iat < ucell.nat; ++iat)
        {
            parts.force_vdw(iat, 0) = vdw_result->force[iat].x;
            parts.force_vdw(iat, 1) = vdw_result->force[iat].y;
            parts.force_vdw(iat, 2) = vdw_result->force[iat].z;
        }
    }
    if (isstress)
    {
        if (!vdw_result->has_stress)
        {
            ModuleBase::WARNING_QUIT("Force_Stress_LCAO::getForceStress",
                                     "The cached vdW stress is unavailable.");
        }
        sparts.stress_vdw = vdw_result->stress.to_matrix();
    }
}

// Compute the external-field force terms (E-field, rt-TDDFT, gate, solvation).
void cal_external_field_forces(UnitCell& ucell,
                               surchem& solvent,
                               ModulePW::PW_Basis* rhopw,
                               const pseudopot_cell_vl& locpp,
                               LCAOForceParts& parts)
{
    //! forces from E-field
    if (PARAM.inp.efield_flag)
    {
        parts.fefield.create(ucell.nat, 3);
        elecstate::Efield::compute_force(ucell, parts.fefield);
    }

    //! atomic forces from E-field of rt-TDDFT
    if (PARAM.inp.esolver_type == "tddft")
    {
        parts.fefield_tddft.create(ucell.nat, 3);
        elecstate::H_TDDFT_pw::compute_force(ucell, parts.fefield_tddft);
    }

    //! atomic forces from gate field
    if (PARAM.inp.gate_flag)
    {
        parts.fgate.create(ucell.nat, 3);
        elecstate::Gatefield::compute_force(ucell, parts.fgate);
    }

    //! atomic forces from implicit solvation model
    if (PARAM.inp.imp_sol)
    {
        parts.fsol.create(ucell.nat, 3);
        solvent.cal_force_sol(ucell, rhopw, locpp.vloc, PARAM.inp.nspin, parts.fsol);
    }
}
} // namespace

void cal_vdw_fields_fs(const vdw::VdwResult* vdw_result,
                       UnitCell& ucell,
                       surchem& solvent,
                       ModulePW::PW_Basis* rhopw,
                       const pseudopot_cell_vl& locpp,
                       const bool isforce,
                       const bool isstress,
                       LCAOForceParts& parts,
                       LCAOStressParts& sparts)
{
    //! forces and stress from vdw
    //  Peize Lin add 2014-04-04, update 2021-03-09
    //  jiyy add 2019-05-18, update 2021-05-02
    copy_vdw_terms(vdw_result, ucell, isforce, isstress, parts, sparts);

    //! external-field forces
    if (isforce)
    {
        cal_external_field_forces(ucell, solvent, rhopw, locpp, parts);
    }
}

template <typename T>
void cal_exx_fs(const UnitCell& ucell,
                const bool isforce,
                const bool isstress,
                const Exx_Info& exx_info,
                Exx_NAO<T>& exx_nao,
                LCAOForceParts& parts,
                LCAOStressParts& sparts)
{
#ifdef __EXX
    bool cal_exx = exx_info.info_global.cal_exx;
    bool real_number = exx_info.info_ri.real_number;
    double hybrid_alpha = exx_info.info_global.hybrid_alpha;

    if (cal_exx)
    {
        if (isforce || isstress)
        {
            std::cout << " >> NOTICE: calculating EXX force/stress, which may be time-consuming" << std::endl;
        }
        if (isforce)
        {
            if (real_number)
            {
                exx_nao.exd->cal_exx_force(ucell.nat);
                parts.force_exx = hybrid_alpha * exx_nao.exd->get_force();
            }
            else
            {
                exx_nao.exc->cal_exx_force(ucell.nat);
                parts.force_exx = hybrid_alpha * exx_nao.exc->get_force();
            }
        }
        if (isstress)
        {
            if (real_number)
            {
                exx_nao.exd->cal_exx_stress(ucell.omega, ucell.lat0);
                sparts.stress_exx = hybrid_alpha * exx_nao.exd->get_stress();
            }
            else
            {
                exx_nao.exc->cal_exx_stress(ucell.omega, ucell.lat0);
                sparts.stress_exx = hybrid_alpha * exx_nao.exc->get_stress();
            }
        }
    }
#endif
}

template <typename T>
void cal_dftu_fs(UnitCell& ucell,
                 const Grid_Driver& gd,
                 Parallel_Orbitals& pv,
                 const LCAO_Orbitals& orb,
                 const K_Vectors& kv,
                 LCAO_domain::Setup_DM<T>& dmat,
                 const TwoCenterBundle& two_center_bundle,
                 Plus_U_Base& dftu,
                 const bool isforce,
                 const bool isstress,
                 LCAOForceParts& parts,
                 LCAOStressParts& sparts)
{
    //! atomic forces from DFT+U (Quxin version)

    if (PARAM.inp.dft_plus_u) // Quxin add for DFT+U on 20201029
    {
        if (isforce)
        {
            parts.force_u.create(ucell.nat, 3);
        }
        if (isstress)
        {
            sparts.stress_u.create(3, 3);
        }
        if (PARAM.inp.dft_plus_u == 2)
        {
            // The legacy dft_plus_u==2 force/stress path is currently broken.
            //
            // Background: Plus_U::force_stress relies on ForceStressArrays
            // members DSloc_x/y/z (gamma_only) or DSloc_Rx/Ry/Rz (multik)
            // and DH_r being pre-allocated and filled with dS/dR data by the
            // main force flow (formerly ForceLcaoGamma::ftable). The DFT+U
            // step 2 refactor (commit 70c54c9d5a, 2026-01-23) removed the
            // main-flow ForceStressArrays because the operator-based force
            // calculation no longer needs it, but the legacy dft_plus_u==2
            // path still depends on it. The local fsr_dftu below is declared
            // without allocating those arrays, so any call into
            // cal_force_gamma / cal_stress_gamma / folding_matrix_k would
            // pass nullptr to pdgemm_ and crash with SIGSEGV.
            //
            // Until the legacy path is restored or re-implemented, we
            // explicitly reject dft_plus_u==2 with cal_force or cal_stress
            // enabled. SCF-only runs (no force/stress) are unaffected
            // because the energy is computed in cal_energy_correction,
            // which does not touch DSloc arrays. Use dft_plus_u=1 for
            // force/stress calculations.
            if (isforce || isstress)
            {
                ModuleBase::WARNING_QUIT("Force_Stress_LCAO::getForceStress",
                    "dft_plus_u==2 with cal_force or cal_stress is currently broken; "
                    "please use dft_plus_u=1 instead. See notes in source/source_lcao/force_stress_lcao.cpp.");
            }
            ForceStressArrays fsr_dftu;
            std::vector<std::vector<double>>* dmk_d = nullptr;
            std::vector<std::vector<std::complex<double>>>* dmk_c = nullptr;
            assign_dmk_ptr<T>(dmat.dm, dmk_d, dmk_c);
            DFTU_LCAO::DftuFsEnv dftu_fs_env(dftu, ucell, gd, pv, fsr_dftu,
                                             orb.cutoffs(), PARAM.inp.ks_solver);
            DFTU_LCAO::force_stress(dftu_fs_env, isforce, isstress,
                                    dmk_d, dmk_c, parts.force_u, sparts.stress_u, kv,
                                    PARAM.globalv.gamma_only_local);
        }
        else
        {
            // Build DFT+U force/stress inputs directly without constructing a
            // full DFTU operator (hsk/hR are irrelevant for this path).
            std::vector<AdjacentAtomInfo> adjs_all = DFTU_LCAO::build_adjacent_atoms(
                &ucell, &dftu, &gd, orb.cutoffs(), PARAM.inp.onsite_radius);

            // The DensityMatrix holds nspin_dm = (nspin==2 ? 2 : 1) real-space DMR
            // channels: nspin=4 (non-collinear) packs all four Pauli components
            // into a single complex DMR, so only one channel exists (cf. setup_dm.cpp
            // and the is0 = nspin==2 ? is : 0 indexing in cal_for/str_IJR_nao_r).
            const int nspin_dm = (PARAM.inp.nspin == 2) ? 2 : 1;
            std::vector<const hamilt::HContainer<double>*> dmR_tmp(nspin_dm, nullptr);
            for (int is = 0; is < nspin_dm; ++is)
            {
                dmR_tmp[is] = dmat.dm->get_DMR_pointer(is + 1);
            }

            DFTU_LCAO::cal_fs_nao_r(&ucell, &dftu,
                                    two_center_bundle.overlap_orb_onsite.get(),
                                    PARAM.inp.nspin,
                                    adjs_all, dmR_tmp,
                                    isforce, isstress, parts.force_u, sparts.stress_u);
        }
    }
}

template <typename T>
void cal_deepks_fs(const UnitCell& ucell,
                   const Grid_Driver& gd,
                   Parallel_Orbitals& pv,
                   const LCAO_Orbitals& orb,
                   const K_Vectors& kv,
                   const bool isforce,
                   const bool isstress,
                   Setup_DeePKS<T>& deepks,
                   LCAOForceParts& parts,
                   LCAOStressParts& sparts)
{
    // Handle DeePKS forces if enabled
#ifdef __MLALGO
    if (PARAM.inp.deepks_scf)
    {
        const int nks = (PARAM.inp.nspin == 1 || PARAM.inp.nspin == 2) ? 1 : kv.get_nks();
        if (PARAM.globalv.gamma_only_local)
        {
            DeePKS_domain::cal_f_delta<double>(
                ucell,
                orb,
                gd,
                pv,
                nks,
                deepks.ld.deepks_param,
                kv.kvec_d,
                deepks.ld.phialpha,
                parts.fvnl_dalpha,
                isstress,
                sparts.svnl_dalpha,
                deepks.ld.dm_r,
                deepks.ld.gedm,
                (PARAM.inp.nspin == 2 && !PARAM.inp.deepks_equiv) ? deepks.ld.dm_r_mag : nullptr,
                (PARAM.inp.nspin == 2 && !PARAM.inp.deepks_equiv) ? deepks.ld.gedm_mag : nullptr);
        }
        else
        {
            DeePKS_domain::cal_f_delta<std::complex<double>>(
                ucell,
                orb,
                gd,
                pv,
                nks,
                deepks.ld.deepks_param,
                kv.kvec_d,
                deepks.ld.phialpha,
                parts.fvnl_dalpha,
                isstress,
                sparts.svnl_dalpha,
                deepks.ld.dm_r,
                deepks.ld.gedm,
                (PARAM.inp.nspin == 2 && !PARAM.inp.deepks_equiv) ? deepks.ld.dm_r_mag : nullptr,
                (PARAM.inp.nspin == 2 && !PARAM.inp.deepks_equiv) ? deepks.ld.gedm_mag : nullptr);
        }

        if (isforce)
        {
            Parallel_Reduce::reduce_pool(parts.fvnl_dalpha.c, parts.fvnl_dalpha.nr * parts.fvnl_dalpha.nc);
        }
        if (isstress)
        {
            Parallel_Reduce::reduce_pool(sparts.svnl_dalpha.c, sparts.svnl_dalpha.nr * sparts.svnl_dalpha.nc);
        }
    }
#endif
}

// Explicit instantiation for the two electronic types used by ABACUS.
template void cal_exx_fs<double>(const UnitCell&, const bool, const bool, const Exx_Info&, Exx_NAO<double>&,
                                 LCAOForceParts&, LCAOStressParts&);
template void cal_exx_fs<std::complex<double>>(const UnitCell&, const bool, const bool, const Exx_Info&,
                                               Exx_NAO<std::complex<double>>&, LCAOForceParts&, LCAOStressParts&);

template void cal_dftu_fs<double>(UnitCell&, const Grid_Driver&, Parallel_Orbitals&, const LCAO_Orbitals&,
                                  const K_Vectors&, LCAO_domain::Setup_DM<double>&, const TwoCenterBundle&,
                                  Plus_U_Base&, const bool, const bool, LCAOForceParts&, LCAOStressParts&);
template void cal_dftu_fs<std::complex<double>>(UnitCell&, const Grid_Driver&, Parallel_Orbitals&,
                                                const LCAO_Orbitals&, const K_Vectors&,
                                                LCAO_domain::Setup_DM<std::complex<double>>&, const TwoCenterBundle&,
                                                Plus_U_Base&, const bool, const bool, LCAOForceParts&,
                                                LCAOStressParts&);

template void cal_deepks_fs<double>(const UnitCell&, const Grid_Driver&, Parallel_Orbitals&, const LCAO_Orbitals&,
                                    const K_Vectors&, const bool, const bool, Setup_DeePKS<double>&, LCAOForceParts&,
                                    LCAOStressParts&);
template void cal_deepks_fs<std::complex<double>>(const UnitCell&, const Grid_Driver&, Parallel_Orbitals&,
                                                  const LCAO_Orbitals&, const K_Vectors&, const bool, const bool,
                                                  Setup_DeePKS<std::complex<double>>&, LCAOForceParts&,
                                                  LCAOStressParts&);

} // namespace LCAO_domain
