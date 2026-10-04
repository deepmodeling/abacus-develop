#include "force_stress_assemble.h"

#include <iomanip>

#include "source_base/constants.h"
#include "source_base/global_function.h"
#include "source_base/global_variable.h"
#include "source_base/tool_quit.h"
#include "source_cell/module_symmetry/symmetry.h"
#include "source_cell/unitcell.h"
#include "source_hamilt/module_vdw/vdw.h"
#include "source_io/module_output/output_log.h"
#include "source_io/module_parameter/parameter.h"
#include "source_lcao/setup_deepks.h" // DeePKS_domain::write_forces/write_stress
#ifdef __MLALGO
#include "source_lcao/module_deepks/lcao_deepks.h"
#include "source_lcao/module_deepks/lcao_deepks_io.h"
#endif

namespace LCAO_domain
{

namespace
{
// Accumulate every active force term into fcs for one Cartesian component.
// Splitting this out of assemble_print_force keeps the latter's cyclomatic
// complexity low; the branch conditions mirror the terms that were computed.
void sum_force_terms(const int iat,
                     const int i,
                     const vdw::VdwResult* vdw_result,
                     const Exx_Info& exx_info,
                     const LCAOForceParts& parts,
                     ModuleBase::matrix& fcs)
{
    fcs(iat, i) += parts.foverlap(iat, i) + parts.ftvnl_dphi(iat, i) + parts.fvnl_dbeta(iat, i) + parts.fvl_dphi(iat, i)
                   + parts.fvl_dvl(iat, i) // derivative of local potential force (pw)
                   + parts.fewalds(iat, i) // ewald force (pw)
                   + parts.fcc(iat, i)     // nonlinear core correction force (pw)
                   + parts.fscc(iat, i)    // self consistent corretion force (pw)
                   + parts.fpothybrid(iat, i); // pulay force for hybrid gauge rt-tddft

    // Force contribution from DFT+U, Quxin add on 20201029
    if (PARAM.inp.dft_plus_u)
    {
        fcs(iat, i) += parts.force_u(iat, i);
    }
    if (PARAM.inp.sc_mag_switch)
    {
        fcs(iat, i) += parts.force_dspin(iat, i);
    }
#ifdef __EXX
    // Force contribution from exx
    if (exx_info.info_global.cal_exx)
    {
        fcs(iat, i) += parts.force_exx(iat, i);
    }
#endif
    // VDW force of vdwd2 or vdwd3
    if (vdw_result != nullptr)
    {
        fcs(iat, i) += parts.force_vdw(iat, i);
    }
    // E-field force
    if (PARAM.inp.efield_flag)
    {
        fcs(iat, i) += parts.fefield(iat, i);
    }
    // E-field force of tddft
    if (PARAM.inp.esolver_type == "tddft")
    {
        fcs(iat, i) += parts.fefield_tddft(iat, i);
    }
    // Gate field force
    if (PARAM.inp.gate_flag)
    {
        fcs(iat, i) += parts.fgate(iat, i);
    }
    // implicit solvation model
    if (PARAM.inp.imp_sol)
    {
        fcs(iat, i) += parts.fsol(iat, i);
    }
#ifdef __MLALGO
    // mohan add 2021-08-04
    if (PARAM.inp.deepks_scf)
    {
        fcs(iat, i) += parts.fvnl_dalpha(iat, i);
    }
#endif
}

// Accumulate every active stress term into scs for one tensor component.
void sum_stress_terms(const int i,
                      const int j,
                      const vdw::VdwResult* vdw_result,
                      const Exx_Info& exx_info,
                      const LCAOStressParts& sparts,
                      ModuleBase::matrix& scs)
{
    scs(i, j) += sparts.soverlap(i, j) + sparts.stvnl_dphi(i, j) + sparts.svnl_dbeta(i, j) + sparts.svl_dphi(i, j)
                 + sparts.sigmadvl(i, j)  // derivative of local potential stress (pw)
                 + sparts.sigmaewa(i, j)  // ewald stress (pw)
                 + sparts.sigmacc(i, j)   // nonlinear core correction stress (pw)
                 + sparts.sigmaxc(i, j)   // exchange corretion stress
                 + sparts.sigmahar(i, j); // hartree stress

    // VDW stress from linpz and jiyy
    if (vdw_result != nullptr)
    {
        scs(i, j) += sparts.stress_vdw(i, j);
    }
    // DFT plus U stress from qux
    if (PARAM.inp.dft_plus_u)
    {
        scs(i, j) += sparts.stress_u(i, j);
    }
    if (PARAM.inp.sc_mag_switch)
    {
        scs(i, j) += sparts.stress_dspin(i, j);
    }
#ifdef __EXX
    // Stress contribution from exx
    if (exx_info.info_global.cal_exx)
    {
        scs(i, j) += sparts.stress_exx(i, j);
    }
#endif
#ifdef __MLALGO
    if (PARAM.inp.deepks_scf)
    {
        scs(i, j) += sparts.svnl_dalpha(i, j);
    }
#endif
}
} // namespace

namespace
{
// Print every individual force term (test output only, istestf == true).
void print_force_parts(const UnitCell& ucell,
                       const vdw::VdwResult* vdw_result,
                       const LCAOForceParts& parts)
{
    const int nat = ucell.nat;
    ModuleBase::matrix ftvnl;
    ftvnl.create(nat, 3);
    for (int iat = 0; iat < nat; iat++)
    {
        for (int i = 0; i < 3; i++)
        {
            ftvnl(iat, i) = parts.ftvnl_dphi(iat, i) + parts.fvnl_dbeta(iat, i);
        }
    }

    GlobalV::ofs_running << "\n PARTS OF FORCE: " << std::endl;
    GlobalV::ofs_running << std::setiosflags(std::ios::showpos);
    GlobalV::ofs_running << std::setiosflags(std::ios::fixed) << std::setprecision(8) << std::endl;
    ModuleIO::print_force(GlobalV::ofs_running, ucell, "OVERLAP    FORCE", parts.foverlap, false);
    ModuleIO::print_force(GlobalV::ofs_running, ucell, "TVNL_DPHI  force", parts.ftvnl_dphi, false);
    ModuleIO::print_force(GlobalV::ofs_running, ucell, "VNL_DBETA  force", parts.fvnl_dbeta, false);
    ModuleIO::print_force(GlobalV::ofs_running, ucell, "T_VNL      FORCE", ftvnl, false);
    ModuleIO::print_force(GlobalV::ofs_running, ucell, "VL_dPHI    FORCE", parts.fvl_dphi, false);
    ModuleIO::print_force(GlobalV::ofs_running, ucell, "VL_dVL     FORCE", parts.fvl_dvl, false);
    ModuleIO::print_force(GlobalV::ofs_running, ucell, "EWALD      FORCE", parts.fewalds, false);
    ModuleIO::print_force(GlobalV::ofs_running, ucell, "NLCC       FORCE", parts.fcc, false);
    ModuleIO::print_force(GlobalV::ofs_running, ucell, "SCC        FORCE", parts.fscc, false);
    if (PARAM.inp.efield_flag)
    {
        ModuleIO::print_force(GlobalV::ofs_running, ucell, "EFIELD     FORCE", parts.fefield, false);
    }
    if (PARAM.inp.esolver_type == "tddft")
    {
        ModuleIO::print_force(GlobalV::ofs_running, ucell, "EFIELD_TDDFT     FORCE", parts.fefield_tddft, false);
    }
    if (PARAM.inp.gate_flag)
    {
        ModuleIO::print_force(GlobalV::ofs_running, ucell, "GATEFIELD     FORCE", parts.fgate, false);
    }
    if (PARAM.inp.imp_sol)
    {
        ModuleIO::print_force(GlobalV::ofs_running, ucell, "IMP_SOL     FORCE", parts.fsol, false);
    }
    if (vdw_result != nullptr)
    {
        ModuleIO::print_force(GlobalV::ofs_running, ucell, "VDW        FORCE", parts.force_vdw, false);
    }
    if (PARAM.inp.dft_plus_u)
    {
        ModuleIO::print_force(GlobalV::ofs_running, ucell, "DFT+U      FORCE", parts.force_u, false);
    }
    if (PARAM.inp.sc_mag_switch)
    {
        ModuleIO::print_force(GlobalV::ofs_running, ucell, "DeltaSpin  FORCE", parts.force_dspin, false);
    }
#ifdef __MLALGO
    // caoyu add 2021-06-03
    if (PARAM.inp.deepks_scf)
    {
        ModuleIO::print_force(GlobalV::ofs_running, ucell, "DeePKS     FORCE", parts.fvnl_dalpha, true);
    }
#endif
}

// Print the per-atom flag table and zero out sub-threshold force components
// (test output only, istestf == true).
void print_force_invalid_table(const UnitCell& ucell,
                               const double force_threshold,
                               ModuleBase::matrix& fcs)
{
    GlobalV::ofs_running << "\n FORCE INVALID TABLE." << std::endl;
    GlobalV::ofs_running << " " << std::setw(8) << "atom" << std::setw(5) << "x" << std::setw(5) << "y"
                         << std::setw(5) << "z" << std::endl;
    for (int iat = 0; iat < ucell.nat; iat++)
    {
        GlobalV::ofs_running << " " << std::setw(8) << iat;
        for (int i = 0; i < 3; i++)
        {
            if (std::abs(fcs(iat, i) * ModuleBase::Ry_to_eV / ModuleBase::BOHR_TO_A) < force_threshold)
            {
                fcs(iat, i) = 0.0;
                GlobalV::ofs_running << std::setw(5) << "1";
            }
            else
            {
                GlobalV::ofs_running << std::setw(5) << "0";
            }
        }
        GlobalV::ofs_running << std::endl;
    }
}
} // namespace

void assemble_print_force(const UnitCell& ucell,
                          const bool istestf,
                          const vdw::VdwResult* vdw_result,
                          const Exx_Info& exx_info,
                          ModuleSymmetry::Symmetry* symm,
                          const std::string& dpks_out_type,
                          const LCAOForceParts& parts,
                          const double force_threshold,
                          ModuleBase::matrix& fcs)
{
    const int nat = ucell.nat;
    //---------------------------------
    // sum all parts of force!
    //---------------------------------
    ModuleBase::Vector3<double> net_force = {0.0, 0.0, 0.0};
    for (int i = 0; i < 3; i++)
    {
        for (int iat = 0; iat < nat; iat++)
        {
            sum_force_terms(iat, i, vdw_result, exx_info, parts, fcs);
        }
    }

    if (PARAM.inp.gate_flag || PARAM.inp.efield_flag)
    {
        GlobalV::ofs_running << "Atomic forces are not shifted if gate_flag or efield_flag == true!" << std::endl;
    }

    // pengfei 2016-12-20
    if (ModuleSymmetry::Symmetry::symm_flag == 1)
    {
        ModuleSymmetry::symmetrize_force_cartesian(symm, nat, ucell.a1, ucell.a2, ucell.a3, fcs);
    }

    // The net force should be evaluated AFTER the symmetrization.
    // With symmetry switched on, the forces assembled above are built from IBZ-reduced
    // quantities and only become physical after the symmetrization, forceSymmetry(). 
    // Force symmetrization is linear, so it commutes with the removal of a
    // uniform shift: the resulting fcs is identical to the previous ordering.
    // Net force is evaluated after symmetrization and before the uniform shift.
    for (int i = 0; i < 3; i++)
    {
        double sum = 0.0;
        for (int iat = 0; iat < nat; iat++)
        {
            sum += fcs(iat, i);
        }
        net_force[i] = sum;
    }
    if (!(PARAM.inp.gate_flag || PARAM.inp.efield_flag))
    {
        ModuleBase::remove_net_force(nat, fcs);
    }

    // compute forces using the DeePKS model
    DeePKS_domain::write_forces(fcs, parts.fvnl_dalpha, dpks_out_type, PARAM.inp);

    if (istestf)
    {
        print_force_parts(ucell, vdw_result, parts);
    }

    GlobalV::ofs_running << std::setiosflags(std::ios::left);

    // this->printforce_total(ry, istestf, fcs);
    ModuleIO::print_force(GlobalV::ofs_running, ucell, "TOTAL-FORCE (eV/Angstrom)", fcs, false);
    net_force*= ModuleBase::Ry_to_eV / ModuleBase::BOHR_TO_A;
    ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running, "Net force vector (eV/Ang)", net_force.x, net_force.y, net_force.z);
    ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running, "Total drift (ev/Ang)", net_force.norm());
    if (istestf)
    {
        print_force_invalid_table(ucell, force_threshold, fcs);
    }
}

void assemble_print_stress(const UnitCell& ucell,
                           const bool istests,
                           const vdw::VdwResult* vdw_result,
                           const Exx_Info& exx_info,
                           ModuleSymmetry::Symmetry* symm,
                           const std::string& dpks_out_type,
                           const LCAOStressParts& sparts,
                           ModuleBase::matrix& scs)
{
    for (int i = 0; i < 3; i++)
    {
        for (int j = 0; j < 3; j++)
        {
            sum_stress_terms(i, j, vdw_result, exx_info, sparts, scs);
        }
    }
    if (ModuleSymmetry::Symmetry::symm_flag == 1)
    {
        symm->symmetrize_mat3(scs, ucell.lat);
    } // end symmetry

    DeePKS_domain::write_stress(scs, sparts.svnl_dalpha, ucell.omega, dpks_out_type, PARAM.inp);

    // print Rydberg stress or not
    bool ry = false;

    // test stress each terms if needed
    if (istests)
    {
        // test
        ModuleBase::matrix svlocal;
        svlocal.create(3, 3);
        ModuleBase::matrix stvnl;
        stvnl.create(3, 3);
        for (int i = 0; i < 3; i++)
        {
            for (int j = 0; j < 3; j++)
            {
                svlocal(i, j) = sparts.svl_dphi(i, j) + sparts.sigmadvl(i, j);
                stvnl(i, j) = sparts.stvnl_dphi(i, j) + sparts.svnl_dbeta(i, j);
            }
        }

        const bool screen = PARAM.inp.test_stress;

        GlobalV::ofs_running << "\n PARTS OF STRESS: " << std::endl;
        GlobalV::ofs_running << std::setiosflags(std::ios::showpos);
        GlobalV::ofs_running << std::setiosflags(std::ios::fixed) << std::setprecision(10) << std::endl;
        ModuleIO::print_stress("OVERLAP  STRESS", sparts.soverlap, screen, ry, GlobalV::ofs_running);
        ModuleIO::print_stress("T        STRESS", sparts.stvnl_dphi, screen, ry, GlobalV::ofs_running);
        ModuleIO::print_stress("VNL      STRESS", sparts.svnl_dbeta, screen, ry, GlobalV::ofs_running);
        ModuleIO::print_stress("T_VNL    STRESS", stvnl, screen, ry, GlobalV::ofs_running);
        ModuleIO::print_stress("VL_dPHI  STRESS", sparts.svl_dphi, screen, ry, GlobalV::ofs_running);
        ModuleIO::print_stress("VL_dVL   STRESS", sparts.sigmadvl, screen, ry, GlobalV::ofs_running);
        ModuleIO::print_stress("HAR      STRESS", sparts.sigmahar, screen, ry, GlobalV::ofs_running);
        ModuleIO::print_stress("EWALD    STRESS", sparts.sigmaewa, screen, ry, GlobalV::ofs_running);
        ModuleIO::print_stress("cc       STRESS", sparts.sigmacc, screen, ry, GlobalV::ofs_running);
        ModuleIO::print_stress("XC       STRESS", sparts.sigmaxc, screen, ry, GlobalV::ofs_running);
        if (vdw_result != nullptr)
        {
            ModuleIO::print_stress("VDW      STRESS", sparts.stress_vdw, screen, ry, GlobalV::ofs_running);
        }
        if (PARAM.inp.dft_plus_u)
        {
            ModuleIO::print_stress("DFTU     STRESS", sparts.stress_u, screen, ry, GlobalV::ofs_running);
        }
        if (PARAM.inp.sc_mag_switch)
        {
            ModuleIO::print_stress("DeltaSpin  STRESS", sparts.stress_dspin, screen, ry, GlobalV::ofs_running);
        }
#ifdef __EXX
        if (exx_info.info_global.cal_exx)
        {
            ModuleIO::print_stress("EXX      STRESS", sparts.stress_exx, screen, ry, GlobalV::ofs_running);
        }
#endif
        ModuleIO::print_stress("TOTAL    STRESS", scs, screen, ry, GlobalV::ofs_running);
    } // end of test

    GlobalV::ofs_running << std::setiosflags(std::ios::left);

    // print total stress
    bool screen_normal = true;
    ModuleIO::print_stress("TOTAL-STRESS", scs, screen_normal, ry, GlobalV::ofs_running);

    double unit_transform = 0.0;
    unit_transform = ModuleBase::RYDBERG_SI / pow(ModuleBase::BOHR_RADIUS_SI, 3) * 1.0e-8;
    double external_stress[3] = {PARAM.inp.press1, PARAM.inp.press2, PARAM.inp.press3};

    for (int i = 0; i < 3; i++)
    {
        scs(i, i) -= external_stress[i] / unit_transform;
    }
}

} // namespace LCAO_domain
