#ifndef WRITE_HS_R_TERMS_H
#define WRITE_HS_R_TERMS_H

#include "source_basis/module_nao/two_center_bundle.h"
#include "source_basis/module_pw/pw_basis.h"
#include "source_cell/klist.h"
#include "source_cell/module_neighbor/sltk_grid_driver.h"
#include "source_estate/module_charge/charge.h"
#include "source_estate/module_pot/potential_new.h"
#include "source_lcao/LCAO_domain.h"
#include "source_lcao/module_hcontainer/hcontainer.h"

#include <vector>

namespace ModuleIO
{

void write_h_t(const UnitCell& ucell,
               const Grid_Driver& gd,
               const Parallel_Orbitals& pv,
               const TwoCenterBundle& two_center_bundle,
               const LCAO_Orbitals& orb,
               const K_Vectors& kv,
               const int nspin,
               const int istep,
               const bool append,
               const int* iat2iwt,
               const int nat,
               const bool also_hk = true);

void write_h_vnl(const UnitCell& ucell,
                 const Grid_Driver& gd,
                 const Parallel_Orbitals& pv,
                 const TwoCenterBundle& two_center_bundle,
                 const LCAO_Orbitals& orb,
                 const K_Vectors& kv,
                 const int nspin,
                 const int istep,
                 const bool append,
                 const int* iat2iwt,
                 const int nat,
                 const bool also_hk = true);

void write_h_vl(const UnitCell& ucell,
                const Grid_Driver& gd,
                const Parallel_Orbitals& pv,
                const LCAO_Orbitals& orb,
                const elecstate::Potential* pot,
                const int nspin,
                const int istep,
                const bool append,
                const int* iat2iwt,
                const int nat,
    const K_Vectors& kv,
                const bool also_hk = true);

void write_h_vh(const UnitCell& ucell,
                const Grid_Driver& gd,
                const Parallel_Orbitals& pv,
                const LCAO_Orbitals& orb,
                const Charge* chg,
                const ModulePW::PW_Basis* rho_basis,
                const int nspin,
                const int istep,
                const bool append,
                const int* iat2iwt,
                const int nat,
    const K_Vectors& kv,
                const bool also_hk = true);

void write_h_vxc(const UnitCell& ucell,
                 const Grid_Driver& gd,
                 const Parallel_Orbitals& pv,
                 const LCAO_Orbitals& orb,
                 const Charge* chg,
                 const int nrxx,
                 const int nspin,
                 const int istep,
                 const bool append,
                 const int* iat2iwt,
                 const int nat,
    const K_Vectors& kv,
                 const bool also_hk = true);

} // namespace ModuleIO

#endif
