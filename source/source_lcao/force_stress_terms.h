#ifndef FORCE_STRESS_TERMS_H
#define FORCE_STRESS_TERMS_H

#include "force_stress_lcao.h"

// Free functions for the individual LCAO force/stress correction terms
// (vdW + external fields, EXX, DFT+U, DeePKS). They are kept outside the
// Force_Stress_LCAO<T> class template; only the ones that genuinely depend on
// the electronic template type T are function templates, the vdW/fields one is
// a plain function.

class UnitCell;
class Grid_Driver;
class Parallel_Orbitals;
class LCAO_Orbitals;
class K_Vectors;
class surchem;
class pseudopot_cell_vl;
namespace ModulePW
{
class PW_Basis;
}
namespace vdw
{
struct VdwResult;
}

namespace LCAO_domain
{

// vdW force/stress and external-field forces: E-field, rt-TDDFT E-field,
// gate field and the implicit solvation model. Does not depend on T.
void cal_vdw_fields_fs(const vdw::VdwResult* vdw_result,
                       UnitCell& ucell,
                       surchem& solvent,
                       ModulePW::PW_Basis* rhopw,
                       const pseudopot_cell_vl& locpp,
                       const bool isforce,
                       const bool isstress,
                       LCAOForceParts& parts,
                       LCAOStressParts& sparts);

// EXX force/stress (only active under __EXX).
template <typename T>
void cal_exx_fs(const UnitCell& ucell,
                const bool isforce,
                const bool isstress,
                const Exx_Info& exx_info,
                Exx_NAO<T>& exx_nao,
                LCAOForceParts& parts,
                LCAOStressParts& sparts);

// DFT+U force/stress.
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
                 LCAOStressParts& sparts);

// DeePKS correction force/stress (only active under __MLALGO). The parallel
// orbitals are passed explicitly instead of going through CalEDM::ParaV.
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
                   LCAOStressParts& sparts);

} // namespace LCAO_domain

#endif
