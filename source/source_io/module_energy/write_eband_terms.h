#ifndef WRITE_EBAND_TERMS_H
#define WRITE_EBAND_TERMS_H

#include "source_base/matrix.h"
#include "source_psi/psi.h"

#include <iosfwd>
#include <vector>

class Grid_Driver;
class Parallel_Orbitals;
class UnitCell;
class Structure_Factor;
class surchem;
class Charge;
class K_Vectors;
class TwoCenterBundle;
class LCAO_Orbitals; // Defined in source/source_basis/module_ao/orb_read.h

namespace ModulePW
{
class PW_Basis;
}

namespace ModuleIO
{

/// @brief Write the band-decomposed energy of each Hamiltonian term
/// (kinetic, local/nonlocal pp, Hartree) in KS orbital representation.
/// The XC term is handled by write_Vxc from ctrl_scf_lcao.cpp instead.
/// @param nspin Number of spin components (1, 2 or 4)
/// @param istep Zero-based ionic step index; written one-based in the header
/// @param out_app_flag Output mode: 0 overwrites, 1 appends, 2 writes one
///        numbered eband_*g{step}.txt file per written step
/// @param etot Total energy of the current frame in Ry
/// @param my_rank World rank; only rank 0 writes the files
/// @tparam TK K-point data type (double or std::complex<double>)
/// @tparam TR Real-space data type (double or std::complex<double>)
template <typename TK, typename TR>
void write_eband_terms(const int nspin,
                       const int istep,
                       const int out_app_flag,
                       const double etot,
                       const int my_rank,
                       const int nbasis,
                       const Parallel_Orbitals* pv,
                       const psi::Psi<TK>& psi,
                       const UnitCell& ucell,
                       Structure_Factor& sf,
                       surchem& solvent,
                       const ModulePW::PW_Basis& rho_basis,
                       const ModulePW::PW_Basis& rhod_basis,
                       const ModuleBase::matrix& vloc,
                       const Charge& chg,
                       const K_Vectors& kv,
                       const ModuleBase::matrix& wg,
                       Grid_Driver& gd,
                       const std::vector<double>& orb_cutoff,
                       const LCAO_Orbitals& orb,
                       const TwoCenterBundle& two_center_bundle,
                       std::ofstream& ofs_running
);

} // namespace ModuleIO

#endif // WRITE_EBAND_TERMS_H
