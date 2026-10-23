#ifndef WRITE_EBAND_TERMS_H
#define WRITE_EBAND_TERMS_H

#include "source_base/matrix.h"
#include "source_psi/psi.h"

#include <array>
#include <complex>
#include <map>
#include <utility>
#include <vector>

class Grid_Driver;
class Parallel_Orbitals;
class UnitCell;
class Structure_Factor;
class surchem;
class Charge;
class K_Vectors;
class TwoCenterBundle;
struct Exx_Info;

namespace ModulePW
{
class PW_Basis;
}

#ifdef __EXX
namespace hamilt
{
using TAC = std::pair<int, std::array<int, 3>>;
}
namespace RI
{
template <typename T>
class Tensor;
}
#endif

namespace ModuleIO
{

/// @brief Write the band-decomposed energy of each Hamiltonian term
/// (kinetic, local/nonlocal pp, Hartree, XC) in KS orbital representation.
/// @tparam TK K-point data type (double or std::complex<double>)
/// @tparam TR Real-space data type (double or std::complex<double>)
template <typename TK, typename TR>
void write_eband_terms(const int nspin,
                       const int nbasis,
                       const int drank,
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
                       const TwoCenterBundle& two_center_bundle,
                       const Exx_Info& exx_info
#ifdef __EXX
                       ,
                       std::vector<std::map<int, std::map<hamilt::TAC, RI::Tensor<double>>>>* Hexxd,
                       std::vector<std::map<int, std::map<hamilt::TAC, RI::Tensor<std::complex<double>>>>>* Hexxc
#endif
);

} // namespace ModuleIO

#endif // WRITE_EBAND_TERMS_H
