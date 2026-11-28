#ifndef ABACUS_SOURCE_LCAO_MODULE_LR_LR_GS_H
#define ABACUS_SOURCE_LCAO_MODULE_LR_LR_GS_H

#include "source_basis/module_ao/parallel_orbitals.h"
#include "source_lcao/module_lr/utils/lr_util.h"

#include <memory>
#include <ostream>
#include <string>
#include <vector>

class Charge;
class K_Vectors;
class Parallel_Grid;
class UnitCell;
struct Input_para;
namespace ModulePW
{
class PW_Basis;
}

namespace LR
{
/// Row/column distribution of the excitation amplitudes and the per-state local length.
struct AmpLayout
{
    std::vector<Parallel_Blkcyclic> para_x;
    int nloc = 0;
};

/// Band layout of the LR window for a `Parallel_Orbitals`, built on the KS process grid.
void set_band_layout(Parallel_Orbitals& pmat, int nbasis, int nbands, const UnitCell& ucell, bool aims);

/// Build the amplitude distribution; `nloc` is the local length of one state.
AmpLayout make_amp_layout(const Parallel_Blkcyclic& para_c, const std::vector<int>& nocc,
                          const std::vector<int>& nvirt, int nk, bool openshell);

/// Fill the initial guess of X: one diagonal occupied-virtual pair per state.
template <typename T>
void set_x_init(std::vector<ct::Tensor>& X, const ModuleBase::matrix& eig_ks, const std::vector<int>& nocc,
                const std::vector<int>& nvirt, const std::vector<int>& npairs,
                const std::vector<Parallel_Blkcyclic>& para_x, int nstates, int nloc, int nk, bool openshell,
                std::ostream& log);

/// AIMS benchmark: eigenvalues from band_out and eigenvectors from KS_eigenvectors.out.
template <typename T>
void read_wfc_aims(psi::Psi<T>& psi, ModuleBase::matrix& eig, const std::string& in_dir, int nocc_in, int nvirt_in,
                   int nbands, int nbasis);

/// Read the [nocc+nvirt] window; `skip_bands` drops the core bands nocc_max-nocc_in.
template <typename T>
bool read_wfc_win(const std::string& in_dir, const Parallel_Orbitals& pmat, psi::Psi<T>& psi, ModuleBase::matrix& eig,
                  ModuleBase::matrix& wg, const K_Vectors& kv, int nspin, bool binary, int skip_bands);

/// Read all ground-state bands (the force path).
template <typename T>
bool read_wfc_all(const std::string& in_dir, const Parallel_Orbitals& pmat, psi::Psi<T>& psi, ModuleBase::matrix& eig,
                  ModuleBase::matrix& wg, const K_Vectors& kv, int nspin, bool binary);

/// The Z-vector (CPSCF) window, widened to every virtual band of the ground state.
template <typename T>
struct ZWin
{
    std::unique_ptr<psi::Psi<T>> psi;
    ModuleBase::matrix eig;
    Parallel_Blkcyclic para_c;
    std::vector<Parallel_Blkcyclic> para_x;
    std::vector<int> nvirt;
    int nbands = 0;
    int nloc = 0;
};

/// (Re)build the Z window from the full ground state. `need` and `psi_all` are checked here so the
/// caller can call this unconditionally; an empty struct is returned when there is nothing to widen.
template <typename T>
ZWin<T> build_z_win(psi::Psi<T>* psi_all, const ModuleBase::matrix& eig_all, const std::vector<int>& nocc,
                    const std::vector<int>& nvirt, int nocc_max, int nbands_x, int nbasis, int nk, bool openshell,
                    const K_Vectors& kv, const Parallel_Blkcyclic& para_c, const Parallel_Orbitals& pmat,
                    const int* desc_src, bool need, bool diagnose, std::ostream& log);

/// Read the ground-state charge from cubes, or from a binary restart when the cubes are absent.
void read_gs_chg(const std::string& in_dir, const Input_para& inp, int nspin, const UnitCell& ucell,
                 Parallel_Grid& pgrid, ModulePW::PW_Basis* pw_rho, Charge& chg, int rank, std::ostream& log);

} // namespace LR

#endif // ABACUS_SOURCE_LCAO_MODULE_LR_LR_GS_H
