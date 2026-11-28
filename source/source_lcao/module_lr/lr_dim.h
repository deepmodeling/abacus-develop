#ifndef ABACUS_SOURCE_LCAO_MODULE_LR_LR_DIM_H
#define ABACUS_SOURCE_LCAO_MODULE_LR_LR_DIM_H

#include <ostream>
#include <string>
#include <vector>

namespace ModuleBase
{
class matrix;
}
class UnitCell;
struct Input_para;

namespace LR
{
/// Resolved dimensions of the LR(-TDDFT) problem.
struct Dim
{
    int nspin = 1;
    int nstates = 1;
    int nbasis = 1;
    int nbands = 2;
    int nk = 1;
    int nupdown = 0;
    int nocc_in = 1;
    int nocc_max = 1;
    int nvirt_in = 1;
    int ks_nbands = 0;
    bool openshell = false;
    std::vector<int> nocc;
    std::vector<int> nvirt;
    std::vector<int> npairs;
};

/// Everything `resolve_dim` needs from the caller; no solver state is read implicitly.
struct DimSrc
{
    const Input_para* inp = nullptr;
    int nlocal = 0;                             ///< PARAM.globalv.nlocal
    int nks = 0;                                ///< K_Vectors::get_nks()
    int nkstot = 0;                             ///< K_Vectors::get_nkstot()
    int rank = 0;                               ///< MPI rank in the world communicator
    bool gamma_only = false;                    ///< real build
    const ModuleBase::matrix* wg_ks = nullptr;  ///< KS occupations; null on the file (`lr`) path
    int kpar = 1;                               ///< K_Vectors::para_k.kpar
    int rank_in_pool = 0;                       ///< K_Vectors::para_k.rank_in_pool
    std::string in_dir;
    UnitCell* ucell = nullptr;     ///< only written on the AIMS benchmark path
    std::ostream* log = nullptr;   ///< `ofs_running`; may be null
};

/// Resolve nspin/nstates/nbasis/nbands/nk/nocc/nvirt/npairs from input and the ground state.
Dim resolve_dim(const DimSrc& src);

/// AIMS / AIMS-LibRPA benchmark: basis size and band count come from `aims_nbasis` + band_out.
void resolve_aims_basis(UnitCell& ucell, const Input_para& inp, Dim& dim);

/// `nspin == 2`: split nocc/nvirt between the spin channels and decide open/closed shell.
void apply_spin_pol(Dim& dim, const Input_para& inp);

} // namespace LR

#endif // ABACUS_SOURCE_LCAO_MODULE_LR_LR_DIM_H
