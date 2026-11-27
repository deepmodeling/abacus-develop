#ifndef CHG_INIT_H
#define CHG_INIT_H

#include "source_base/complexmatrix.h"
#include "source_base/parallel_grid.h"
#include "source_basis/module_pw/pw_basis.h"
#include "source_cell/module_symmetry/symmetry.h"

#include <string>

class Charge;
class UnitCell;

namespace module_charge
{

/// Configuration for charge-density initialization, replacing direct
/// PARAM reads in init_rho. Callers fill this from the parsed input
/// once per run.
struct InitRhoCfg
{
    std::string init_chg;          ///< charge initialization mode (PARAM.inp.init_chg)
    std::string suffix;            ///< system suffix for file names (PARAM.inp.suffix)
    std::string esolver_type;      ///< esolver type; "sdft" selects the band-group rank (PARAM.inp.esolver_type)
    std::string global_readin_dir; ///< directory to read files from (PARAM.globalv.global_readin_dir)
    double nelec = 0.0;            ///< target total electron number (PARAM.inp.nelec)
    int nbands = 0;                ///< number of bands for wfc-based init (PARAM.inp.nbands)
    int test_charge = 0;           ///< verbosity flag (PARAM.inp.test_charge)
    bool domag = false;            ///< whether to compute magnetization (PARAM.globalv.domag)
    bool domag_z = false;          ///< whether to compute z-only magnetization (PARAM.globalv.domag_z)
    bool meta_gga = false;         ///< whether the functional is meta-GGA (XC_Functional::get_ked_flag())
    int npol = 1;                  ///< number of polarization components (PARAM.globalv.npol)
};

/**
 * @brief Initialize charge density from file, atomic superposition,
 *        restart cache, or wavefunctions, depending on cfg.init_chg.
 *
 * @param chr [inout] Charge object whose rho/rhog/kin_r buffers are filled.
 * @param rhopw [in] plane-wave basis bound to chr (grid sizes and FFT backend)
 * @param ucell [in] unit cell
 * @param pgrid [in] parallel grid descriptor
 * @param strucFac [in] structure factor for atomic-charge superposition
 * @param symm [in] symmetry operations (used by wfc-based init)
 * @param klist [in] k-point list pointer (K_Vectors*), needed only for wfc init
 * @param wfcpw [in] PW_Basis_K pointer, needed only for wfc init
 * @param cfg [in] INPUT values for charge initialization
 */
void init_rho(Charge& chr,
              const ModulePW::PW_Basis& rhopw,
              const UnitCell& ucell,
              const Parallel_Grid& pgrid,
              const ModuleBase::ComplexMatrix& strucFac,
              ModuleSymmetry::Symmetry& symm,
              const void* klist,
              const void* wfcpw,
              const InitRhoCfg& cfg);

} // namespace module_charge

#endif // CHG_INIT_H
