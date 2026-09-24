#ifndef CHARGE_H
#define CHARGE_H

#include <vector>

#include "source_base/complexmatrix.h"
#include "source_base/global_function.h"
#include "source_base/global_variable.h"
#include "source_basis/module_pw/pw_basis.h"
#include "source_cell/module_symmetry/symmetry.h"
// #include "source_estate/fp_energy.h"
#include "source_base/parallel_grid.h"

//a forward declaration of UnitCell
class UnitCell;

namespace module_charge
{
struct InitRhoCfg;
}

// Electron Charge Density
class Charge
{

  public:

    Charge();
    ~Charge();

    // rho/rhog/kin_r views alias the vector-backed _space_* storage, so
    // copying a Charge would duplicate dangling pointers into another
    // object's vector buffer. Forbid copies until a deep copy is needed.
    Charge(const Charge&) = delete;
    Charge& operator=(const Charge&) = delete;

    //==========================================================
    // MEMBER VARIABLES :
    // init_chg : "atomic" or "file"
    // NAME : total number of electrons
    // NAME : rho (nspin,ncxyz), the charge density in real space
    // NAME : rho_save (nspin,ncxyz), for charge mixing
    // NAME : rhog, charge density in G space
    // NAME : rhog_save, chage density in G space
    // NAME : rho_core [nrxx], the core charge in real space
    // NAME : rhog_core [ngm], the core charge in reciprocal space
    //==========================================================

    double **rho = nullptr;
    double **rho_save = nullptr;

    std::complex<double> **rhog = nullptr;
    std::complex<double> **rhog_save = nullptr;

    double **kin_r = nullptr; // kinetic energy density in real space, for meta-GGA
    double **kin_r_save = nullptr; // same as kin_r, kept for mixing
    const Parallel_Grid* pgrid = nullptr;

  private:

    // Underlying contiguous storage backing the public rho/rhog/kin_r views.
    // Each buffer holds nspin rows; rho[is] points at _space_rho.data()+is*nrxx.
    // Owned here as std::vector so the storage self-manages (no raw new/delete).
    std::vector<double> _space_rho;
    std::vector<double> _space_rho_save;
    std::vector<std::complex<double>> _space_rhog;
    std::vector<std::complex<double>> _space_rhog_save;
    std::vector<double> _space_kin_r;
    std::vector<double> _space_kin_r_save;

    // Pointer arrays backing the public double** views (rho, rhog, etc.)
    std::vector<double*> _ptrs_rho;
    std::vector<std::complex<double>*> _ptrs_rhog;
    std::vector<double*> _ptrs_rho_save;
    std::vector<std::complex<double>*> _ptrs_rhog_save;
    std::vector<double*> _ptrs_kin_r;
    std::vector<double*> _ptrs_kin_r_save;

    // Contiguous storage for rho_core and rhog_core
    std::vector<double> _space_rho_core;
    std::vector<std::complex<double>> _space_rhog_core;

  public:

    double *rho_core = nullptr;
    std::complex<double> *rhog_core = nullptr;

    void set_rhopw(ModulePW::PW_Basis* rhopw_in);

    /**
     * @brief Init charge density from file or atomic pseudo-wave-functions
     *
     * @param ucell [in] unit cell
     * @param pgrid [in] parallel grid descriptor
     * @param strucFac [in] structure factor
     * @param symm [in] symmetry
     * @param klist [in] k points list if needed
     * @param wfcpw [in] PW basis for wave function if needed
     * @param cfg [in] INPUT values for charge initialization
     */
    void init_rho(const UnitCell& ucell,
                  const Parallel_Grid& pgrid,
                  const ModuleBase::ComplexMatrix& strucFac,
                  ModuleSymmetry::Symmetry& symm,
                  const void* klist,
                  const void* wfcpw,
                  const module_charge::InitRhoCfg& cfg);

    /**
     * @brief Allocate the rho/rhog/kin_r buffers
     *
     * @param nspin_in number of spins
     * @param kin_den whether to allocate the kinetic-energy density buffers
     * @param meta_gga whether the functional is meta-GGA (kin_r carries XC
     *        physics, not just ELF output); stored for tau handling
     * @param test_charge verbosity flag (PARAM.inp.test_charge)
     */
    void allocate(const int &nspin_in, const bool kin_den, const bool meta_gga,
                  const int test_charge);

    /**
     * @brief Renormalize rho so that its integral equals the electron number
     *
     * @param nelec target total electron number (PARAM.inp.nelec)
     * @param omega current unit-cell volume. Must be ucell.omega, NOT
     *        rhopw->omega, because in variable-cell calculations (e.g. NPT)
     *        rhopw->omega is stale (pw_rho/pw_rhod are not rebuilt on cell
     *        change) while ucell.omega is updated every step. Using the stale
     *        volume gives a wrong electron count and a wrong renormalization
     *        factor, which corrupts the stress.
     */
    void renormalize_rho(const double nelec, const double omega);

    double sum_rho(const double omega) const;

    void save_rho_before_sum_band();

    /// @brief whether allocate() has already claimed the density buffers
    /// @return true once allocate() has run; allocate() is idempotent and
    ///         returns early when this is already true
    bool get_allocate_rho() const { return allocate_rho; }

    // mohan add 2021-02-20
    int nrxx=0; // number of r vectors in this processor
    int nxyz = 0; // total number of r vectors
    int ngmc=0; // number of g vectors in this processor
    int nspin=0; // number of spins
    bool meta_gga = false; // whether the functional is meta-GGA (set by allocate)
    ModulePW::PW_Basis* rhopw = nullptr;// When double_grid is used, rhopw = rhodpw (dense grid)

  private:

    void destroy();    // free arrays  liuyu 2023-03-12

    bool allocate_rho;
};

#endif // charge
