#ifndef ESOLVER_KS_LCAO_H
#define ESOLVER_KS_LCAO_H

#include "esolver_ks.h"

// for adjacent atoms
#include "source_lcao/hamilt_lcaodft/record_adj.h"

// for NAO basis
#include "source_basis/module_nao/two_center_bundle.h"

// for grid integration
#include "source_lcao/module_gint/gint_gamma.h"
#include "source_lcao/module_gint/gint_k.h"
#include "source_lcao/module_gint/temp_gint/gint.h"
#include "source_lcao/module_gint/temp_gint/gint_info.h"

// for DeePKS
#ifdef __MLALGO
#include "source_lcao/module_deepks/LCAO_deepks.h"
#endif

// for EXX
#ifdef __EXX
#include "module_ri/Exx_LRI_interface.h"
#include "module_ri/Mix_DMk_2D.h"
#endif

// for RDMFT
#include "module_rdmft/rdmft.h"

#include <memory>


// for Linear Response
namespace LR
{
template <typename T, typename TR>
class ESolver_LR;
}

//-----------------------------------
// ESolver for LCAO
//-----------------------------------
namespace ModuleESolver
{

template <typename TK, typename TR>
class ESolver_KS_LCAO : public ESolver_KS<TK>
{
  public:
    ESolver_KS_LCAO();
    ~ESolver_KS_LCAO();

    void before_all_runners(UnitCell& ucell, const Input_para& inp) override;

    double cal_energy() override;

    void cal_force(UnitCell& ucell, ModuleBase::matrix& force) override;

    void cal_stress(UnitCell& ucell, ModuleBase::matrix& stress) override;

    void after_all_runners(UnitCell& ucell) override;

  protected:
    virtual void before_scf(UnitCell& ucell, const int istep) override;

    virtual void iter_init(UnitCell& ucell, const int istep, const int iter) override;

    virtual void hamilt2rho_single(UnitCell& ucell, const int istep, const int iter, const double ethr) override;

    virtual void update_pot(UnitCell& ucell, const int istep, const int iter, const bool conv_esolver) override;

    virtual void iter_finish(UnitCell& ucell, const int istep, int& iter, bool& conv_esolver) override;

    virtual void after_scf(UnitCell& ucell, const int istep, const bool conv_esolver) override;

    virtual void others(UnitCell& ucell, const int istep) override;

    //! Store information about Adjacent Atoms 
    Record_adj RA;

    //! Store information about Adjacent Atoms 
    Grid_Driver gd;

    //! NAO orbitals: 2d block-cyclic distribution info
    Parallel_Orbitals pv;

    //! Grid integration: used for k-point-dependent algorithm
    Gint_k GK;

    //! Grid integration: used for gamma only algorithms.
    Gint_Gamma GG;

    //! Grid integration: used to store some basic information
    Grid_Technique GridT;

#ifndef __OLD_GINT
    //! GintInfo: used to store some basic infomation about module_gint
    std::unique_ptr<ModuleGint::GintInfo> gint_info_;
#endif

    //! NAO orbitals: two-center integrations
    TwoCenterBundle two_center_bundle_;

    //! For RDMFT calculations, added by jghan, 2024-03-16 
    rdmft::RDMFT<TK, TR> rdmft_solver;

    //! NAO: store related information 
    LCAO_Orbitals orb_;

    // Temporarily store the stress to unify the interface with PW,
    // because it's hard to seperate force and stress calculation in LCAO.
    ModuleBase::matrix scs;
    bool have_force = false;

#ifdef __MLALGO
    LCAO_Deepks<TK> ld;
#endif

#ifdef __EXX
    std::shared_ptr<Exx_LRI_Interface<TK, double>> exd = nullptr;
    std::shared_ptr<Exx_LRI_Interface<TK, std::complex<double>>> exc = nullptr;
#endif

    friend class LR::ESolver_LR<double, double>;
    friend class LR::ESolver_LR<std::complex<double>, double>;
};
} // namespace ModuleESolver
#endif
