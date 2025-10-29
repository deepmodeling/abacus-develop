#ifndef ESOLVER_KS_LCAO_H
#define ESOLVER_KS_LCAO_H

#include "esolver_ks.h"
#include "source_lcao/record_adj.h" // adjacent atoms
#include "source_basis/module_nao/two_center_bundle.h" // nao basis
#include "source_lcao/module_gint/gint.h" // gint
#include "source_lcao/module_gint/gint_info.h"
#include "source_lcao/setup_deepks.h" // for deepks, mohan add 20251008
#include "source_lcao/setup_exx.h" // for exx, mohan add 20251008
#include "source_lcao/module_rdmft/rdmft.h" // rdmft

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

    virtual void iter_finish(UnitCell& ucell, const int istep, int& iter, bool& conv_esolver) override;

    virtual void after_scf(UnitCell& ucell, const int istep, const bool conv_esolver) override;

    virtual void others(UnitCell& ucell, const int istep) override;

    //! Store information about Adjacent Atoms 
    Record_adj RA;

    //! Store information about Adjacent Atoms 
    Grid_Driver gd;

    //! NAO orbitals: 2d block-cyclic distribution info
    Parallel_Orbitals pv;

    //! GintInfo: used to store some basic infomation about module_gint
    std::unique_ptr<ModuleGint::GintInfo> gint_info_;

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

    // deepks method, mohan add 2025-10-08
    Setup_DeePKS<TK> deepks;

    // exact-exchange energy, mohan add 2025-10-08
    Exx_NAO<TK> exx_nao;

    friend class LR::ESolver_LR<double, double>;
    friend class LR::ESolver_LR<std::complex<double>, double>;


  public:
    const Record_adj & get_RA() const { return RA; }
    const Grid_Driver & get_gd() const { return gd; }
    const Parallel_Orbitals & get_pv() const { return pv; }
    const std::unique_ptr<ModuleGint::GintInfo> & get_gint_info() const { return gint_info_; }
    const TwoCenterBundle & get_two_center_bundle() const { return two_center_bundle_; }
    const rdmft::RDMFT<TK, TR> & get_rdmft_solver() const { return rdmft_solver; }
    const LCAO_Orbitals & get_orb() const { return orb_; }
    const ModuleBase::matrix & get_scs() const { return scs; }
    const Setup_DeePKS<TK> & get_deepks() const { return deepks; }
    const Exx_NAO<TK> & get_exx_nao() const { return exx_nao; }
};
} // namespace ModuleESolver
#endif
