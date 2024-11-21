#ifndef ESOLVER_KS_LCAO_H
#define ESOLVER_KS_LCAO_H
#include "esolver_ks.h"
#include "module_hamilt_lcao/hamilt_lcaodft/record_adj.h"
// for grid integration
#include "module_hamilt_lcao/module_gint/gint_gamma.h"
#include "module_hamilt_lcao/module_gint/gint_k.h"
#ifdef __EXX
#include "module_ri/Exx_LRI_interface.h"
#include "module_ri/Mix_DMk_2D.h"
#endif
#include "module_basis/module_nao/two_center_bundle.h"
#include "module_io/output_mat_sparse.h"

#include <memory>

namespace LR
{
    template<typename T, typename TR>
    class ESolver_LR;
}
namespace ModuleESolver
{
template <typename TK, typename TR>
class ESolver_KS_LCAO : public ESolver_KS<TK> {
  public:
    ESolver_KS_LCAO();
    ~ESolver_KS_LCAO();

    void before_all_runners(const Input_para& inp, UnitCell& cell) override;

    double cal_energy() override;

    void cal_force(ModuleBase::matrix& force) override;

    void cal_stress(ModuleBase::matrix& stress) override;

    void after_all_runners() override;

  protected:
    virtual void before_scf(const int istep) override;

    virtual void iter_init(const int istep, const int iter) override;

    virtual void hamilt2density_single(const int istep, const int iter, const double ethr) override;

    virtual void update_pot(const int istep, const int iter) override;

    virtual void iter_finish(const int istep, int& iter) override;

    virtual void after_scf(const int istep) override;

    virtual void others(const int istep) override;

    // we will get rid of this class soon, don't use it, mohan 2024-03-28
    Record_adj RA;

    // 2d block-cyclic distribution info
    Parallel_Orbitals pv;

    // used for k-dependent grid integration.
    Gint_k GK;

    // used for gamma only algorithms.
    Gint_Gamma GG;

    Grid_Technique GridT;

    TwoCenterBundle two_center_bundle_;

    // temporary introduced during removing GlobalC::ORB
    LCAO_Orbitals orb_;

    // Temporarily store the stress to unify the interface with PW,
    // because it's hard to seperate force and stress calculation in LCAO.
    // The copy costs memory and time !
    // Are there any better way to solve this problem ?
    ModuleBase::matrix scs;
    bool have_force = false;

    //--------------common for all calculation, not only scf-------------
    // set matrix and grid integral
    void set_matrix_grid(Record_adj& ra);

    void beforesolver(const int istep);
    //---------------------------------------------------------------------

#ifdef __EXX
    std::shared_ptr<Exx_LRI_Interface<TK, double>> exd = nullptr;
    std::shared_ptr<Exx_LRI_Interface<TK, std::complex<double>>> exc = nullptr;
    std::shared_ptr<Exx_LRI<double>> exx_lri_double = nullptr;
    std::shared_ptr<Exx_LRI<std::complex<double>>> exx_lri_complex = nullptr;
#endif

  private:
#ifdef __DEEPKS
    void dpks_cal_e_delta_band(const std::vector<std::vector<TK>>& dm) const;

    void dpks_cal_projected_DM(
        const elecstate::DensityMatrix<TK, double>* dm) const;
#endif
    friend class LR::ESolver_LR<double, double>;
    friend class LR::ESolver_LR<std::complex<double>, double>;
};
} // namespace ModuleESolver
#endif
