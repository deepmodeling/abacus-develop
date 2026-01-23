#ifndef POT_ML_EXX_H
#define POT_ML_EXX_H
#include "source_io/module_parameter/parameter.h"

#ifdef __MLALGO

#include "pot_base.h"
#include "source_pw/module_ofdft/ml_tools/nn_of.h"
#include "source_io/cal_mlkedf_descriptors.h"
#include "source_pw/module_ofdft/ml_base.h"

namespace elecstate
{

class ML_EXX : public ML_Base
{
public:
    ML_EXX();
    virtual ~ML_EXX();

    void set_para(const Input_para& inp, const UnitCell* ucell_in, const ModulePW::PW_Basis* rho_basis_in);

    void ml_potential(const double * const * prho, const ModulePW::PW_Basis *pw_rho, ModuleBase::matrix &rpotential);

    // output all parameters
    void generateTrainData(const double * const *prho, const ModulePW::PW_Basis *pw_rho, const double *veff);
    void localTest(const double * const *prho, const ModulePW::PW_Basis *pw_rho);

    void init_data(
        const int &nkernel,
        const bool &of_ml_gamma,
        const bool &of_ml_p,
        const bool &of_ml_q,
        const bool &of_ml_tanhp,
        const bool &of_ml_tanhq,
        const std::vector<int> &of_ml_gammanl_,
        const std::vector<int> &of_ml_pnl,
        const std::vector<int> &of_ml_qnl,
        const std::vector<int> &of_ml_xi,
        const std::vector<int> &of_ml_tanhxi,
        const std::vector<int> &of_ml_tanhxi_nl,
        const std::vector<int> &of_ml_tanh_pnl,
        const std::vector<int> &of_ml_tanh_qnl,
        const std::vector<int> &of_ml_tanhp_nl,
        const std::vector<int> &of_ml_tanhq_nl
    );

    double ml_exx_energy = 0.0;
};


class PotML_EXX : public PotBase
{
  public:
    PotML_EXX(const ModulePW::PW_Basis* rho_basis_in, const UnitCell* ucell_in)
    {
        this->rho_basis_ = rho_basis_in;
        this->dynamic_mode = true;
        this->fixed_mode = false;

        this->ml_exx.set_para(PARAM.inp, ucell_in, rho_basis_in);
    }
    ~PotML_EXX() {};

    void cal_v_eff(const Charge*const chg, const UnitCell*const ucell, ModuleBase::matrix& v_eff)
    {
        if (PARAM.inp.of_ml_local_test) this->ml_exx.localTest(chg->rho, this->rho_basis_);
        this->ml_exx.ml_potential(chg->rho, this->rho_basis_, v_eff);
    }

    double get_energy() const override { return this->ml_exx.ml_exx_energy; }

private:
   ML_EXX ml_exx;
};


}
#endif
#endif
