#ifndef GINT_PRECISION_CONTROLLER_H
#define GINT_PRECISION_CONTROLLER_H

#include "source_lcao/module_gint/gint_precision.h"

class GintPrecisionController
{
  public:
    GintPrecisionController() = default;

    void reset_for_new_scf();

    void update_after_iteration(int iter,
                                double drho,
                                double scf_thr,
                                bool conv_esolver,
                                bool is_restart_step);

    ModuleGint::GintExecConfig current_config() const;

  private:
    void apply_runtime_override_();

    ModuleGint::GintExecConfig current_cfg_{
        ModuleGint::GintRealPrecision::fp32
    };
    int qualified_small_drho_iters_ = 0;
    bool locked_fp64_ = false;
    bool force_precision_enabled_ = false;
    ModuleGint::GintRealPrecision forced_precision_ = ModuleGint::GintRealPrecision::fp32;
};

#endif
