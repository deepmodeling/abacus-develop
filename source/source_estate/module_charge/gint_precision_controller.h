#ifndef GINT_PRECISION_CONTROLLER_H
#define GINT_PRECISION_CONTROLLER_H

#include "source_lcao/module_gint/gint_precision.h"

class GintPrecisionController
{
  public:
    GintPrecisionController() = default;

    void reset_for_new_scf();

    void update_after_iteration(double drho, double scf_thr);

    ModuleGint::GintExecConfig current_config() const;

  private:
    void apply_runtime_override_();

    ModuleGint::GintExecConfig current_cfg_{
        ModuleGint::GintRealPrecision::fp32
    };
    bool locked_fp64_ = false;
    bool force_precision_enabled_ = false;
    ModuleGint::GintRealPrecision forced_precision_ = ModuleGint::GintRealPrecision::fp32;
};

#endif
