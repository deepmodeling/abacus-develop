#ifndef SOURCE_MD_TEST_LJ_TEST_SOLVER_H
#define SOURCE_MD_TEST_LJ_TEST_SOLVER_H

#include "source_base/global_variable.h"
#include "source_esolver/esolver_forcefields.h"
#include "source_io/module_parameter/parameter.h"
#include "source_md/module_lj/lj_forcefield.h"

#include <memory>

inline ModuleESolver::ESolver *make_lj_test_solver()
{
    std::unique_ptr<ModuleMD::ForceField> potential(
        new ModuleMD::LJ(PARAM.globalv.global_out_dir, GlobalV::ofs_running));
    return new ModuleESolver::esolver_forcefields(std::move(potential), GlobalV::ofs_running);
}

#endif
