#ifndef ABACUS_RELAX_ELECTRONIC_CONTROL_H
#define ABACUS_RELAX_ELECTRONIC_CONTROL_H

class UnitCell;
class Parameter;
namespace ModuleESolver { class ESolver; }
namespace elecstate { struct ElectronicControl; }

namespace relax
{
// Apply one external electron-number request and evaluate the current geometry.
void run_electronic_step(ModuleESolver::ESolver& solver, UnitCell& cell,
                         Parameter& parameters, const elecstate::ElectronicControl& control,
                         int istep);
}
#endif
