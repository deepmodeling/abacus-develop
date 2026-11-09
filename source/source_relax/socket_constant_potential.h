#ifndef ABACUS_SOCKET_CONSTANT_POTENTIAL_H
#define ABACUS_SOCKET_CONSTANT_POTENTIAL_H

#include "source_estate/electronic_control.h"
#include <string>

class UnitCell;
struct Input_para;
struct System_para;
namespace ModuleESolver { class ESolver; }

namespace SocketCP
{
struct Configuration
{
    bool enabled = false;
    elecstate::PotentialWindow window;
};

struct Observables
{
    double fermi_eV = 0.0;
    double nelect = 0.0;
    bool has_workfunction = false;
    double reference_eV = 0.0;
    double workfunction_eV = 0.0;
};

// Borrow one rank-local request through the existing cell argument for this evaluation.
class ScopedElectronicStep
{
  public:
    ScopedElectronicStep(UnitCell& cell, elecstate::ElectronicStep* step);
    ~ScopedElectronicStep();
    ScopedElectronicStep(const ScopedElectronicStep&) = delete;
    ScopedElectronicStep& operator=(const ScopedElectronicStep&) = delete;

  private:
    UnitCell& cell_;
};

// Apply one external electron-number request and evaluate the current geometry.
void run_electronic_step(ModuleESolver::ESolver& solver, UnitCell& cell,
                         const Input_para& input, const System_para& system,
                         int functional_type, int istep);

Configuration parse_init(const std::string& parameters);
elecstate::ElectronicControl make_control(const Configuration& config, double nelect, int flags);
// Converts only the physical result completed after this frame's force calculation.
Observables collect_observables(const Configuration& config, const elecstate::ElectronicStep& step);
// Serializes only the completed frame, without calling the solver or MPI.
std::string serialize_observables(const Observables& result);
}

#endif
