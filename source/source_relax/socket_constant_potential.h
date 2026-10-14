#ifndef ABACUS_SOCKET_CONSTANT_POTENTIAL_H
#define ABACUS_SOCKET_CONSTANT_POTENTIAL_H

#include "source_esolver/electronic_control.h"
#include <string>

class Parameter;
class UnitCell;
namespace ModuleESolver { class ESolver; }

namespace SocketCP
{
struct Configuration
{
    bool enabled = false;
    ModuleESolver::PotentialWindow window;
};

struct Observables
{
    double fermi_eV = 0.0;
    double nelect = 0.0;
    bool has_workfunction = false;
    double reference_eV = 0.0;
    double workfunction_eV = 0.0;
};

Configuration parse_init(const std::string& parameters);
ModuleESolver::ElectronicControl make_control(const Configuration& config, double nelect, int flags);
// Called on every rank while evaluating POSDATA; the solver may use collectives.
Observables collect_observables(const Configuration& config, const ModuleESolver::ESolver& solver,
                               const UnitCell& cell, const Parameter& parameters);
// Serializes only the completed frame, without calling the solver or MPI.
std::string serialize_observables(const Observables& result);
}

#endif
