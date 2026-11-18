#ifndef ABACUS_ELECTRONIC_CONTROL_H
#define ABACUS_ELECTRONIC_CONTROL_H

class Parameter;
class UnitCell;
class surchem;
struct Input_para;
namespace ModulePW { class PW_Basis; }

namespace elecstate
{
class ElecState;
struct Efermi;

// Physical requests and results; no socket or serialization dependencies.
struct ElectronicControl
{
    double nelect = 0.0;
};

struct PotentialWindow
{
    int axis = -1;
    double minimum_A = 0.0;
    double maximum_A = 0.0;
};

struct ElectronicObservables
{
    double nelect = 0.0;
    double fermi_ry = 0.0;
    double reference_ry = 0.0;
};

// Owned by the caller for one evaluation; all referenced objects outlive it.
// Each MPI rank has its own step. Only a completed result may be published.
struct ElectronicStep
{
    Parameter* parameters = nullptr;
    ElectronicControl control;
    PotentialWindow window;
    bool two_fermi = false;
    ElectronicObservables result;
    bool request_applied = false;
    bool result_ready = false;
};

void prepare_electronic_step(const UnitCell& cell, ElecState* state, ElectronicStep& step, const Input_para& input);
// Call after nuclear forces: constructing the reference updates electric-field state.
void complete_electronic_step(const UnitCell& cell, ElecState* state, ModulePW::PW_Basis* basis,
                              const surchem& solvent, const Input_para& input,
                              bool converged, ElectronicStep& step);

void update_electron_count(ElecState& state, double nelect, int nspin, double nupdown, bool two_fermi);
double electronic_chemical_potential(const Efermi& fermi, bool two_fermi);
} // namespace elecstate

#endif
