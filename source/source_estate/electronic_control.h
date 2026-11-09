#ifndef ABACUS_ELECTRONIC_CONTROL_H
#define ABACUS_ELECTRONIC_CONTROL_H

class Parameter;

namespace elecstate
{
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
}

#endif
