#ifndef ABACUS_ELECTRONIC_CONTROL_H
#define ABACUS_ELECTRONIC_CONTROL_H

namespace ModuleESolver
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
}

#endif
