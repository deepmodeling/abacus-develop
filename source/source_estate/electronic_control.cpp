#include "electronic_control.h"

#include "elecstate.h"
#include "electronic_occupations.h"
#include "module_pot/elpot.h"
#include "param_update.h"
#include "source_io/module_parameter/input_parameter.h"

#include <cmath>
#include <stdexcept>

namespace elecstate
{
void prepare_electronic_step(
    const UnitCell &cell,
    ElecState *state,
    ElectronicStep &step,
    const Input_para &input)
{
    if (!state || !state->charge || !state->pot)
    {
        throw std::runtime_error("electronic state is not initialized for electron control");
    }
    if (!step.parameters)
    {
        throw std::runtime_error("electronic request requires parameter storage");
    }
    if (step.request_applied || step.result_ready)
    {
        throw std::runtime_error("electronic request has already been applied");
    }

    if (step.window.axis >= 0)
    {
        if (!state->charge->rhopw)
        {
            throw std::runtime_error("electrostatic window requires the density FFT grid");
        }
        potential_window_planes(
            cell, *state->charge->rhopw, step.window.axis, step.window.minimum_A, step.window.maximum_A);
    }

    ParamUpdater::update_electron_count(step.control.nelect, *step.parameters);
    update_electron_count(*state, step.control.nelect, input.nspin, input.nupdown, step.two_fermi);
    step.request_applied = true;
}

void complete_electronic_step(
    const UnitCell &cell,
    ElecState *state,
    ModulePW::PW_Basis *basis,
    const surchem &solvent,
    const Input_para &input,
    bool converged,
    ElectronicStep &step)
{
    if (!step.request_applied || !converged || !state || !state->charge || !state->pot || !basis)
    {
        throw std::runtime_error("electronic observables require an applied request and a converged SCF state");
    }
    if (step.result_ready)
    {
        throw std::runtime_error("electronic observables have already been collected");
    }

    ElectronicObservables result;
    const int npool = input.kpar * input.bndpar;
    result.nelect = checked_occupation_electron_count(state->wg, step.control.nelect, npool, basis->poolnproc);
    result.fermi_ry = electronic_chemical_potential(state->eferm, step.two_fermi);
    if (step.window.axis >= 0)
    {
        result.reference_ry =
            electrostatic_reference(cell, *basis, *state->charge, *state->pot, solvent, input, step.window);
    }
    step.result = result;
    step.result_ready = true;
}

void update_electron_count(
    ElecState &state,
    double nelect,
    int nspin,
    double nupdown,
    bool two_fermi)
{
    state.eferm.two_efermi = two_fermi;
    state.init_nelec_spin(nspin, nelect, nupdown);
}

double electronic_chemical_potential(const Efermi &fermi, bool two_fermi)
{
    double chemical_potential = fermi.ef;
    if (two_fermi)
    {
        // Total-N derivative at fixed N_up-N_down.
        chemical_potential = 0.5 * (fermi.ef_up + fermi.ef_dw);
    }
    if (!std::isfinite(chemical_potential))
    {
        throw std::runtime_error("non-finite SCF chemical potential");
    }
    return chemical_potential;
}
} // namespace elecstate
