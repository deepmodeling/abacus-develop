#include "lat_change_method.h"

#include "source_base/global_function.h"

Lattice_Change_Methods::Lattice_Change_Methods()
{
}
Lattice_Change_Methods::~Lattice_Change_Methods()
{
}

void Lattice_Change_Methods::allocate()
{
    Lattice_Change_Basic::dim = 9;
    lccg.allocate();

    return;
}

void Lattice_Change_Methods::cal_lattice_change(const int& istep,
                                                const int& stress_step,
                                                const ModuleBase::matrix& stress,
                                                const double& etot,
                                                UnitCell& ucell,
                                                std::ofstream& ofs,
                                                const Relax_Criteria& criteria)
{
    ModuleBase::TITLE("Lattice_Change_Methods", "lattice_change_init");
    Lattice_Change_Basic::stress_step = stress_step;

    converged_ = lccg.start(ucell, stress, etot, ofs, etot_info_, criteria);

    return;
}

void Lattice_Change_Methods::check_converged(const ModuleBase::matrix& stress,
                                             const UnitCell& ucell,
                                             std::ofstream& ofs,
                                             const Relax_Criteria& criteria)
{
    std::vector<double> lat(9, 0.0);
    std::vector<double> grad(9, 0.0);
    ModuleBase::matrix stress_copy(stress);
    Lattice_Change_Basic::setup_gradient(ucell, lat.data(), grad.data(), stress_copy);
    converged_ = Lattice_Change_Basic::check_converged(ucell, stress_copy, grad.data(), ofs, criteria.stress_thr);
}
