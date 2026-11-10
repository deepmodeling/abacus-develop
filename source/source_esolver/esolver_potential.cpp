#include "esolver_potential.h"

#include "source_base/constants.h"
#include "source_cell/mdcell.h"
#include "source_io/module_output/output_log.h"
#include "source_io/module_parameter/parameter.h"
#include "source_md/potential.h"

#include <cmath>
#include <iomanip>
#include <fstream>
#include <utility>

namespace ModuleESolver
{
ESolver_Potential::ESolver_Potential(std::unique_ptr<ModuleMD::Potential> potential, std::ofstream& running_log)
    : potential_(std::move(potential)), running_log_(running_log)
{
    classname = "ESolver_Potential";
}

ESolver_Potential::~ESolver_Potential() = default;

void ESolver_Potential::before_all_runners(BaseCell& cell, const Input_para& inp)
{
    inp_ = &inp;
    stress_.create(3, 3);
    potential_->initialize(cell, inp);
}

void ESolver_Potential::runner(BaseCell& cell, int istep)
{
    static_cast<void>(istep);
    energy_ = 0.0;
    stress_.zero_out();
    potential_->compute(cell, energy_, force_, stress_);
    if (cell.kind() == BaseCell::Kind::unitcell)
    {
        running_log_ << " #TOTAL ENERGY# " << std::setprecision(11) << energy_ * ModuleBase::Ry_to_eV
                     << " eV" << std::endl;
    }
}

double ESolver_Potential::cal_energy()
{
    return energy_;
}

void ESolver_Potential::cal_force(BaseCell& cell, ModuleBase::matrix& force)
{
    if (cell.kind() == BaseCell::Kind::mdcell)
    {
        const MDCell& mdcell = static_cast<const MDCell&>(cell);
        const std::vector<LocalAtom>& atoms = mdcell.owned_atoms();
        force.create(atoms.size(), 3);
        for (std::size_t i = 0; i < atoms.size(); ++i)
        {
            force(i, 0) = atoms[i].force.x;
            force(i, 1) = atoms[i].force.y;
            force(i, 2) = atoms[i].force.z;
        }
        return;
    }
    force = force_;
    ModuleIO::print_force(running_log_, static_cast<UnitCell&>(cell), "TOTAL-FORCE (eV/Angstrom)", force, false);
}

void ESolver_Potential::cal_stress(BaseCell& cell, ModuleBase::matrix& stress)
{
    stress = stress_;
    bool print_stress = true;
#ifdef __MPI
    if (cell.kind() == BaseCell::Kind::mdcell)
    {
        print_stress = static_cast<MDCell&>(cell).mpi_rank() == 0;
    }
#endif
    if (print_stress)
    {
        ModuleIO::print_stress("TOTAL-STRESS", stress, true, false, running_log_);
    }
    const double unit_transform = ModuleBase::RYDBERG_SI / std::pow(ModuleBase::BOHR_RADIUS_SI, 3) * 1.0e-8;
    stress(0, 0) -= inp_->press1 / unit_transform;
    stress(1, 1) -= inp_->press2 / unit_transform;
    stress(2, 2) -= inp_->press3 / unit_transform;
}

void ESolver_Potential::after_all_runners(BaseCell& cell)
{
    static_cast<void>(cell);
    running_log_ << "\n --------------------------------------------" << std::endl;
    running_log_ << std::setprecision(16);
    running_log_ << " !FINAL_ETOT_IS " << energy_ * ModuleBase::Ry_to_eV << " eV" << std::endl;
    running_log_ << " --------------------------------------------\n\n" << std::endl;
}
} // namespace ModuleESolver
