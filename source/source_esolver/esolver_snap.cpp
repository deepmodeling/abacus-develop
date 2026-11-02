#include "esolver_snap.h"

#include "source_base/constants.h"
#include "source_base/global_variable.h"
#include "source_base/parallel_reduce.h"
#include "source_cell/mdcell.h"
#include "source_cell/module_neighlist/neighbor_search.h"
#include "source_io/module_output/output_log.h"
#include "source_io/module_parameter/parameter.h"

#include <array>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <stdexcept>

namespace ModuleESolver
{

namespace
{
void fail(const std::string& message)
{
    ModuleBase::WARNING_QUIT("ESolver_SNAP", message);
}

std::string full_path(const std::string& path)
{
    char* resolved = ::realpath(path.c_str(), nullptr);
    if (resolved == nullptr)
    {
        fail("Cannot resolve SNAP file path: " + path);
    }
    const std::string result(resolved);
    std::free(resolved);
    return result;
}

ModuleMD::SNAP::Neighbor make_neighbor(const NeighborAtom& center,
                                       const NeighborAtom& neighbor,
                                       const double lat0,
                                       const ModuleMD::SNAP::ModelParameters& parameters)
{
    ModuleMD::SNAP::Neighbor result;
    result.x = (neighbor.position_x - center.position_x) * lat0 * ModuleBase::BOHR_TO_A;
    result.y = (neighbor.position_y - center.position_y) * lat0 * ModuleBase::BOHR_TO_A;
    result.z = (neighbor.position_z - center.position_z) * lat0 * ModuleBase::BOHR_TO_A;
    result.weight = parameters.element_weights[static_cast<std::size_t>(neighbor.atom_type)];
    result.cutoff_radius = (parameters.element_radii[static_cast<std::size_t>(center.atom_type)]
                            + parameters.element_radii[static_cast<std::size_t>(neighbor.atom_type)])
                           * parameters.rcutfac;
    result.element_index = static_cast<std::size_t>(neighbor.atom_type);
    if (parameters.switchinner_flag)
    {
        const std::size_t center_type = static_cast<std::size_t>(center.atom_type);
        const std::size_t neighbor_type = static_cast<std::size_t>(neighbor.atom_type);
        result.sinner = 0.5 * (parameters.sinner[center_type] + parameters.sinner[neighbor_type]);
        result.dinner = 0.5 * (parameters.dinner[center_type] + parameters.dinner[neighbor_type]);
    }
    return result;
}
}

ESolver_SNAP::ESolver_SNAP(const std::string& coefficient_file,
                           const std::string& parameter_file,
                           std::ofstream& running_log)
    : coefficient_file_(coefficient_file), parameter_file_(parameter_file), running_log_(running_log)
{
    classname = "ESolver_SNAP";
}

void ESolver_SNAP::initialize(const std::vector<std::string>& labels)
{
    std::string error;
    if (!model_.load(coefficient_file_, parameter_file_, labels, error))
    {
        fail(error);
    }
    const std::string coefficient_path = full_path(coefficient_file_);
    const std::string parameter_path = full_path(parameter_file_);
    running_log_ << " Read SNAP coefficient file: " << coefficient_path << std::endl;
    running_log_ << " Read SNAP parameter file: " << parameter_path << std::endl;
    if (GlobalV::MY_RANK == 0)
    {
        std::cout << " Read SNAP coefficient file: " << coefficient_path << std::endl;
        std::cout << " Read SNAP parameter file: " << parameter_path << std::endl;
    }
    calculator_.reset(new ModuleMD::SNAP::Calculator(model_));
    snap_virial_.create(3, 3);
}

void ESolver_SNAP::before_all_runners(BaseCell& cell, const Input_para& inp)
{
    inp_ = &inp;
    if (cell.kind() == BaseCell::Kind::mdcell)
    {
        MDCell& mdcell = static_cast<MDCell&>(cell);
        initialize(mdcell.type_labels());
        mdcell.set_neighbor_cutoff(model_.parameters.rcut * ModuleBase::ANGSTROM_AU);
        return;
    }

    UnitCell& ucell = static_cast<UnitCell&>(cell);
    std::vector<std::string> labels;
    for (int it = 0; it < ucell.ntype; ++it)
    {
        labels.push_back(ucell.atoms[it].label);
    }
    initialize(labels);
    snap_force_.create(ucell.nat, 3);
}

void ESolver_SNAP::runner(BaseCell& cell, const int istep) try
{
    static_cast<void>(istep);
    if (!calculator_)
    {
        fail("SNAP solver was not initialized.");
    }
    snap_potential_ = 0.0;
    snap_virial_.zero_out();

    if (cell.kind() == BaseCell::Kind::mdcell)
    {
        MDCell& mdcell = static_cast<MDCell&>(cell);
        if (!mdcell.has_neighbor_search())
        {
            fail("MDCell neighbors must be prepared before ESolver_SNAP::runner().");
        }
        std::vector<LocalAtom>& owned_atoms = mdcell.owned_atoms();
        const std::vector<NeighborAtom>& all_atoms = mdcell.neighbor_search().get_all_atoms();
        const NeighborList& list = mdcell.neighbor_search().get_neighbor_list();
        const int nowned = static_cast<int>(owned_atoms.size());
        ModuleMD::SNAP::CellAccumulator contribution(all_atoms.size());
        for (int i = 0; i < nowned; ++i)
        {
            const NeighborAtom center(owned_atoms[static_cast<std::size_t>(i)].cart.x,
                                      owned_atoms[static_cast<std::size_t>(i)].cart.y,
                                      owned_atoms[static_cast<std::size_t>(i)].cart.z,
                                      owned_atoms[static_cast<std::size_t>(i)].type,
                                      i,
                                      i);
            std::vector<ModuleMD::SNAP::Neighbor> neighbors;
            std::vector<std::size_t> neighbor_indices;
            for (int k = 0; k < list.get_numneigh(i); ++k)
            {
                const int neighbor_index = list.get_firstneigh(i)[k];
                const NeighborAtom& neighbor = all_atoms[static_cast<std::size_t>(neighbor_index)];
                if (neighbor.atom_type < 0 || static_cast<std::size_t>(neighbor.atom_type) >= model_.parameters.elements.size())
                {
                    fail("MDCell atom type is outside the SNAP element map.");
                }
                neighbors.push_back(make_neighbor(center, neighbor, mdcell.lat0(), model_.parameters));
                neighbor_indices.push_back(static_cast<std::size_t>(neighbor_index));
            }
            const std::size_t center_element = static_cast<std::size_t>(center.atom_type);
            const ModuleMD::SNAP::ForceResult result = calculator_->compute_with_analytic_force(neighbors, center_element);
            const std::size_t center_index = static_cast<std::size_t>(i);
            contribution.add(center_index, neighbor_indices, neighbors, result);
        }
        snap_potential_ = contribution.energy() / ModuleBase::Ry_to_eV;
        const double virial_scale = 1.0 / (ModuleBase::Ry_to_eV * mdcell.omega());
        for (int row = 0; row < 3; ++row)
        {
            for (int column = 0; column < 3; ++column)
            {
                const std::size_t index = static_cast<std::size_t>(3 * row + column);
                snap_virial_(row, column) = contribution.virial()[index] * virial_scale;
            }
        }
        const double force_scale = 1.0 / (ModuleBase::Ry_to_eV * ModuleBase::ANGSTROM_AU);
        for (std::size_t i = 0; i < mdcell.ghost_atoms().size(); ++i)
        {
            const std::array<double, 3>& value = contribution.forces()[nowned + i];
            mdcell.ghost_atoms()[i].force.set(value[0] * force_scale,
                                               value[1] * force_scale,
                                               value[2] * force_scale);
        }
        Parallel_Reduce::reduce_all(snap_potential_);
        Parallel_Reduce::reduce_all(snap_virial_.c, 9);
        for (int i = 0; i < nowned; ++i)
        {
            const std::array<double, 3>& value = contribution.forces()[static_cast<std::size_t>(i)];
            owned_atoms[static_cast<std::size_t>(i)].force.set(
                value[0] * force_scale,
                value[1] * force_scale,
                value[2] * force_scale);
        }
        return;
    }

    UnitCell& ucell = static_cast<UnitCell&>(cell);
    NeighborSearch search;
    search.init(ucell, model_.parameters.rcut * ModuleBase::ANGSTROM_AU);
    search.build_neighbors();
    const NeighborList& list = search.get_neighbor_list();
    const std::vector<NeighborAtom>& inside = search.get_inside_atoms();
    const std::vector<NeighborAtom>& all = search.get_all_atoms();
    snap_force_.zero_out();
    const std::size_t atom_count = static_cast<std::size_t>(ucell.nat);
    ModuleMD::SNAP::CellAccumulator contribution(atom_count);
    std::vector<int> offsets(static_cast<std::size_t>(ucell.ntype + 1), 0);
    for (int it = 0; it < ucell.ntype; ++it)
    {
        offsets[static_cast<std::size_t>(it + 1)] = offsets[static_cast<std::size_t>(it)] + ucell.atoms[it].na;
    }
    for (int local_i = 0; local_i < list.get_ncentral_atoms(); ++local_i)
    {
        const NeighborAtom& center = inside[static_cast<std::size_t>(local_i)];
        std::vector<ModuleMD::SNAP::Neighbor> neighbors;
        std::vector<std::size_t> atom_indices;
        for (int k = 0; k < list.get_numneigh(local_i); ++k)
        {
            const NeighborAtom& neighbor = all[static_cast<std::size_t>(list.get_firstneigh(local_i)[k])];
            neighbors.push_back(make_neighbor(center, neighbor, ucell.lat0, model_.parameters));
            const int atom_index = offsets[static_cast<std::size_t>(neighbor.atom_type)] + neighbor.atom_index;
            const std::size_t index = static_cast<std::size_t>(atom_index);
            atom_indices.push_back(index);
        }
        const std::size_t center_element = static_cast<std::size_t>(center.atom_type);
        const ModuleMD::SNAP::ForceResult result = calculator_->compute_with_analytic_force(neighbors, center_element);
        const int center_index = offsets[static_cast<std::size_t>(center.atom_type)] + center.atom_index;
        const std::size_t index = static_cast<std::size_t>(center_index);
        contribution.add(index, atom_indices, neighbors, result);
    }
    snap_potential_ = contribution.energy() / ModuleBase::Ry_to_eV;
    const double virial_scale = 1.0 / (ModuleBase::Ry_to_eV * ucell.omega);
    for (int row = 0; row < 3; ++row)
    {
        for (int column = 0; column < 3; ++column)
        {
            const std::size_t index = static_cast<std::size_t>(3 * row + column);
            snap_virial_(row, column) = contribution.virial()[index] * virial_scale;
        }
    }
    const double force_scale = 1.0 / (ModuleBase::Ry_to_eV * ModuleBase::ANGSTROM_AU);
    for (int atom_index = 0; atom_index < ucell.nat; ++atom_index)
    {
        for (int direction = 0; direction < 3; ++direction)
        {
            const std::size_t index = static_cast<std::size_t>(atom_index);
            const std::size_t component = static_cast<std::size_t>(direction);
            snap_force_(atom_index, direction) = contribution.forces()[index][component] * force_scale;
        }
    }
}
catch (const std::exception& error)
{
    const std::string message = std::string("SNAP evaluation failed: ") + error.what();
    fail(message);
}

double ESolver_SNAP::cal_energy()
{
    return snap_potential_;
}

void ESolver_SNAP::cal_force(BaseCell& cell, ModuleBase::matrix& force)
{
    if (cell.kind() == BaseCell::Kind::mdcell)
    {
        const MDCell& mdcell = static_cast<const MDCell&>(cell);
        force.create(mdcell.owned_atoms().size(), 3);
        for (std::size_t i = 0; i < mdcell.owned_atoms().size(); ++i)
        {
            force(i, 0) = mdcell.owned_atoms()[i].force.x;
            force(i, 1) = mdcell.owned_atoms()[i].force.y;
            force(i, 2) = mdcell.owned_atoms()[i].force.z;
        }
        return;
    }
    force = snap_force_;
    ModuleIO::print_force(running_log_,
                          static_cast<UnitCell&>(cell),
                          "TOTAL-FORCE (eV/Angstrom)",
                          force,
                          false);
}

void ESolver_SNAP::cal_stress(BaseCell& cell, ModuleBase::matrix& stress)
{
    stress = snap_virial_;

    if (cell.kind() == BaseCell::Kind::unitcell)
    {
        ModuleIO::print_stress("TOTAL-STRESS", stress, true, false, running_log_);
    }
    else
    {
#ifdef __MPI
        if (static_cast<MDCell&>(cell).mpi_rank() == 0)
#endif
        {
            ModuleIO::print_stress("TOTAL-STRESS", stress, true, false, running_log_);
        }
    }

    const double unit_transform = ModuleBase::RYDBERG_SI / pow(ModuleBase::BOHR_RADIUS_SI, 3) * 1.0e-8;
    const double external_stress[3] = {inp_->press1, inp_->press2, inp_->press3};
    for (int i = 0; i < 3; ++i)
    {
        stress(i, i) -= external_stress[i] / unit_transform;
    }
}

void ESolver_SNAP::after_all_runners(BaseCell& cell)
{
    static_cast<void>(cell);
    running_log_ << "\n !FINAL_ETOT_IS " << snap_potential_ * ModuleBase::Ry_to_eV << " eV\n";
}

} // namespace ModuleESolver
