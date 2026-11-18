#include "source_cell/mdcell_adapter.h"

#include "source_base/parallel_cell.h"
#include "source_base/matrix.h"
#include "source_cell/dom_decomp.h"
#include "source_cell/mdcell.h"
#include "source_cell/mdcell_reader.h"
#include "source_cell/print_cell.h"
#include "source_cell/unitcell.h"

#include <cstdint>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <vector>

namespace MDCellAdapter
{
void initialize_from_stru(MDCell& mdcell,
                          DomainDecomposition& decomp,
                          const std::string& stru_file,
                          const std::vector<int>& cell_replica,
                          double skin,
                          const ModuleBase::CommunicationDomain& domain)
{
    mdcell = MDCellReader::read_stru(stru_file, cell_replica, skin, domain, decomp);
}

void initialize_from_unitcell(MDCell& mdcell,
                              DomainDecomposition& decomp,
                              const UnitCell& ucell,
                              double skin,
                              const ModuleBase::CommunicationDomain& domain)
{
    decomp.init(domain, ucell.latvec, ucell.lat0, 0.0, skin);
    const std::vector<LocalAtom> owned_atoms = decomp.split_owned_atoms_from_ucell(ucell);
    std::vector<std::string> type_labels;
    std::vector<double> type_masses;
    std::vector<std::int64_t> type_atom_counts;
    for (int it = 0; it < ucell.ntype; ++it)
    {
        type_labels.push_back(ucell.atoms[it].label);
        type_masses.push_back(ucell.atoms[it].mass);
        type_atom_counts.push_back(ucell.atoms[it].na);
    }
    mdcell.initialize_from_owned_atoms(ucell.latvec,
                                       ucell.GT,
                                       ucell.lat0,
                                       ucell.omega,
                                       ucell.nat,
                                       owned_atoms,
                                       type_labels,
                                       type_masses,
                                       type_atom_counts,
                                       skin,
                                       domain);
    mdcell.mutable_stru_meta() = unitcell::make_stru_meta(ucell);
}

void sync_geometry_from_unitcell(MDCell& mdcell, const UnitCell& ucell)
{
    if (mdcell.nat() != ucell.nat || mdcell.type_atom_counts().size() != static_cast<std::size_t>(ucell.ntype))
    {
        throw std::runtime_error("UnitCell and MDCell atom counts differ.");
    }
    for (int it = 0; it < ucell.ntype; ++it)
    {
        if (mdcell.type_atom_counts()[static_cast<std::size_t>(it)] != ucell.atoms[it].na)
        {
            throw std::runtime_error("UnitCell and MDCell type atom counts differ.");
        }
    }

    const ModuleBase::Matrix3& old_lattice = mdcell.latvec();
    const ModuleBase::Matrix3& new_lattice = ucell.latvec;
    const bool changed = mdcell.lat0() != ucell.lat0
                         || old_lattice.e11 != new_lattice.e11 || old_lattice.e12 != new_lattice.e12
                         || old_lattice.e13 != new_lattice.e13 || old_lattice.e21 != new_lattice.e21
                         || old_lattice.e22 != new_lattice.e22 || old_lattice.e23 != new_lattice.e23
                         || old_lattice.e31 != new_lattice.e31 || old_lattice.e32 != new_lattice.e32
                         || old_lattice.e33 != new_lattice.e33;
    if (changed)
    {
        mdcell.set_lattice_geometry(new_lattice, ucell.lat0);
    }

    std::vector<LocalAtom>& atoms = mdcell.owned_atoms();
    for (std::size_t i = 0; i < atoms.size(); ++i)
    {
        LocalAtom& atom = atoms[i];
        if (atom.type < 0 || atom.type >= ucell.ntype || atom.type_index < 0
            || atom.type_index >= ucell.atoms[atom.type].na)
        {
            throw std::runtime_error("MDCell atom identity is outside the UnitCell.");
        }
        const int ia = static_cast<int>(atom.type_index);
        const ModuleBase::Vector3<double>& position = ucell.atoms[atom.type].tau[ia];
        ModuleBase::Vector3<double> frac = position * mdcell.GT();
        for (int axis = 0; axis < 3; ++axis)
        {
            frac[axis] -= std::floor(frac[axis]);
        }
        atom.frac = frac;
        atom.cart = frac * mdcell.latvec();
    }
}

void collect_full_forces(MDCell& mdcell,
                         DomainDecomposition& decomp,
                         const UnitCell& ucell,
                         ModuleBase::matrix& full_force)
{
    decomp.accumulate_ghost_forces(mdcell);
    if (mdcell.nat() != ucell.nat || ucell.nat < 0 || ucell.ntype < 0
        || mdcell.type_atom_counts().size() != static_cast<std::size_t>(ucell.ntype)
        || ucell.nat > std::numeric_limits<int>::max() / 3)
    {
        throw std::runtime_error("Invalid UnitCell atom count for force collection.");
    }

    std::vector<int> offsets(static_cast<std::size_t>(ucell.ntype + 1), 0);
    for (int it = 0; it < ucell.ntype; ++it)
    {
        offsets[static_cast<std::size_t>(it + 1)] = offsets[static_cast<std::size_t>(it)] + ucell.atoms[it].na;
        if (mdcell.type_atom_counts()[static_cast<std::size_t>(it)] != ucell.atoms[it].na)
        {
            throw std::runtime_error("UnitCell and MDCell type atom counts differ.");
        }
    }
    if (offsets.back() != ucell.nat)
    {
        throw std::runtime_error("Invalid UnitCell atom offsets.");
    }

    std::vector<double> values(static_cast<std::size_t>(3 * ucell.nat), 0.0);
    std::vector<int> owners(static_cast<std::size_t>(ucell.nat), 0);
    const std::vector<LocalAtom>& atoms = mdcell.owned_atoms();
    for (std::size_t i = 0; i < atoms.size(); ++i)
    {
        const LocalAtom& atom = atoms[i];
        if (atom.type < 0 || atom.type >= ucell.ntype || atom.type_index < 0
            || atom.type_index >= ucell.atoms[atom.type].na)
        {
            throw std::runtime_error("MDCell force atom identity is outside the UnitCell.");
        }
        const int iat = offsets[static_cast<std::size_t>(atom.type)] + static_cast<int>(atom.type_index);
        owners[static_cast<std::size_t>(iat)] += 1;
        values[static_cast<std::size_t>(3 * iat)] += atom.force.x;
        values[static_cast<std::size_t>(3 * iat + 1)] += atom.force.y;
        values[static_cast<std::size_t>(3 * iat + 2)] += atom.force.z;
    }

    ModuleBase::CommunicationDomain domain;
#ifdef __MPI
    domain.initialize(mdcell.communicator());
#endif
    domain.sum(owners.data(), ucell.nat);
    domain.sum(values.data(), 3 * ucell.nat);
    full_force.create(ucell.nat, 3);
    for (int iat = 0; iat < ucell.nat; ++iat)
    {
        if (owners[static_cast<std::size_t>(iat)] != 1)
        {
            throw std::runtime_error("UnitCell atom does not have exactly one MDCell owner.");
        }
        for (int axis = 0; axis < 3; ++axis)
        {
            full_force(iat, axis) = values[static_cast<std::size_t>(3 * iat + axis)];
        }
    }
}
}
