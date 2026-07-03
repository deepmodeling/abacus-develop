#include "source_lcao/module_ri/sternheimer_abacus_fd_nonlocal.h"

#include "source_cell/atom_spec.h"
#include "source_cell/unitcell.h"
#include "source_lcao/module_ri/sternheimer_fd_projector_sampler.h"

#include <stdexcept>
#include <vector>

namespace
{

void validate_grid(const ModuleRI::SternheimerFDHamiltonian::Grid& grid)
{
    if (grid.nx <= 0 || grid.ny <= 0 || grid.nz <= 0)
    {
        throw std::invalid_argument("Sternheimer ABACUS nonlocal projector requires positive grid dimensions.");
    }
    if (grid.hx <= 0.0 || grid.hy <= 0.0 || grid.hz <= 0.0)
    {
        throw std::invalid_argument("Sternheimer ABACUS nonlocal projector requires positive grid spacings.");
    }
}

void validate_unitcell(const UnitCell& ucell)
{
    if (ucell.ntype < 0)
    {
        throw std::invalid_argument("Sternheimer ABACUS nonlocal projector received negative ntype.");
    }
    if (ucell.lat0 <= 0.0)
    {
        throw std::invalid_argument("Sternheimer ABACUS nonlocal projector requires a positive lat0.");
    }
    if (ucell.ntype > 0 && ucell.atoms == nullptr)
    {
        throw std::invalid_argument("Sternheimer ABACUS nonlocal projector received a null atom list.");
    }
}

void validate_atom_pseudo(const Atom& atom)
{
    if (atom.na < 0)
    {
        throw std::invalid_argument("Sternheimer ABACUS nonlocal projector received negative atom count.");
    }
    if (static_cast<int>(atom.tau.size()) < atom.na)
    {
        throw std::invalid_argument("Sternheimer ABACUS nonlocal projector atom position count mismatch.");
    }
    if (atom.ncpp.nbeta == 0)
    {
        return;
    }
    if (atom.ncpp.tvanp)
    {
        throw std::invalid_argument(
            "Sternheimer ABACUS FD nonlocal projector currently supports norm-conserving pseudopotentials only.");
    }
    if (atom.ncpp.has_so)
    {
        throw std::invalid_argument(
            "Sternheimer ABACUS FD nonlocal projector currently does not support spin-orbit projectors.");
    }
    if (atom.ncpp.mesh <= 1 || static_cast<int>(atom.ncpp.r.size()) != atom.ncpp.mesh)
    {
        throw std::invalid_argument("Sternheimer ABACUS nonlocal projector radial mesh size mismatch.");
    }
}

int grid_size(const ModuleRI::SternheimerFDHamiltonian::Grid& grid)
{
    return grid.nx * grid.ny * grid.nz;
}

} // namespace

namespace ModuleRI
{

std::shared_ptr<SternheimerFDNonlocalProjector>
make_sternheimer_fd_nonlocal_projector_from_unitcell(const UnitCell& ucell,
                                                     const SternheimerFDHamiltonian::Grid& grid,
                                                     const double volume_element)
{
    validate_grid(grid);
    validate_unitcell(ucell);
    if (volume_element <= 0.0)
    {
        throw std::invalid_argument("Sternheimer ABACUS nonlocal projector requires a positive volume element.");
    }

    std::vector<SternheimerFDNonlocalProjector::ProjectorBlock> blocks;
    for (int it = 0; it != ucell.ntype; ++it)
    {
        const Atom& atom = ucell.atoms[it];
        validate_atom_pseudo(atom);
        if (atom.ncpp.nbeta == 0)
        {
            continue;
        }

        const SternheimerFDRadialProjectorSet radial_set
            = make_sternheimer_fd_radial_projector_set_from_abacus_matrices(atom.ncpp.r,
                                                                            atom.ncpp.betar,
                                                                            atom.ncpp.lll,
                                                                            atom.ncpp.dion);
        for (int ia = 0; ia != atom.na; ++ia)
        {
            blocks.push_back(sample_sternheimer_fd_projector_block(radial_set, grid, atom.tau[ia] * ucell.lat0));
        }
    }

    if (blocks.empty())
    {
        return nullptr;
    }
    return std::make_shared<SternheimerFDNonlocalProjector>(grid_size(grid), volume_element, blocks);
}

} // namespace ModuleRI
