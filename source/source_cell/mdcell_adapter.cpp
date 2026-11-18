#include "source_cell/mdcell_adapter.h"

#include "source_base/parallel_cell.h"
#include "source_cell/dom_decomp.h"
#include "source_cell/mdcell.h"
#include "source_cell/mdcell_reader.h"
#include "source_cell/print_cell.h"
#include "source_cell/unitcell.h"

#include <cstdint>
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
}
