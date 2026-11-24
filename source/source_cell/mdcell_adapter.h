#ifndef SOURCE_CELL_MDCELL_ADAPTER_H
#define SOURCE_CELL_MDCELL_ADAPTER_H

#include <iosfwd>
#include <string>
#include <vector>

class DomainDecomposition;
class MDCell;
class UnitCell;

namespace ModuleBase
{
class CommunicationDomain;
class matrix;
}

namespace MDCellAdapter
{
/// Initialize the direct MDCell workflow from a STRU frame. Restart frames
/// already contain the replicated cell; the skin is supplied in Angstrom.
void initialize_from_stru(MDCell& mdcell,
                          DomainDecomposition& decomp,
                          const std::string& stru_file,
                          const std::vector<int>& cell_replica,
                          bool md_restart,
                          double neighbor_skin_angstrom,
                          std::ofstream& running_log);

/// Initialize MD for a solver that continues to use a long-lived UnitCell.
void initialize_backed_from_unitcell(MDCell& mdcell,
                                     DomainDecomposition& decomp,
                                     UnitCell& ucell);

void initialize_from_unitcell(MDCell& mdcell,
                              DomainDecomposition& decomp,
                              const UnitCell& ucell,
                              double skin,
                              const ModuleBase::CommunicationDomain& domain);

void sync_geometry_from_unitcell(MDCell& mdcell, const UnitCell& ucell);

void collect_full_forces(MDCell& mdcell,
                         DomainDecomposition& decomp,
                         const UnitCell& ucell,
                         ModuleBase::matrix& full_force);
}

#endif
