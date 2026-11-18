#ifndef SOURCE_CELL_MDCELL_ADAPTER_H
#define SOURCE_CELL_MDCELL_ADAPTER_H

#include <string>
#include <vector>

class DomainDecomposition;
class MDCell;
class UnitCell;

namespace ModuleBase
{
class CommunicationDomain;
}

namespace MDCellAdapter
{
void initialize_from_stru(MDCell& mdcell,
                          DomainDecomposition& decomp,
                          const std::string& stru_file,
                          const std::vector<int>& cell_replica,
                          double skin,
                          const ModuleBase::CommunicationDomain& domain);

void initialize_from_unitcell(MDCell& mdcell,
                              DomainDecomposition& decomp,
                              const UnitCell& ucell,
                              double skin,
                              const ModuleBase::CommunicationDomain& domain);
}

#endif
