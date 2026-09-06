#ifndef MDCELL_READER_H
#define MDCELL_READER_H

#include <string>
#include <vector>

class MDCell;
namespace ModuleBase
{
class CommunicationDomain;
}

class MDCellReader
{
public:
    static MDCell read_stru(const std::string& stru_file,
                            const std::vector<int>& cell_replica,
                            double cutoff,
                            double skin,
                            const ModuleBase::CommunicationDomain& comm_domain);
};

#endif
