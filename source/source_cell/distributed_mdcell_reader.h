#ifndef DISTRIBUTED_MDCELL_READER_H
#define DISTRIBUTED_MDCELL_READER_H

#include <string>
#include <vector>

class MDCell;
namespace ModuleBase
{
class CommunicationDomain;
}

class DistributedMDCellReader
{
public:
    static MDCell read_stru(const std::string& stru_file,
                            const std::vector<int>& replicate,
                            double cutoff,
                            double skin,
                            const ModuleBase::CommunicationDomain& communication_domain);
};

#endif
