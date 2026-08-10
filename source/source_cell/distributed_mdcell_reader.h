#ifndef DISTRIBUTED_MDCELL_READER_H
#define DISTRIBUTED_MDCELL_READER_H

#include "source_cell/md_stru_file_metadata.h"

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
                            MdStruFileMetadata& stru_metadata,
                            const ModuleBase::CommunicationDomain& communication_domain);
};

#endif
