#ifndef DISTRIBUTED_MDCELL_READER_H
#define DISTRIBUTED_MDCELL_READER_H

#include <string>

class MdCell;

class DistributedMdCellReader
{
public:
    static MdCell read_lj_stru(const std::string& stru_file,
                               double cutoff_bohr,
                               double skin_bohr);
};

#endif
