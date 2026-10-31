#ifndef DFPT_WRITE_DYN_H
#define DFPT_WRITE_DYN_H

#include <iosfwd>
#include <string>

namespace ModuleDFPT
{
class DFPT_PW_Data;
}

class UnitCell;

namespace ModuleDFPT
{
/// Write the raw (pre-ASR, pre-LO-TO) dynamical matrices in the QE dyn
/// file format, one file per q point into out_dir/dyn.q<i>. See the .cpp
/// header comment for the exact layout and conventions.
///
/// running_log (nullable): the running log stream for the created-file
/// confirmations (AGENTS rule 16); when null the confirmation goes to
/// stdout only.
void DFPT_Write_Dyn(const UnitCell& ucell,
                    const DFPT_PW_Data& data,
                    int nq,
                    const std::string& out_dir,
                    int my_rank,
                    std::ostream* running_log);
} // namespace ModuleDFPT

#endif // DFPT_WRITE_DYN_H
