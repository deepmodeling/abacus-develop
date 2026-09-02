#ifndef RUN_MD_H
#define RUN_MD_H

#include "source_cell/md_cell.h"
#include "source_cell/md_stru_file_metadata.h"
#include "source_esolver/esolver.h"
#include "source_io/module_parameter/parameter.h"

/**
 * @brief the md loop line
 *
 */
namespace Run_MD
{
/**
 * @brief the md loop line
 *
 * @param cell cell information
 * @param p_esolver energy solver
 * @param md_para input parameters used in md
 */
void md_line(MDCell& mdcell,
             ModuleESolver::ESolver* p_esolver,
             const Parameter& param_in,
             const MdStruFileMetadata& stru_metadata);
} // namespace Run_MD

#endif
