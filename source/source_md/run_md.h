#ifndef RUN_MD_H
#define RUN_MD_H

class MDCell;
class UnitCell;
struct Parameter;

namespace ModuleESolver
{
class ESolver;
}

/**
 * @brief the md loop line
 *
 */
namespace Run_MD
{
void prepare_mdcell(MDCell& mdcell,
                    ModuleESolver::ESolver* p_esolver,
                    const Parameter& param_in);

void prepare_mdcell(MDCell& mdcell, UnitCell& ucell);

/**
 * @brief the md loop line
 *
 * @param cell cell information
 * @param p_esolver energy solver
 * @param md_para input parameters used in md
 */
void md_line(MDCell& mdcell,
             ModuleESolver::ESolver* p_esolver,
             const Parameter& param_in);
} // namespace Run_MD

#endif
