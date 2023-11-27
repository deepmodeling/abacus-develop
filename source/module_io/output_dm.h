#ifndef OUTPUT_DM_H
#define OUTPUT_DM_H

#include <string>

#include "module_cell/unitcell.h"
#include "module_hamilt_lcao/module_gint/grid_technique.h"
#include "module_io/output_interface.h"

namespace ModuleIO
{
  
/// @brief the output interface to write the density matrix
class Output_DM : public Output_Interface
{
  public:
    Output_DM(const Grid_Technique& GridT,
              int is,
              int iter,
              int precision,
              int out_dm,
              double*** DM,
              const double& ef,
              const UnitCell* ucell,
              const std::string& directory,
              bool gamma_only_local);
    void write() override;

  private:
    const Grid_Technique& _GridT;
    int _is;
    int _iter;
    std::string _fn;
    int _precision;
    int _out_dm;
    double*** _DM;
    const double& _ef;
    const UnitCell* _ucell;
    const std::string& _directory;
    bool _gamma_only_local;
};
} // namespace ModuleIO

#endif // OUTPUT_DM_H