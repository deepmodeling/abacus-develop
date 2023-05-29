#ifndef OUTPUT_RHO_H
#define OUTPUT_RHO_H
#include "module_basis/module_pw/pw_basis.h"
#include "module_basis/module_pw/pw_basis_big.h"
#include "module_cell/unitcell.h"
#include "module_io/output_interface.h"

namespace ModuleIO
{

class Output_Rho : public Output_Interface
{
  public:
    Output_Rho(const ModulePW::PW_Basis_Big* pw_big,
               const ModulePW::PW_Basis* pw_rho,
               int is,
               int nspin,
               const double* data,
               int iter,
               const double& ef,
               const UnitCell* ucell,
               const std::string directory,
               int precision,
               const std::string& prefix);
    void write() override;

  public:
    const ModulePW::PW_Basis_Big* _pw_big;
    const ModulePW::PW_Basis* _pw_rho;
    int _is;
    int _nspin;
    const double* _data;
    int _iter;
    const double& _ef;
    const UnitCell* _ucell;
    const std::string _directory;
    const std::string _prefix;
    std::string _fn;
    int _precision;
};

} // namespace ModuleIO

#endif