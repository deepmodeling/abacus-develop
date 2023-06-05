#ifndef OUTPUT_POTENTIAL_H
#define OUTPUT_POTENTIAL_H

#include <string>

#include "module_io/output_interface.h"
#include "module_basis/module_pw/pw_basis.h"
#include "module_cell/unitcell.h"
#include "module_elecstate/module_charge/charge.h"

namespace ModuleIO
{

/// @brief the output interface to write the potential
class Output_Potential : public Output_Interface
{
  public:
    Output_Potential(const ModulePW::PW_Basis_Big* pw_big,
                     ModulePW::PW_Basis* pw_rho,
                     int nspin,
                     int iter,
                     int out_pot,
                     const ModuleBase::matrix& v_effective,
                     const double* v_effective_fixed,
                     const UnitCell* ucell,
                     const Charge* const charge,
                     int precision,
                     const std::string directory,
                     const std::string tag,
                     const std::string prefix = "None");
    void write() override;

  private:
    const ModulePW::PW_Basis_Big* _pw_big;
    ModulePW::PW_Basis* _pw_rho;
    int _nspin;
    int _iter;
    int _out_pot;
    const ModuleBase::matrix& _v_effective;
    const double* _v_effective_fixed;
    const UnitCell* _ucell;
    const Charge* const _charge;
    int _precision;
    const std::string _directory;
    const std::string _prefix;
    const std::string _tag;
    std::string _fn_Pot;
    std::string _fn_ElecStatic;
};

} // namespace ModuleIO

#endif