#ifndef PW_HYBRID_NSCF_H
#define PW_HYBRID_NSCF_H

#include <string>

struct Input_para;
struct General_Exx_Info;
class UnitCell;
class Exx_HelperBase;
namespace ModulePW
{
class PW_Basis_K;
}
namespace ModuleESolver
{
void validate_hybrid_nscf(const Input_para& inp, const General_Exx_Info& info);

// Read the frozen SCF ensemble in readin_dir and install it into the EXX operator
// held by exx_helper, which takes ownership for its whole remaining lifetime.
template <typename T, typename Device>
void prepare_exx_nscf(const Input_para& inp,
                      const UnitCell& ucell,
                      const ModulePW::PW_Basis_K* pw_wfc,
                      Exx_HelperBase* exx_helper,
                      const std::string& readin_dir);
}
#endif
