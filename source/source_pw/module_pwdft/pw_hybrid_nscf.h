#ifndef PW_HYBRID_NSCF_H
#define PW_HYBRID_NSCF_H

#include <string>

struct Input_para;
struct General_Exx_Info;
class UnitCell;
namespace hamilt
{
template <typename T, typename Device>
class OperatorEXXPW;
}
namespace ModulePW
{
class PW_Basis_K;
}
namespace pw
{
void validate_hybrid_nscf(const Input_para& inp, const General_Exx_Info& info);

// Read the frozen SCF ensemble in readin_dir and install it into the EXX operator,
// which takes ownership for its whole remaining lifetime.
template <typename T, typename Device>
void prepare_exx_nscf(const Input_para& inp,
                      const UnitCell& ucell,
                      const ModulePW::PW_Basis_K* pw_wfc,
                      hamilt::OperatorEXXPW<T, Device>& op_exx,
                      const std::string& readin_dir);
}
#endif
