#ifndef PW_HYBRID_NSCF_H
#define PW_HYBRID_NSCF_H

#include <string>

struct Input_para;
struct General_Exx_Info;
class K_Vectors;
namespace ModuleBase
{
class matrix;
}
namespace ModulePW
{
class PW_Basis_K;
}

namespace ModuleESolver
{
void validate_hybrid_nscf(const Input_para& inp, const General_Exx_Info& info);
void invalidate_exx_source(const Input_para& inp,
                           const General_Exx_Info& info,
                           const ModulePW::PW_Basis_K& basis,
                           const std::string& out_dir);
void save_exx_source(const Input_para& inp,
                     const General_Exx_Info& info,
                     const K_Vectors& points,
                     const ModuleBase::matrix& weights,
                     const ModulePW::PW_Basis_K& basis,
                     const std::string& out_dir);
}
#endif
