#ifndef PW_HYBRID_NSCF_H
#define PW_HYBRID_NSCF_H

struct Input_para;
struct General_Exx_Info;
namespace ModuleESolver
{
void validate_hybrid_nscf(const Input_para& inp, const General_Exx_Info& info);
}
#endif
