#ifndef RPA_LRI_INPUT_H
#define RPA_LRI_INPUT_H

class Input_para;
class K_Vectors;

namespace ModuleRI
{
void validate_librpa_v1_kpoints(const Input_para& input, const K_Vectors& kv);
}

#endif
