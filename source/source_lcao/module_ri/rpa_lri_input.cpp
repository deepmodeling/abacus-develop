#include "rpa_lri_input.h"

#include "source_base/tool_quit.h"
#include "source_cell/klist.h"
#include "source_io/module_parameter/input_parameter.h"

namespace ModuleRI
{
void validate_librpa_v1_kpoints(const Input_para& input, const K_Vectors& kv)
{
    if (input.rpa && input.out_librpa_ver == 1 && !kv.is_uniform_mp())
    {
        ModuleBase::WARNING_QUIT("INPUT", "out_librpa_ver=1 requires a uniform Monkhorst-Pack k-point grid.");
    }
}
} // namespace ModuleRI
