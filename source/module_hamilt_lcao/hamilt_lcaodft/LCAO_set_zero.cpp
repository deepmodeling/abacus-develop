#include "module_hamilt_lcao/hamilt_lcaodft/LCAO_domain.h"
#include "module_base/global_variable.h"
#include "module_base/timer.h"
#include "module_base/tool_threading.h"

namespace LCAO_domain
{

void zeros_HSR(const char &mtype, LCAO_HS_Arrays& HS_arrays)
{
    auto zeros_HSR_ker = [&](int num_threads, int thread_id)
    {
        long long beg, len;
        if(GlobalV::NSPIN!=4)
        {
            if (mtype=='T')
            {
                ModuleBase::BLOCK_TASK_DIST_1D(num_threads, thread_id, (long long)HS_arrays.Hloc_fixedR.size(), (long long)512, beg, len);
                ModuleBase::GlobalFunc::ZEROS(HS_arrays.Hloc_fixedR.data() + beg, len);
            }
        }
    };
    ModuleBase::OMP_PARALLEL(zeros_HSR_ker);
    return;
}

}
