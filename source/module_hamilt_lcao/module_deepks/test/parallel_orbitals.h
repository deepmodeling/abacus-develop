///adapted from parallel_orbitals from module_basis/module_ao
///deals with the parallelization of atomic basis

#include "module_base/global_function.h"
#include "module_base/global_variable.h"
#include "module_base/memory.h"

namespace Test_Deepks
{
    class Parallel_Orbitals
    {
        public:

        Parallel_Orbitals();
        ~Parallel_Orbitals();

        int* trace_loc_row;
        int* trace_loc_col;
        void set_global2local(void);

        int ncol;
        int nrow;
        int nloc;
    };
}
