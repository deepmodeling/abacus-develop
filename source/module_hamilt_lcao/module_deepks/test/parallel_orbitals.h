<<<<<<< HEAD
///adapted from parallel_orbitals from module_orbital
=======
///adapted from parallel_orbitals from module_basis/module_ao
>>>>>>> 597d101b5e2f0979645e60b803172ecac0895b52
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
        void set_trace(void);

        int ncol;
        int nrow;
        int nloc;
    };
}
