/// adapted from parallel_orbitals from source_basis/module_ao
/// deals with the parallelization of atomic basis

#include "source_base/global_function.h"
#include "source_base/global_variable.h"
#include "source_base/memory.h"

namespace Test_Deepks
{
class Parallel_Orbitals
{
  public:
    Parallel_Orbitals();
    ~Parallel_Orbitals();

    int* global2local_row;
    int* global2local_col;
    void set_global2local(void);

    int ncol;
    int nrow;
    int nloc;
};
} // namespace Test_Deepks
