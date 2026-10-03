#ifndef MODULE_IO_OUTPUT_DMK_H
#define MODULE_IO_OUTPUT_DMK_H
#include "source_basis/module_ao/parallel_orbitals.h"
#include "source_estate/module_dm/density_matrix.h"

namespace ModuleIO
{

template <typename TK>
class Output_DMK
{
  public:
    Output_DMK(module_dm::DensityMatrix<TK, double>* p_DM, 
            Parallel_Orbitals* ParaV, 
            int nspin, 
            int nks);

    TK* get_dmk(int ik);

  private:
    module_dm::DensityMatrix<TK, double>* p_DM_ = nullptr;
    Parallel_Orbitals* ParaV_ = nullptr;
    int nks_;
    int nspin_;
    std::vector<TK> DMK;
};

} // namespace ModuleIO

#endif // MODULE_IO_OUTPUT_DMK_H
