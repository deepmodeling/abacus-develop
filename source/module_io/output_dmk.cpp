#include "module_io/output_dmk.h"

namespace ModuleIO
{

template <typename TK>
Output_DMK<TK>::Output_DMK(elecstate::DensityMatrix<TK, double>* p_DM, Parallel_Orbitals* ParaV, int nspin, int nks)
    : p_DM_(p_DM), ParaV_(ParaV), nspin_(nspin), nks_(nks)
{
}

template <typename TK>
TK* Output_DMK<TK>::get_DMK(int ik)
{
    return p_DM_->get_DMK_vector()[ik].data();
}

template class Output_DMK<double>;
template class Output_DMK<std::complex<double>>;

} // namespace ModuleIO