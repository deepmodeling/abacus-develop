#include "density_matrix.h"

#include "source_io/module_parameter/parameter.h"
#include "source_base/libm/libm.h"
#include "source_base/memory_recorder.h"
#include "source_base/timer.h"
#include "source_base/tool_title.h"
#include "source_cell/klist.h"

#include <algorithm>
#include <cstddef>
#include <memory>
#include <stdexcept>

namespace module_dm
{

// set DMK using a pointer
template <typename TK, typename TR>
void DensityMatrix<TK, TR>::set_dmk_ptr(const int ik, TK* DMK_in)
{
#ifdef __DEBUG
    assert(ik < this->_nk * this->spin_mult);
#endif
    this->dmk[ik].assign(DMK_in, DMK_in + this->pv->nrow * this->pv->ncol);
}

// set dmk element
template <typename TK, typename TR>
void DensityMatrix<TK, TR>::set_dmk(const int ispin, const int ik, const int i, const int j, const TK value)
{
#ifdef __DEBUG
    assert(ispin > 0 && ispin <= this->spin_mult);
    assert(ik >= 0 && ik < this->_nk);
#endif
    // consider transpose col=>row
    this->dmk[ik + this->_nk * (ispin - 1)][i * this->pv->nrow + j] = value;
}

// set dmk element
template <typename TK, typename TR>
void DensityMatrix<TK, TR>::set_dmk_zero()
{
    for (int ik = 0; ik < spin_mult * _nk; ik++)
    {
        std::fill(this->dmk[ik].begin(), this->dmk[ik].end(), TK{});
    }
}

template <typename TK, typename TR>
void DensityMatrix<TK, TR>::save_dmr()
{
    ModuleBase::TITLE("DensityMatrix", "save_dmr");
    ModuleBase::timer::start("DensityMatrix", "save_dmr");

    const int nnr = this->dmr[0]->get_nnr();
    // allocate if dmr_save is empty
    if (dmr_save.size() == 0)
    {
        dmr_save.resize(this->dmr.size());
    }
    // resize if dmr_save[is].size is not equal to dmr.size
    for (int is = 0; is < dmr_save.size(); is++)
    {
        if (dmr_save[is].size() != nnr)
        {
            dmr_save[is].resize(nnr);
        }
    }
    // save dmr to dmr_save
    for (int is = 0; is < this->dmr.size(); is++)
    {
        TR* DMR_pointer = this->dmr[is]->get_wrapper();
        TR* DMR_save_pointer = dmr_save[is].data();
        // copy DMR_pointer into DMR_save_pointer over [0, nnr); the
        // destination is fully overwritten, so prior zeroing is a dead
        // store, and resize above value-initializes newly added elements.
        std::copy(DMR_pointer, DMR_pointer + nnr, DMR_save_pointer);
    }

    ModuleBase::timer::end("DensityMatrix", "save_dmr");
}

// T of HContainer can be double or std::complex<double>
template class DensityMatrix<double, double>;               // Gamma-Only case
template class DensityMatrix<std::complex<double>, double>; // Multi-k case
template class DensityMatrix<std::complex<double>, std::complex<double>>; // For EXX in future

} // namespace module_dm
