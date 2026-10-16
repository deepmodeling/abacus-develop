#include "density_matrix.h"

#include "source_io/module_parameter/parameter.h"
#include "source_base/libm/libm.h"
#include "source_base/memory_recorder.h"
#include "source_base/timer.h"
#include "source_base/tool_title.h"
#include "source_base/tool_quit.h"
#include "source_base/constants.h"
#include "source_cell/klist.h"

namespace module_dm
{

//----------------------------------------------------
// density matrix class
//----------------------------------------------------

// destructor
template <typename TK, typename TR>
DensityMatrix<TK, TR>::~DensityMatrix()
{
    this->clear_dmr();
}

template <typename TK, typename TR>
void DensityMatrix<TK, TR>::clear_dmr()
{
    for (hamilt::HContainer<TR>*& it: this->dmr)
    {
        delete it;
    }
    this->dmr.clear();
    this->_dmr_ready = false;
}

template <typename TK, typename TR>
DensityMatrix<TK, TR>::DensityMatrix(const Parallel_Orbitals* paraV_in,
                                     const int spin_mult,
                                     const std::vector<ModuleBase::Vector3<double>>& kvec_d,
                                     const int nk,
                                     const int nspin)
    : pv(paraV_in), spin_mult(spin_mult), nspin(nspin > 0 ? nspin : spin_mult),
      _kvec_d(kvec_d), _nk((nk > 0 && nk <= _kvec_d.size()) ? nk : _kvec_d.size())
{
    ModuleBase::TITLE("DensityMatrix", "resize_DMK");
    const int nks = _nk * this->spin_mult;
    this->dmk.resize(nks);
    for (int ik = 0; ik < nks; ik++)
    {
        this->dmk[ik].resize(this->pv->get_row_size() * this->pv->get_col_size());
    }
    ModuleBase::Memory::record("DensityMatrix::DMK", this->dmk.size() * this->dmk[0].size() * sizeof(TK));
}

template <typename TK, typename TR>
DensityMatrix<TK, TR>::DensityMatrix(const Parallel_Orbitals* paraV_in, const int spin_mult, const int nspin)
    : pv(paraV_in), spin_mult(spin_mult), nspin(nspin > 0 ? nspin : spin_mult),
      _kvec_d({ModuleBase::Vector3<double>(0, 0, 0)}), _nk(1)
{
    ModuleBase::TITLE("DensityMatrix", "resize_gamma");
    this->dmk.resize(this->spin_mult);
    for (int ik = 0; ik < this->spin_mult; ik++)
    {
        this->dmk[ik].resize(this->pv->get_row_size() * this->pv->get_col_size());
    }
    ModuleBase::Memory::record("DensityMatrix::DMK", this->dmk.size() * this->dmk[0].size() * sizeof(TK));
}


// switch_dmr
template <typename TK, typename TR>
void DensityMatrix<TK, TR>::switch_dmr(const int mode)
{
    ModuleBase::TITLE("DensityMatrix", "switch_dmr");
    ModuleBase::timer::start("DensityMatrix", "switch_dmr");
    if (this->spin_mult == 2)
    {
        switch(mode)
        {
        case 0:
            // switch to original density matrix
            if (!this->dmr_tmp.empty() && this->dmr_origin.size() != 0)
            {
                this->dmr[0]->allocate(this->dmr_origin.data(), false);
                this->dmr_tmp.clear();
            }
            // else: do nothing
            break;
        case 1:
            // switch to total magnetization density matrix, dmr_up + dmr_down
            if(this->dmr_tmp.empty())
            {
                const size_t size = this->dmr[0]->get_nnr();
                this->dmr_tmp.resize(size);
                this->dmr_origin.resize(size);
                for (int i = 0; i < size; ++i)
                {
                    this->dmr_origin[i] = this->dmr[0]->get_wrapper()[i];
                    this->dmr_tmp[i] = this->dmr_origin[i] + this->dmr[1]->get_wrapper()[i];
                }
                this->dmr[0]->allocate(this->dmr_tmp.data(), false);
            }
            else
            {
                const size_t size = this->dmr[0]->get_nnr();
                for (int i = 0; i < size; ++i)
                {
                    this->dmr_tmp[i] = this->dmr_origin[i] + this->dmr[1]->get_wrapper()[i];
                }
            }
            break;
        case 2:
            // switch to magnetization density matrix, dmr_up - dmr_down
            if(this->dmr_tmp.empty())
            {
                const size_t size = this->dmr[0]->get_nnr();
                this->dmr_tmp.resize(size);
                this->dmr_origin.resize(size);
                for (int i = 0; i < size; ++i)
                {
                    this->dmr_origin[i] = this->dmr[0]->get_wrapper()[i];
                    this->dmr_tmp[i] = this->dmr_origin[i] - this->dmr[1]->get_wrapper()[i];
                }
                this->dmr[0]->allocate(this->dmr_tmp.data(), false);
            }
            else
            {
                const size_t size = this->dmr[0]->get_nnr();
                for (int i = 0; i < size; ++i)
                {
                    this->dmr_tmp[i] = this->dmr_origin[i] - this->dmr[1]->get_wrapper()[i];
                }
            }
            break;
        default:
            ModuleBase::WARNING_QUIT("density_matrix.cpp", "Unknown mode in switch_dmr");
        }
    }
    ModuleBase::timer::end("DensityMatrix", "switch_dmr");
}



// T of HContainer can be double or complex<double>
template class DensityMatrix<double, double>;               // Gamma-Only case
template class DensityMatrix<std::complex<double>, double>; // Multi-k case
template class DensityMatrix<std::complex<double>, std::complex<double>>; // For EXX in future

} // namespace module_dm
