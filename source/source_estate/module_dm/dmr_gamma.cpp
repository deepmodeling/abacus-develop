#include "density_matrix.h"

#include "source_base/libm/libm.h"
#include "source_base/timer.h"
#include "source_base/tool_title.h"

namespace module_dm
{

// calculate DMR from DMK using blas for gamma-only calculation
template <>
void DensityMatrix<double, double>::cal_dmr(const int ik_in)
{
    ModuleBase::TITLE("DensityMatrix", "cal_dmr");
    using TK = double;
    using TR = double;

    assert(ik_in == -1 || ik_in == 0);
    assert(this->_nk == 1);

    assert(this->dmr.size()==this->spin_mult && "DMR has not been initialized!");

    ModuleBase::timer::start("DensityMatrix", "cal_dmr");
    const int ld_hk = this->pv->nrow;
    for (int is = 1; is <= this->spin_mult; ++is)
    {
        const int ik_begin = this->_nk * (is - 1); // jump this->_nk for spin_down if nspin==2
        hamilt::HContainer<TR>*const target_DMR = this->dmr[is - 1];
        target_DMR->set_zero();
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic)
#endif
        for (int i = 0; i < target_DMR->size_atom_pairs(); ++i)
        {
            hamilt::AtomPair<TR>& target_ap = target_DMR->get_atom_pair(i);
            const int iat1 = target_ap.get_atom_i();
            const int iat2 = target_ap.get_atom_j();
            const int row_ap = this->pv->atom_begin_row[iat1];
            const int col_ap = this->pv->atom_begin_col[iat2];
            const int row_size = this->pv->get_nrow_atom(iat1);
            const int col_size = this->pv->get_ncol_atom(iat2);
            const int R_size = target_ap.get_R_size();
            assert(row_ap != -1 && col_ap != -1 && "Atom-pair not belong this process");
            assert(R_size == 1);
            const ModuleBase::Vector3<int> R_index = target_ap.get_R_index(0);
            assert(R_index.x == 0 && R_index.y == 0 && R_index.z == 0);
            hamilt::BaseMatrix<TR>*const target_mat = target_ap.find_matrix(R_index);
#ifdef __DEBUG
            if (target_mat == nullptr)
            {
                std::cout << "target_mat is nullptr" << std::endl;
                continue;
            }
#endif
            // k index
            constexpr TK kphase = 1;
            // transpose DMK col=>row
            const TK* DMK_mat_ptr
                = this->dmk[0 + ik_begin].data()
                  + col_ap * this->pv->nrow + row_ap;
            // set DMR element
            TR* target_DMR_ptr = target_mat->get_pointer();
            for (int mu = 0; mu < row_size; ++mu)
            {
                BlasConnector::axpy(col_size,
                                    kphase,
                                    DMK_mat_ptr,
                                    ld_hk,
                                    target_DMR_ptr,
                                    1);
                DMK_mat_ptr += 1;
                target_DMR_ptr += col_size;
            }
        }
    }
    this->_dmr_ready = true;
    ModuleBase::timer::end("DensityMatrix", "cal_dmr");
}

} // namespace module_dm
