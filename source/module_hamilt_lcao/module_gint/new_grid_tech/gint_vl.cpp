#include "module_base/array_pool.h"
#include "module_hamilt_lcao/module_hcontainer/hcontainer_funcs.h"
#include "gint_vl.h"
#include "phi_operator.h"

namespace ModuleGint
{

void Gint_vl::cal_gint()
{
    init_hRGint_();
    cal_hRGint_();
    compose_hRGint_();
    transfer_hRGint_to_hR_();
}

//========================
// Private functions
//========================

void Gint_vl::init_hRGint_()
{
    hRGint_ = gint_info_->get_hr<double>();
}

void Gint_vl::cal_hRGint_()
{
    for(const auto& biggrid: gint_info_->get_biggrids())
    {
        if(biggrid->get_atoms().size() == 0)
        {
            continue;
        }
        PhiOperator phi_op(biggrid);
        ModuleBase::Array_Pool<double> phi(phi_op.get_rows(), phi_op.get_cols());
        ModuleBase::Array_Pool<double> phi_vldr3(phi_op.get_rows(), phi_op.get_cols());
        phi_op.set_phi(phi.get_ptr_1D());
        phi_op.phi_mul_vldr3(vr_eff_, dr3_, phi.get_ptr_2D(), phi_vldr3.get_ptr_2D());
        phi_op.phi_mul_phi_vldr3(phi.get_ptr_2D(), phi_vldr3.get_ptr_2D(), hRGint_.get());
    }
}

void Gint_vl::compose_hRGint_()
{
    for (int iap = 0; iap < hRGint_->size_atom_pairs(); iap++)
    {
        auto& ap = hRGint_->get_atom_pair(iap);
        const int iat1 = ap.get_atom_i();
        const int iat2 = ap.get_atom_j();
        if (iat1 > iat2)
        {
            // fill lower triangle matrix with upper triangle matrix
            // the upper <IJR> is <iat2, iat1>
            const hamilt::AtomPair<double>* upper_ap = hRGint_->find_pair(iat2, iat1);
            const hamilt::AtomPair<double>* lower_ap = hRGint_->find_pair(iat1, iat2);
#ifdef __DEBUG
            assert(upper_ap != nullptr);
#endif
            for (int ir = 0; ir < ap.get_R_size(); ir++)
            {   
                auto R_index = ap.get_R_index(ir);
                auto upper_mat = upper_ap->find_matrix(-R_index);
                auto lower_mat = lower_ap->find_matrix(R_index);
                for (int irow = 0; irow < upper_mat->get_row_size(); ++irow)
                {
                    for (int icol = 0; icol < upper_mat->get_col_size(); ++icol)
                    {
                        lower_mat->get_value(icol, irow) = upper_ap->get_value(irow, icol);
                    }
                }
            }
        }
    }
}

void Gint_vl::transfer_hRGint_to_hR_()
{
#ifdef __MPI
    int size = 0;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    if (size == 1)
    {
        hR_->add(*hRGint_);
    }
    else
    {
        hamilt::transferSerials2Parallels(*hRGint_, hR_);
    }
#else
    hR_->add(*hRGint_);
#endif
}

}