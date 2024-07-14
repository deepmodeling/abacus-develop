#include "gint_k.h"

void Gint::gint_kernel_rho(const int na_grid,
                           const int grid_index,
                           const double delta_r,
                           int* vindex,
                           const int LD_pool,
                           const UnitCell& ucell,
                           Gint_inout* inout)
{
    // prepare block information
    int *block_iw, *block_index, *block_size;
    bool** cal_flag;
    Gint_Tools::get_block_info(*this->gridt,
                               this->bxyz,
                               na_grid,
                               grid_index,
                               block_iw,
                               block_index,
                               block_size,
                               cal_flag);

    // evaluate psi on grids
    ModuleBase::Array_Pool<double> psir_ylm(this->bxyz, LD_pool);
    Gint_Tools::cal_psir_ylm(*this->gridt,
                             this->bxyz,
                             na_grid,
                             grid_index,
                             delta_r,
                             block_index,
                             block_size,
                             cal_flag,
                             psir_ylm.get_ptr_2D());

    for (int is = 0; is < this->gridt->nspin; ++is)
    {
        ModuleBase::Array_Pool<double> psir_DM(this->bxyz, LD_pool);
        Gint_Func::ZEROS(psir_DM.get_ptr_1D(), this->bxyz * LD_pool);
        if (this->gridt->gamma_only_local)

        {
            Gint_Tools::mult_psi_DM_new(*this->gridt,
                                        this->bxyz,
                                        grid_index,
                                        na_grid,
                                        LD_pool,
                                        block_iw,
                                        block_size,
                                        block_index,
                                        cal_flag,
                                        psir_ylm.get_ptr_2D(),
                                        psir_DM.get_ptr_2D(),
                                        this->DMRGint[is],
                                        inout->if_symm);
        }
        else
        {
            // calculating g_mu(r) = sum_nu rho_mu,nu psi_nu(r)
            Gint_Tools::mult_psi_DMR(*this->gridt,
                                     this->bxyz,
                                     grid_index,
                                     na_grid,
                                     block_index,
                                     block_size,
                                     cal_flag,
                                     psir_ylm.get_ptr_2D(),
                                     psir_DM.get_ptr_2D(),
                                     this->DMRGint[is],
                                     inout->if_symm);
        }

        // do sum_mu g_mu(r)psi_mu(r) to get electron density on grid
        this->cal_meshball_rho(na_grid, block_index, vindex, psir_ylm.get_ptr_2D(), psir_DM.get_ptr_2D(), inout->rho[is]);
    }
    delete[] block_iw;
    delete[] block_index;
    delete[] block_size;
    for (int ib = 0; ib < this->bxyz; ++ib)
    {
        delete[] cal_flag[ib];
    }
    delete[] cal_flag;
}

void Gint::cal_meshball_rho(const int na_grid,
                            int* block_index,
                            int* vindex,
                            double** psir_ylm,
                            double** psir_DMR,
                            double* rho)
{
    const int inc = 1;
    // sum over mu to get density on grid
    for (int ib = 0; ib < this->bxyz; ++ib)
    {
        double r = ddot_(&block_index[na_grid], psir_ylm[ib], &inc, psir_DMR[ib], &inc);
        const int grid = vindex[ib];
        rho[grid] += r;
    }
}
