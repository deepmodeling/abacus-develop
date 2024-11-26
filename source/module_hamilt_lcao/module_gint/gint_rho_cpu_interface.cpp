#include "gint.h"
#include "module_base/memory.h"
#include "module_parameter/parameter.h"
#include "module_base/timer.h"

void Gint::gint_kernel_rho(Gint_inout* inout) {
    ModuleBase::TITLE("Gint_interface", "cal_gint_rho");
    ModuleBase::timer::tick("Gint_interface", "cal_gint_rho");
    const int max_size = this->gridt->max_atom;
    const int ncyz = this->ny * this->nplane;
    const double delta_r = this->gridt->dr_uniform;

#pragma omp parallel
{
    std::vector<int> block_iw(max_size, 0);
    std::vector<int> block_index(max_size+1, 0);
    std::vector<int> block_size(max_size, 0);
    std::vector<int> vindex(this->bxyz, 0);
#pragma omp for
    for (int grid_index = 0; grid_index < this->nbxx; grid_index++)
    {
        const int na_grid = this->gridt->how_many_atoms[grid_index];
        if (na_grid == 0) {
            continue;
        }
        Gint_Tools::get_vindex(this->bxyz,
                                    this->bx,
                                    this->by,
                                    this->bz,
                                    this->nplane,
                                    this->gridt->start_ind[grid_index],
                                    ncyz,
                                    vindex.data());
         // prepare block information
        ModuleBase::Array_Pool<bool> cal_flag(this->bxyz,max_size);
        Gint_Tools::get_block_info(*this->gridt,
                                this->bxyz,
                                na_grid,
                                grid_index,
                                block_iw.data(),
                                block_index.data(),
                                block_size.data(),
                                cal_flag.get_ptr_2D());

        // evaluate psi on grids
        const int LD_pool = block_index[na_grid];
        ModuleBase::Array_Pool<double> psir_ylm(this->bxyz, LD_pool);
        Gint_Tools::cal_psir_ylm(*this->gridt,
                                this->bxyz,
                                na_grid,
                                grid_index,
                                delta_r,
                                block_index.data(),
                                block_size.data(),
                                cal_flag.get_ptr_2D(),
                                psir_ylm.get_ptr_2D());

        for (int is = 0; is < inout->nspin_rho; ++is)
        {
            // psir_ylm_new = psir_func(psir_ylm)
            // psir_func==nullptr means psir_ylm_new=psir_ylm
            const ModuleBase::Array_Pool<double> &psir_ylm_1 = (!this->psir_func_1) ? psir_ylm : this->psir_func_1(psir_ylm, *this->gridt, grid_index, is, block_iw, block_size, block_index, cal_flag);
            const ModuleBase::Array_Pool<double> &psir_ylm_2 = (!this->psir_func_2) ? psir_ylm : this->psir_func_2(psir_ylm, *this->gridt, grid_index, is, block_iw, block_size, block_index, cal_flag);

            ModuleBase::Array_Pool<double> psir_DM(this->bxyz, LD_pool);
            ModuleBase::GlobalFunc::ZEROS(psir_DM.get_ptr_1D(), this->bxyz * LD_pool);

            // calculating g_mu(r) = sum_nu rho_mu,nu psi_nu(r)
            Gint_Tools::mult_psi_DMR(*this->gridt,
                                    this->bxyz,
                                    LD_pool,
                                    grid_index,
                                    na_grid,
                                    block_index.data(),
                                    block_size.data(),
                                    cal_flag.get_ptr_2D(),
                                    psir_ylm_1.get_ptr_2D(),
                                    psir_DM.get_ptr_2D(),
                                    this->DMRGint[is],
                                    inout->if_symm);

            // do sum_mu g_mu(r)psi_mu(r) to get electron density on grid
            this->cal_meshball_rho(na_grid, block_index.data(), vindex.data(), psir_ylm_2.get_ptr_2D(), psir_DM.get_ptr_2D(), inout->rho[is]);
        }
    }
}
    ModuleBase::TITLE("Gint_interface", "cal_gint_rho");
    ModuleBase::timer::tick("Gint_interface", "cal_gint_rho");
}

void Gint::gint_kernel_tau(Gint_inout* inout) {
    ModuleBase::TITLE("Gint_interface", "cal_gint_tau");
    ModuleBase::timer::tick("Gint_interface", "cal_gint_tau");
    const int max_size = this->gridt->max_atom;
    const int ncyz = this->ny * this->nplane;
    const double delta_r = this->gridt->dr_uniform;


#pragma omp parallel
{
    std::vector<int> block_iw(max_size, 0);
    std::vector<int> block_index(max_size+1, 0);
    std::vector<int> block_size(max_size, 0);
    std::vector<int> vindex(bxyz, 0);
#pragma omp for
    for (int grid_index = 0; grid_index < this->nbxx; grid_index++)
    {
        const int na_grid = this->gridt->how_many_atoms[grid_index];
        if (na_grid == 0) {
            continue;
        }
        Gint_Tools::get_vindex(this->bxyz,
                                this->bx,
                                this->by,
                                this->bz,
                                this->nplane,
                                this->gridt->start_ind[grid_index],
                                ncyz,
                                vindex.data());
        //prepare block information
        ModuleBase::Array_Pool<bool> cal_flag(this->bxyz,max_size);
        Gint_Tools::get_block_info(*this->gridt, this->bxyz, na_grid, grid_index,
                                            block_iw.data(), block_index.data(), block_size.data(), cal_flag.get_ptr_2D());

        //evaluate psi and dpsi on grids
        const int LD_pool = block_index[na_grid];
        ModuleBase::Array_Pool<double> psir_ylm(this->bxyz, LD_pool);
        ModuleBase::Array_Pool<double> dpsir_ylm_x(this->bxyz, LD_pool);
        ModuleBase::Array_Pool<double> dpsir_ylm_y(this->bxyz, LD_pool);
        ModuleBase::Array_Pool<double> dpsir_ylm_z(this->bxyz, LD_pool);

        Gint_Tools::cal_dpsir_ylm(*this->gridt,
            this->bxyz, na_grid, grid_index, delta_r,
            block_index.data(), block_size.data(),
            cal_flag.get_ptr_2D(),
            psir_ylm.get_ptr_2D(),
            dpsir_ylm_x.get_ptr_2D(),
            dpsir_ylm_y.get_ptr_2D(),
            dpsir_ylm_z.get_ptr_2D());

        for(int is=0; is<PARAM.inp.nspin; ++is)
        {
            ModuleBase::Array_Pool<double> dpsix_DM(this->bxyz, LD_pool);
            ModuleBase::Array_Pool<double> dpsiy_DM(this->bxyz, LD_pool);
            ModuleBase::Array_Pool<double> dpsiz_DM(this->bxyz, LD_pool);
            ModuleBase::GlobalFunc::ZEROS(dpsix_DM.get_ptr_1D(), this->bxyz*LD_pool);
            ModuleBase::GlobalFunc::ZEROS(dpsiy_DM.get_ptr_1D(), this->bxyz*LD_pool);
            ModuleBase::GlobalFunc::ZEROS(dpsiz_DM.get_ptr_1D(), this->bxyz*LD_pool);

            //calculating g_i,mu(r) = sum_nu rho_mu,nu d/dx_i psi_nu(r), x_i=x,y,z
            Gint_Tools::mult_psi_DMR(
                *this->gridt, this->bxyz,
                LD_pool,
                grid_index, na_grid,
                block_index.data(), block_size.data(),
                cal_flag.get_ptr_2D(),
                dpsir_ylm_x.get_ptr_2D(),
                dpsix_DM.get_ptr_2D(),
                this->DMRGint[is],
                true);
            Gint_Tools::mult_psi_DMR(
                *this->gridt, this->bxyz,
                LD_pool,
                grid_index, na_grid,
                block_index.data(), block_size.data(),
                cal_flag.get_ptr_2D(),
                dpsir_ylm_y.get_ptr_2D(),
                dpsiy_DM.get_ptr_2D(),
                this->DMRGint[is],
                true);
            Gint_Tools::mult_psi_DMR(
                *this->gridt, this->bxyz,
                LD_pool,
                grid_index, na_grid,
                block_index.data(), block_size.data(),
                cal_flag.get_ptr_2D(),
                dpsir_ylm_z.get_ptr_2D(),
                dpsiz_DM.get_ptr_2D(),
                this->DMRGint[is],
                true);

            //do sum_i,mu g_i,mu(r) * d/dx_i psi_mu(r) to get kinetic energy density on grid
            if(inout->job==Gint_Tools::job_type::tau)
            {
                this->cal_meshball_tau(
                    na_grid, block_index.data(),
                    vindex.data(),
                    dpsir_ylm_x.get_ptr_2D(), dpsir_ylm_y.get_ptr_2D(), dpsir_ylm_z.get_ptr_2D(),
                    dpsix_DM.get_ptr_2D(), dpsiy_DM.get_ptr_2D(), dpsiz_DM.get_ptr_2D(),
                    inout->rho[is]);
            }
        }
    }
}
    ModuleBase::TITLE("Gint_interface", "cal_gint_tau");
    ModuleBase::timer::tick("Gint_interface", "cal_gint_tau");
}
