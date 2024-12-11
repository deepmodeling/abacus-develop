#include "gint.h"
#include "module_base/memory.h"
#include "module_base/timer.h"

void Gint::gint_kernel_force(Gint_inout* inout) {
    ModuleBase::TITLE("Gint_interface", "cal_gint_force");
    ModuleBase::timer::tick("Gint_interface", "cal_gint_force");
    const UnitCell& ucell = *this->ucell;
    const int max_size = this->gridt->max_atom;
    const int ncyz = this->ny * this->nplane;
    const double dv = ucell.omega / this->ncxyz;
    const double delta_r = this->gridt->dr_uniform;


#pragma omp parallel 
{
    ModuleBase::matrix* fvl_dphi_thread=inout->fvl_dphi;
    ModuleBase::matrix* svl_dphi_thread=inout->svl_dphi;
    if (inout->isforce) {
        fvl_dphi_thread=new ModuleBase::matrix(*inout->fvl_dphi);
        fvl_dphi_thread->zero_out();
    }
    if (inout->isstress) {
        svl_dphi_thread=new ModuleBase::matrix(*inout->svl_dphi);
        svl_dphi_thread->zero_out();
    }
    std::vector<int> block_iw(max_size,0);
    std::vector<int> block_index(max_size+1,0);
    std::vector<int> block_size(max_size,0);
    std::vector<double> vldr3(this->bxyz,0.0);
#pragma omp for
    for (int grid_index = 0; grid_index < this->nbxx; grid_index++) {
        const int na_grid = this->gridt->how_many_atoms[grid_index];
        if (na_grid == 0) {
            continue;
        }
        Gint_Tools::get_gint_vldr3(vldr3.data(),
                                    inout->vl,
                                    this->bxyz,
                                    this->bx,
                                    this->by,
                                    this->bz,
                                    this->nplane,
                                    this->gridt->start_ind[grid_index],
                                    ncyz,
                                    dv);
         //prepare block information
        ModuleBase::Array_Pool<bool> cal_flag(this->bxyz,max_size);
        Gint_Tools::get_block_info(*this->gridt, this->bxyz, na_grid, grid_index,
                                            block_iw.data(), block_index.data(), block_size.data(), 
                                            cal_flag.get_ptr_2D());
        const int LD_pool = block_index[na_grid];

    //evaluate psi and dpsi on grids
        ModuleBase::Array_Pool<double> psir_ylm(this->bxyz, LD_pool);
        ModuleBase::Array_Pool<double> dpsir_ylm_x(this->bxyz, LD_pool);
        ModuleBase::Array_Pool<double> dpsir_ylm_y(this->bxyz, LD_pool);
        ModuleBase::Array_Pool<double> dpsir_ylm_z(this->bxyz, LD_pool);

        Gint_Tools::cal_dpsir_ylm(*this->gridt, this->bxyz, na_grid, grid_index, delta_r,	
                                    block_index.data(), block_size.data(),
                                    cal_flag.get_ptr_2D(),psir_ylm.get_ptr_2D(),
                                    dpsir_ylm_x.get_ptr_2D(), dpsir_ylm_y.get_ptr_2D(), dpsir_ylm_z.get_ptr_2D());

    //calculating f_mu(r) = v(r)*psi_mu(r)*dv
        const ModuleBase::Array_Pool<double> psir_vlbr3 = 
                Gint_Tools::get_psir_vlbr3(this->bxyz, na_grid, LD_pool, block_index.data(), 
                cal_flag.get_ptr_2D(), vldr3.data(), psir_ylm.get_ptr_2D());

        ModuleBase::Array_Pool<double> psir_vlbr3_DM(this->bxyz, LD_pool);
        ModuleBase::GlobalFunc::ZEROS(psir_vlbr3_DM.get_ptr_1D(), this->bxyz*LD_pool);

	//calculating g_mu(r) = sum_nu rho_mu,nu f_nu(r)
        Gint_Tools::mult_psi_DMR(
                *this->gridt, 
                this->bxyz,
                LD_pool, 
                grid_index, 
                na_grid, 
                block_index.data(), 
                block_size.data(), 
                cal_flag.get_ptr_2D(),
                psir_vlbr3.get_ptr_2D(), 
                psir_vlbr3_DM.get_ptr_2D(), 
                this->DMRGint[inout->ispin], 
                false);

        if(inout->isforce)
        {
            //do integration to get force
            this-> cal_meshball_force(grid_index, na_grid, block_size.data(), block_index.data(),
                                        psir_vlbr3_DM.get_ptr_2D(), dpsir_ylm_x.get_ptr_2D(),
                                        dpsir_ylm_y.get_ptr_2D(), dpsir_ylm_z.get_ptr_2D(),
                                        fvl_dphi_thread);
        }
        if(inout->isstress)
        {
            //calculating g_mu(r)*(r-R) where R is the location of atom

            // The array dpsirr contains derivatives of psir in the xx, xy, xz, yy, yz, zz directions,
            // with each set of six numbers representing the derivatives in these respective directions.
            ModuleBase::Array_Pool<double> dpsirr_ylm(this->bxyz, LD_pool * 6);
            Gint_Tools::cal_dpsirr_ylm(*this->gridt, this->bxyz, na_grid, grid_index, block_index.data(), 
                                        block_size.data(), cal_flag.get_ptr_2D(),dpsir_ylm_x.get_ptr_2D(), 
                                        dpsir_ylm_y.get_ptr_2D(),dpsir_ylm_z.get_ptr_2D(),
                                        dpsirr_ylm.get_ptr_2D());

            //do integration to get stress
            this-> cal_meshball_stress(na_grid, block_index.data(), psir_vlbr3_DM.get_ptr_1D(), 
                                       dpsirr_ylm.get_ptr_1D(), svl_dphi_thread);
        }
    }
#pragma omp critical(gint)
    {
        if (inout->isforce) {
            inout->fvl_dphi[0] += fvl_dphi_thread[0];
            delete fvl_dphi_thread;
        }
        if (inout->isstress) {
            inout->svl_dphi[0] += svl_dphi_thread[0];
            delete svl_dphi_thread;
        }
    }
}
    ModuleBase::TITLE("Gint_interface", "cal_gint_force");
    ModuleBase::timer::tick("Gint_interface", "cal_gint_force");
}

void Gint::gint_kernel_force_meta(Gint_inout* inout) {
    ModuleBase::TITLE("Gint_interface", "cal_gint_force_meta");
    ModuleBase::timer::tick("Gint_interface", "cal_gint_force_meta");
    const UnitCell& ucell = *this->ucell;
    const int max_size = this->gridt->max_atom;
    const int ncyz = this->ny * this->nplane;
    const double dv = ucell.omega / this->ncxyz;
    const double delta_r = this->gridt->dr_uniform;


#pragma omp parallel 
{
    ModuleBase::matrix* fvl_dphi_thread=inout->fvl_dphi;
    ModuleBase::matrix* svl_dphi_thread=inout->svl_dphi;
    if (inout->isforce) {
        fvl_dphi_thread=new ModuleBase::matrix(*inout->fvl_dphi);
        fvl_dphi_thread->zero_out();
    }
    if (inout->isstress) {
        svl_dphi_thread=new ModuleBase::matrix(*inout->svl_dphi);
        svl_dphi_thread->zero_out();
    }
    std::vector<int> block_iw(max_size,0);
    std::vector<int> block_index(max_size+1,0);
    std::vector<int> block_size(max_size,0);
    std::vector<double> vldr3(this->bxyz,0.0);
    std::vector<double> vkdr3(this->bxyz,0.0);
#pragma omp for
    for (int grid_index = 0; grid_index < this->nbxx; grid_index++) {
        const int na_grid = this->gridt->how_many_atoms[grid_index];
        if (na_grid == 0) {
            continue;
        }
        Gint_Tools::get_gint_vldr3(vldr3.data(),
                                    inout->vl,
                                    this->bxyz,
                                    this->bx,
                                    this->by,
                                    this->bz,
                                    this->nplane,
                                    this->gridt->start_ind[grid_index],
                                    ncyz,
                                    dv);

        Gint_Tools::get_gint_vldr3(vkdr3.data(),
                                    inout->vofk,
                                    this->bxyz,
                                    this->bx,
                                    this->by,
                                    this->bz,
                                    this->nplane,
                                    this->gridt->start_ind[grid_index],
                                    ncyz,
                                    dv);
         //prepare block information
        ModuleBase::Array_Pool<bool> cal_flag(this->bxyz,max_size);
        Gint_Tools::get_block_info(*this->gridt, this->bxyz, na_grid, grid_index, 
                                            block_iw.data(), block_index.data(), block_size.data(), cal_flag.get_ptr_2D());
        const int LD_pool = block_index[na_grid];

    //evaluate psi and dpsi on grids
        ModuleBase::Array_Pool<double> psir_ylm(this->bxyz, LD_pool);
        ModuleBase::Array_Pool<double> dpsir_ylm_x(this->bxyz, LD_pool);
        ModuleBase::Array_Pool<double> dpsir_ylm_y(this->bxyz, LD_pool);
        ModuleBase::Array_Pool<double> dpsir_ylm_z(this->bxyz, LD_pool);
        ModuleBase::Array_Pool<double> ddpsir_ylm_xx(this->bxyz, LD_pool);
        ModuleBase::Array_Pool<double> ddpsir_ylm_xy(this->bxyz, LD_pool);
        ModuleBase::Array_Pool<double> ddpsir_ylm_xz(this->bxyz, LD_pool);
        ModuleBase::Array_Pool<double> ddpsir_ylm_yy(this->bxyz, LD_pool);
        ModuleBase::Array_Pool<double> ddpsir_ylm_yz(this->bxyz, LD_pool);
        ModuleBase::Array_Pool<double> ddpsir_ylm_zz(this->bxyz, LD_pool);

	//psi and gradient of psi
        Gint_Tools::cal_dpsir_ylm(*this->gridt, this->bxyz, na_grid, grid_index, delta_r,	block_index.data(), block_size.data(), cal_flag.get_ptr_2D(),
            psir_ylm.get_ptr_2D(), dpsir_ylm_x.get_ptr_2D(), dpsir_ylm_y.get_ptr_2D(), dpsir_ylm_z.get_ptr_2D());

	//hessian of psi
        Gint_Tools::cal_ddpsir_ylm(*this->gridt, this->bxyz, na_grid, grid_index, delta_r, block_index.data(), block_size.data(), cal_flag.get_ptr_2D(),
            ddpsir_ylm_xx.get_ptr_2D(), ddpsir_ylm_xy.get_ptr_2D(), ddpsir_ylm_xz.get_ptr_2D(),
            ddpsir_ylm_yy.get_ptr_2D(), ddpsir_ylm_yz.get_ptr_2D(), ddpsir_ylm_zz.get_ptr_2D());

    //calculating f_mu(r) = v(r)*psi_mu(r)*dv 
        const ModuleBase::Array_Pool<double> psir_vlbr3 
            = Gint_Tools::get_psir_vlbr3(this->bxyz, na_grid, LD_pool, block_index.data(), cal_flag.get_ptr_2D(), vldr3.data(), psir_ylm.get_ptr_2D());
        const ModuleBase::Array_Pool<double> dpsir_x_vlbr3 
            = Gint_Tools::get_psir_vlbr3(this->bxyz, na_grid, LD_pool, block_index.data(), cal_flag.get_ptr_2D(), vkdr3.data(), dpsir_ylm_x.get_ptr_2D());
        const ModuleBase::Array_Pool<double> dpsir_y_vlbr3 
            = Gint_Tools::get_psir_vlbr3(this->bxyz, na_grid, LD_pool, block_index.data(), cal_flag.get_ptr_2D(), vkdr3.data(), dpsir_ylm_y.get_ptr_2D());
        const ModuleBase::Array_Pool<double> dpsir_z_vlbr3 
            = Gint_Tools::get_psir_vlbr3(this->bxyz, na_grid, LD_pool, block_index.data(), cal_flag.get_ptr_2D(), vkdr3.data(), dpsir_ylm_z.get_ptr_2D());

        ModuleBase::Array_Pool<double> psir_vlbr3_DM(this->bxyz, LD_pool);
        ModuleBase::Array_Pool<double> dpsirx_v_DM(this->bxyz, LD_pool);
        ModuleBase::Array_Pool<double> dpsiry_v_DM(this->bxyz, LD_pool);
        ModuleBase::Array_Pool<double> dpsirz_v_DM(this->bxyz, LD_pool);

        ModuleBase::GlobalFunc::ZEROS(psir_vlbr3_DM.get_ptr_1D(), this->bxyz*LD_pool);
        ModuleBase::GlobalFunc::ZEROS(dpsirx_v_DM.get_ptr_1D(), this->bxyz*LD_pool);
        ModuleBase::GlobalFunc::ZEROS(dpsiry_v_DM.get_ptr_1D(), this->bxyz*LD_pool);
        ModuleBase::GlobalFunc::ZEROS(dpsirz_v_DM.get_ptr_1D(), this->bxyz*LD_pool);

	//calculating g_mu(r) = sum_nu rho_mu,nu f_nu(r)
        Gint_Tools::mult_psi_DMR(*this->gridt, this->bxyz, LD_pool, grid_index, 
            na_grid, block_index.data(), block_size.data(), cal_flag.get_ptr_2D(),
            psir_vlbr3.get_ptr_2D(), psir_vlbr3_DM.get_ptr_2D(), this->DMRGint[inout->ispin], false);

        Gint_Tools::mult_psi_DMR(*this->gridt, this->bxyz, LD_pool, grid_index, 
            na_grid, block_index.data(), block_size.data(), cal_flag.get_ptr_2D(),
            dpsir_x_vlbr3.get_ptr_2D(), dpsirx_v_DM.get_ptr_2D(), this->DMRGint[inout->ispin], false);

        Gint_Tools::mult_psi_DMR(*this->gridt, this->bxyz, LD_pool, grid_index,
            na_grid, block_index.data(), block_size.data(), cal_flag.get_ptr_2D(),
            dpsir_y_vlbr3.get_ptr_2D(), dpsiry_v_DM.get_ptr_2D(), this->DMRGint[inout->ispin], false);

        Gint_Tools::mult_psi_DMR(*this->gridt, this->bxyz, LD_pool, grid_index, 
            na_grid, block_index.data(), block_size.data(), cal_flag.get_ptr_2D(),
            dpsir_z_vlbr3.get_ptr_2D(), dpsirz_v_DM.get_ptr_2D(), this->DMRGint[inout->ispin], false);

        if(inout->isforce)
        {
            //do integration to get force
            this-> cal_meshball_force(grid_index, na_grid, block_size.data(), block_index.data(),
                psir_vlbr3_DM.get_ptr_2D(), dpsir_ylm_x.get_ptr_2D(), dpsir_ylm_y.get_ptr_2D(), dpsir_ylm_z.get_ptr_2D(), 
                fvl_dphi_thread);
                
            this-> cal_meshball_force(grid_index, na_grid, block_size.data(), block_index.data(),
                dpsirx_v_DM.get_ptr_2D(), ddpsir_ylm_xx.get_ptr_2D(), ddpsir_ylm_xy.get_ptr_2D(), ddpsir_ylm_xz.get_ptr_2D(), 
                fvl_dphi_thread);
            this-> cal_meshball_force(grid_index, na_grid, block_size.data(), block_index.data(),
                dpsiry_v_DM.get_ptr_2D(), ddpsir_ylm_xy.get_ptr_2D(), ddpsir_ylm_yy.get_ptr_2D(), ddpsir_ylm_yz.get_ptr_2D(), 
                fvl_dphi_thread);
            this-> cal_meshball_force(grid_index, na_grid, block_size.data(), block_index.data(),
                dpsirz_v_DM.get_ptr_2D(), ddpsir_ylm_xz.get_ptr_2D(), ddpsir_ylm_yz.get_ptr_2D(), ddpsir_ylm_zz.get_ptr_2D(), 
                fvl_dphi_thread);		
            
        }
        if(inout->isstress)
        {
            //calculating g_mu(r)*(r-R) where R is the location of atom
            ModuleBase::Array_Pool<double> array(this->bxyz, LD_pool * 6);

            //the vxc part
            Gint_Tools::cal_dpsirr_ylm(*this->gridt, this->bxyz, na_grid, grid_index, block_index.data(), block_size.data(), cal_flag.get_ptr_2D(),
                dpsir_ylm_x.get_ptr_2D(), dpsir_ylm_y.get_ptr_2D(),	dpsir_ylm_z.get_ptr_2D(), array.get_ptr_2D());
            //do integration to get stress
            this-> cal_meshball_stress(na_grid, block_index.data(), psir_vlbr3_DM.get_ptr_1D(),
                array.get_ptr_1D(), svl_dphi_thread);

            //partial x of vtau part
            Gint_Tools::cal_dpsirr_ylm(*this->gridt, this->bxyz, na_grid, grid_index, block_index.data(), block_size.data(), cal_flag.get_ptr_2D(),
                ddpsir_ylm_xx.get_ptr_2D(), ddpsir_ylm_xy.get_ptr_2D(),	ddpsir_ylm_xz.get_ptr_2D(), array.get_ptr_2D());
            //do integration to get stress
            this-> cal_meshball_stress(na_grid, block_index.data(), dpsirx_v_DM.get_ptr_1D(),
                array.get_ptr_1D(), svl_dphi_thread);

            //partial y of vtau part
            Gint_Tools::cal_dpsirr_ylm(*this->gridt, this->bxyz, na_grid, grid_index, block_index.data(), block_size.data(), cal_flag.get_ptr_2D(),
                ddpsir_ylm_xy.get_ptr_2D(), ddpsir_ylm_yy.get_ptr_2D(),	ddpsir_ylm_yz.get_ptr_2D(), array.get_ptr_2D());
            //do integration to get stress
            this-> cal_meshball_stress(na_grid, block_index.data(), dpsiry_v_DM.get_ptr_1D(),
                array.get_ptr_1D(), svl_dphi_thread);

            //partial z of vtau part
            Gint_Tools::cal_dpsirr_ylm(*this->gridt, this->bxyz, na_grid, grid_index, block_index.data(), block_size.data(), cal_flag.get_ptr_2D(),
                ddpsir_ylm_xz.get_ptr_2D(), ddpsir_ylm_yz.get_ptr_2D(), ddpsir_ylm_zz.get_ptr_2D(), array.get_ptr_2D());
            //do integration to get stress
            this-> cal_meshball_stress(na_grid, block_index.data(), dpsirz_v_DM.get_ptr_1D(),
                array.get_ptr_1D(), svl_dphi_thread);
        }
    }
#pragma omp critical(gint)
    {
        if (inout->isforce) {
            inout->fvl_dphi[0] += fvl_dphi_thread[0];
            delete fvl_dphi_thread;
        }
        if (inout->isstress) {
            inout->svl_dphi[0] += svl_dphi_thread[0];
            delete svl_dphi_thread;
        }
    }
}
    ModuleBase::TITLE("Gint_interface", "cal_gint_force_meta");
    ModuleBase::timer::tick("Gint_interface", "cal_gint_force_meta");
}
