#include "gint.h"
#include "module_base/memory.h"
#include "module_base/timer.h"

void Gint::cpu_force_interface(Gint_inout* inout) {
    ModuleBase::TITLE("Gint_interface", "cal_gint_force");
    ModuleBase::timer::tick("Gint_interface", "cal_gint_force");
    const UnitCell& ucell = *this->ucell;
    const int max_size = this->gridt->max_atom;
    const int LD_pool = max_size * ucell.nwmax;
    const int ncyz = this->ny * this->nplane;
    const double dv = ucell.omega / this->ncxyz;
    const double delta_r = this->gridt->dr_uniform;
    ModuleBase::matrix fvl_dphi_thread=*inout->fvl_dphi;
    ModuleBase::matrix svl_dphi_thread=*inout->svl_dphi;
#ifdef _OPENMP
#pragma omp parallel private(fvl_dphi_thread, svl_dphi_thread)
{
    if (inout->isforce) {
        fvl_dphi_thread.create(inout->fvl_dphi->nr, inout->fvl_dphi->nc);
        fvl_dphi_thread.zero_out();
    }
    if (inout->isstress) {
        svl_dphi_thread.create(inout->svl_dphi->nr, inout->svl_dphi->nc);
        svl_dphi_thread.zero_out();
    }
#pragma omp for
#endif
    for (int grid_index = 0; grid_index < this->nbxx; grid_index++) {
        const int na_grid = this->gridt->how_many_atoms[grid_index];
        if (na_grid == 0) {
            continue;
        }
        double* vldr3 = Gint_Tools::get_vldr3(inout->vl,
                                    this->bxyz,
                                    this->bx,
                                    this->by,
                                    this->bz,
                                    this->nplane,
                                    this->gridt->start_ind[grid_index],
                                    ncyz,
                                    dv);
         //prepare block information
        int* block_iw=nullptr;
        int* block_index=nullptr;
        int* block_size=nullptr;
        bool** cal_flag;
        Gint_Tools::get_block_info(*this->gridt, this->bxyz, na_grid, grid_index, block_iw, block_index, block_size, cal_flag);
        int LD_pool = block_index[na_grid];

    //evaluate psi and dpsi on grids
        ModuleBase::Array_Pool<double> psir_ylm(this->bxyz, LD_pool);
        ModuleBase::Array_Pool<double> dpsir_ylm_x(this->bxyz, LD_pool);
        ModuleBase::Array_Pool<double> dpsir_ylm_y(this->bxyz, LD_pool);
        ModuleBase::Array_Pool<double> dpsir_ylm_z(this->bxyz, LD_pool);

        Gint_Tools::cal_dpsir_ylm(*this->gridt, this->bxyz, na_grid, grid_index, delta_r,	block_index, block_size, cal_flag,
            psir_ylm.get_ptr_2D(), dpsir_ylm_x.get_ptr_2D(), dpsir_ylm_y.get_ptr_2D(), dpsir_ylm_z.get_ptr_2D());

    //calculating f_mu(r) = v(r)*psi_mu(r)*dv
        const ModuleBase::Array_Pool<double> psir_vlbr3 
            = Gint_Tools::get_psir_vlbr3(this->bxyz, na_grid, LD_pool, block_index, cal_flag, vldr3, psir_ylm.get_ptr_2D());

        ModuleBase::Array_Pool<double> psir_vlbr3_DM(this->bxyz, LD_pool);
        ModuleBase::GlobalFunc::ZEROS(psir_vlbr3_DM.get_ptr_1D(), this->bxyz*LD_pool);

	//calculating g_mu(r) = sum_nu rho_mu,nu f_nu(r)
        if(GlobalV::GAMMA_ONLY_LOCAL)
        {
		//Gint_Tools::mult_psi_DM(*this->gridt, this->bxyz, na_grid, LD_pool, block_iw, block_size, block_index, cal_flag,
		//	psir_vlbr3.get_ptr_2D(), psir_vlbr3_DM.get_ptr_2D(), DM_in, 2);
            Gint_Tools::mult_psi_DM_new(
                    *this->gridt, 
                    this->bxyz, 
                    grid_index, 
                    na_grid, 
                    LD_pool, 
                    block_iw, 
                    block_size, 
                    block_index, 
                    cal_flag,
                    psir_vlbr3.get_ptr_2D(), 
                    psir_vlbr3_DM.get_ptr_2D(), 
                    this->DMRGint[inout->ispin], 
                    false);
        }
        else
        {
            Gint_Tools::mult_psi_DMR(*this->gridt, this->bxyz, LD_pool, grid_index, na_grid, block_index, 
                                        block_size, cal_flag, psir_vlbr3.get_ptr_2D(), psir_vlbr3_DM.get_ptr_2D(), 
                                        this->DMRGint[inout->ispin], false);
        }

        if(inout->isforce)
        {
            //do integration to get force
            this-> cal_meshball_force(grid_index, na_grid, block_size, block_index,psir_vlbr3_DM.get_ptr_2D(), 
                                        dpsir_ylm_x.get_ptr_2D(), dpsir_ylm_y.get_ptr_2D(), dpsir_ylm_z.get_ptr_2D(),
                                        &fvl_dphi_thread);
        }
        if(inout->isstress)
        {
            //calculating g_mu(r)*(r-R) where R is the location of atom

            // The array dpsirr contains derivatives of psir in the xx, xy, xz, yy, yz, zz directions,
            // with each set of six numbers representing the derivatives in these respective directions.
            ModuleBase::Array_Pool<double> dpsirr_ylm(this->bxyz, LD_pool * 6);
            Gint_Tools::cal_dpsirr_ylm(*this->gridt, this->bxyz, na_grid, grid_index, block_index, 
                                        block_size, cal_flag,dpsir_ylm_x.get_ptr_2D(), dpsir_ylm_y.get_ptr_2D(),
                                        dpsir_ylm_z.get_ptr_2D(),dpsirr_ylm.get_ptr_2D());

            //do integration to get stress
            this-> cal_meshball_stress(na_grid, block_index, psir_vlbr3_DM.get_ptr_1D(), 
                                        dpsirr_ylm.get_ptr_1D(), &svl_dphi_thread);
        }

    //release memories
        delete[] block_iw;
        delete[] block_index;
        delete[] block_size;
        for(int ib=0; ib<this->bxyz; ++ib)
        {
            delete[] cal_flag[ib];
        }
        delete[] cal_flag;
            delete[] vldr3;
        }
#ifdef _OPENMP
#pragma omp critical(gint)
    {
        if (inout->isforce) {
            inout->fvl_dphi[0] += fvl_dphi_thread;
        }
        if (inout->isstress) {
            inout->svl_dphi[0] += svl_dphi_thread;
        }
    }
}
#endif
    ModuleBase::TITLE("Gint_interface", "cal_gint_force");
    ModuleBase::timer::tick("Gint_interface", "cal_gint_force");
}

void Gint::cpu_force_meta_interface(Gint_inout* inout) {
    ModuleBase::TITLE("Gint_interface", "cal_gint_force_meta");
    ModuleBase::timer::tick("Gint_interface", "cal_gint_force_meta");
    const UnitCell& ucell = *this->ucell;
    const int max_size = this->gridt->max_atom;
    const int LD_pool = max_size * ucell.nwmax;
    const int ncyz = this->ny * this->nplane;
    const double dv = ucell.omega / this->ncxyz;
    const double delta_r = this->gridt->dr_uniform;
#ifdef _OPENMP
    ModuleBase::matrix fvl_dphi_thread;
    ModuleBase::matrix svl_dphi_thread;
    if (inout->isforce) {
        fvl_dphi_thread.create(inout->fvl_dphi->nr, inout->fvl_dphi->nc);
        fvl_dphi_thread.zero_out();
    }
    if (inout->isstress) {
        svl_dphi_thread.create(inout->svl_dphi->nr, inout->svl_dphi->nc);
        svl_dphi_thread.zero_out();
    }
#pragma omp for
#endif
    for (int grid_index = 0; grid_index < this->nbxx; grid_index++) {
        const int na_grid = this->gridt->how_many_atoms[grid_index];
        if (na_grid == 0) {
            continue;
        }
        double* vldr3
            = Gint_Tools::get_vldr3(inout->vl,
                                    this->bxyz,
                                    this->bx,
                                    this->by,
                                    this->bz,
                                    this->nplane,
                                    this->gridt->start_ind[grid_index],
                                    ncyz,
                                    dv);

        double* vkdr3
            = Gint_Tools::get_vldr3(inout->vofk,
                                    this->bxyz,
                                    this->bx,
                                    this->by,
                                    this->bz,
                                    this->nplane,
                                    this->gridt->start_ind[grid_index],
                                    ncyz,
                                    dv);
#ifdef _OPENMP
        this->gint_kernel_force_meta(na_grid,
                                     grid_index,
                                     delta_r,
                                     vldr3,
                                     vkdr3,
                                     inout->ispin,
                                     inout->isforce,
                                     inout->isstress,
                                     &fvl_dphi_thread,
                                     &svl_dphi_thread,
                                     ucell);
#else
        this->gint_kernel_force_meta(na_grid,
                                     grid_index,
                                     delta_r,
                                     vldr3,
                                     vkdr3,
                                     inout->ispin,
                                     inout->isforce,
                                     inout->isstress,
                                     inout->fvl_dphi,
                                     inout->svl_dphi,
                                     ucell);
#endif
        delete[] vldr3;
        delete[] vkdr3;
    }
#ifdef _OPENMP
#pragma omp critical(gint)
    {
        if (inout->isforce) {
            inout->fvl_dphi[0] += fvl_dphi_thread;
        }
        if (inout->isstress) {
            inout->svl_dphi[0] += svl_dphi_thread;
        }
    }
#endif
    ModuleBase::TITLE("Gint_interface", "cal_gint_force_meta");
    ModuleBase::timer::tick("Gint_interface", "cal_gint_force_meta");
}
