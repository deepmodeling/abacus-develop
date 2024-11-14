#include "gint.h"
#include "module_base/memory.h"
#include "module_parameter/parameter.h"
#include "module_base/timer.h"

void Gint::gint_kernel_vlocal(Gint_inout* inout) {
    ModuleBase::TITLE("Gint_interface", "cal_gint_vlocal");
    ModuleBase::timer::tick("Gint_interface", "cal_gint_vlocal");
    const UnitCell& ucell = *this->ucell;
    const int max_size = this->gridt->max_atom;
    const int lgd = this->gridt->lgd;
    const int ncyz = this->ny * this->nplane;
    const double dv = ucell.omega / this->ncxyz;
    const double delta_r = this->gridt->dr_uniform;
    hamilt::HContainer<double>* hRGint_kernel = PARAM.inp.nspin != 4 ? this->hRGint : this->hRGint_tmp[inout->ispin];
    hRGint_kernel->set_zero();

#pragma omp parallel 
    {   /**
        * @brief When in OpenMP, it points to a newly allocated memory,
        */
        hamilt::HContainer<double> hRGint_thread(*hRGint_kernel);
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
            /**
             * @brief Prepare block information
            */
            ModuleBase::Array_Pool<bool> cal_flag(this->bxyz,max_size);

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

            Gint_Tools::get_block_info(*this->gridt, this->bxyz, na_grid, grid_index, 
                                                block_iw.data(), block_index.data(), block_size.data(), cal_flag.get_ptr_2D());

        /**
         * @brief Evaluate psi and dpsi on grids
        */
        const int LD_pool = block_index[na_grid];
        ModuleBase::Array_Pool<double> psir_ylm(this->bxyz, LD_pool);
	    Gint_Tools::cal_psir_ylm(*this->gridt, 
            this->bxyz, na_grid, grid_index, delta_r,
            block_index.data(), block_size.data(), 
            cal_flag.get_ptr_2D(),psir_ylm.get_ptr_2D());

        // psir_ylm_new=psir_func(psir_ylm)
        // psir_func==nullptr means psir_ylm_new=psir_ylm
        const ModuleBase::Array_Pool<double> &psir_ylm_1 = (!this->psir_func_1) ? psir_ylm : this->psir_func_1(psir_ylm, *this->gridt, grid_index, 0, block_iw, block_size, block_index, cal_flag);
        const ModuleBase::Array_Pool<double> &psir_ylm_2 = (!this->psir_func_2) ? psir_ylm : this->psir_func_2(psir_ylm, *this->gridt, grid_index, 0, block_iw, block_size, block_index, cal_flag);

	//calculating f_mu(r) = v(r)*psi_mu(r)*dv
        const ModuleBase::Array_Pool<double> psir_vlbr3 = Gint_Tools::get_psir_vlbr3(
                this->bxyz, na_grid, LD_pool, block_index.data(), 
                cal_flag.get_ptr_2D(), vldr3.data(), psir_ylm_1.get_ptr_2D());

            //integrate (psi_mu*v(r)*dv) * psi_nu on grid
            //and accumulates to the corresponding element in Hamiltonian
            this->cal_meshball_vlocal(
                na_grid, LD_pool, block_iw.data(), block_size.data(), block_index.data(), grid_index, 
                cal_flag.get_ptr_2D(),psir_ylm.get_ptr_2D(), psir_vlbr3.get_ptr_2D(),
                &hRGint_thread);
        }

    #pragma omp critical
        {
            BlasConnector::axpy(hRGint_thread.get_nnr(),
                                1.0,
                                hRGint_thread.get_wrapper(),
                                1,
                                hRGint_kernel->get_wrapper(),
                                1);
        }

        ModuleBase::TITLE("Gint_interface", "cal_gint_vlocal");
        ModuleBase::timer::tick("Gint_interface", "cal_gint_vlocal");
    }
}

void Gint::gint_kernel_dvlocal(Gint_inout* inout) {
    ModuleBase::TITLE("Gint_interface", "cal_gint_dvlocal");
    ModuleBase::timer::tick("Gint_interface", "cal_gint_dvlocal");
    const UnitCell& ucell = *this->ucell;
    const int max_size = this->gridt->max_atom;
    const int lgd = this->gridt->lgd;
    const int nnrg = pvdpRx_reduced[inout->ispin].get_nnr();
    const int ncyz = this->ny * this->nplane;
    const double dv = ucell.omega / this->ncxyz;
    const double delta_r = this->gridt->dr_uniform;

    if (PARAM.globalv.gamma_only_local) {
        ModuleBase::WARNING_QUIT("Gint_interface::cal_gint","dvlocal only for k point!");
    }
    pvdpRx_reduced[inout->ispin].set_zero();
    pvdpRy_reduced[inout->ispin].set_zero();
    pvdpRz_reduced[inout->ispin].set_zero();

#pragma omp parallel 
{
    hamilt::HContainer<double> pvdpRx_thread(pvdpRx_reduced[inout->ispin]);
    hamilt::HContainer<double> pvdpRy_thread(pvdpRy_reduced[inout->ispin]);
    hamilt::HContainer<double> pvdpRz_thread(pvdpRz_reduced[inout->ispin]);
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
                                    block_iw.data(), block_index.data(), block_size.data(), cal_flag.get_ptr_2D());
        
	//evaluate psi and dpsi on grids
        const int LD_pool = block_index[na_grid];

        ModuleBase::Array_Pool<double> psir_ylm(this->bxyz, LD_pool);
        ModuleBase::Array_Pool<double> dpsir_ylm_x(this->bxyz, LD_pool);
        ModuleBase::Array_Pool<double> dpsir_ylm_y(this->bxyz, LD_pool);
        ModuleBase::Array_Pool<double> dpsir_ylm_z(this->bxyz, LD_pool);
        Gint_Tools::cal_dpsir_ylm(*this->gridt, this->bxyz, na_grid, grid_index, delta_r, 
                                    block_index.data(), block_size.data(), cal_flag.get_ptr_2D(),psir_ylm.get_ptr_2D(),
                                    dpsir_ylm_x.get_ptr_2D(), dpsir_ylm_y.get_ptr_2D(), dpsir_ylm_z.get_ptr_2D());

	//calculating f_mu(r) = v(r)*psi_mu(r)*dv
        const ModuleBase::Array_Pool<double> psir_vlbr3 = Gint_Tools::get_psir_vlbr3(
                this->bxyz, na_grid, LD_pool, block_index.data(), cal_flag.get_ptr_2D(), vldr3.data(), psir_ylm.get_ptr_2D());

	//integrate (psi_mu*v(r)*dv) * psi_nu on grid
	//and accumulates to the corresponding element in Hamiltonian
        this->cal_meshball_vlocal(na_grid, LD_pool, block_size.data(), block_index.data(),
                                    block_iw.data(), grid_index, cal_flag.get_ptr_2D(),psir_vlbr3.get_ptr_2D(),
                                    dpsir_ylm_x.get_ptr_2D(), &pvdpRx_thread);
        this->cal_meshball_vlocal(na_grid, LD_pool, block_size.data(), block_index.data(),
                                    block_iw.data(), grid_index, cal_flag.get_ptr_2D(),psir_vlbr3.get_ptr_2D(),
                                    dpsir_ylm_y.get_ptr_2D(), &pvdpRy_thread);
        this->cal_meshball_vlocal(na_grid, LD_pool, block_size.data(), block_index.data(),
                                    block_iw.data(), grid_index, cal_flag.get_ptr_2D(),psir_vlbr3.get_ptr_2D(),
                                    dpsir_ylm_z.get_ptr_2D(), &pvdpRz_thread);
    }
    #pragma omp critical(gint_k)
    {
        BlasConnector::axpy(nnrg,
                            1.0,
                            pvdpRx_thread.get_wrapper(),
                            1,
                            this->pvdpRx_reduced[inout->ispin].get_wrapper(),
                            1);
        BlasConnector::axpy(nnrg,
                            1.0,
                            pvdpRy_thread.get_wrapper(),
                            1,
                            this->pvdpRy_reduced[inout->ispin].get_wrapper(),
                            1);
        BlasConnector::axpy(nnrg,
                            1.0,
                            pvdpRz_thread.get_wrapper(),
                            1,
                            this->pvdpRz_reduced[inout->ispin].get_wrapper(),
                            1);
    }
}
    ModuleBase::TITLE("Gint_interface", "cal_gint_dvlocal");
    ModuleBase::timer::tick("Gint_interface", "cal_gint_dvlocal");
}

void Gint::gint_kernel_vlocal_meta(Gint_inout* inout) {
    ModuleBase::TITLE("Gint_interface", "cal_gint_vlocal_meta");
    ModuleBase::timer::tick("Gint_interface", "cal_gint_vlocal_meta");
    const UnitCell& ucell = *this->ucell;
    const int max_size = this->gridt->max_atom;
    const int lgd = this->gridt->lgd;
    const int ncyz = this->ny * this->nplane;
    const double dv = ucell.omega / this->ncxyz;
    const double delta_r = this->gridt->dr_uniform;
    hamilt::HContainer<double>* hRGint_kernel = PARAM.inp.nspin != 4 ? this->hRGint : this->hRGint_tmp[inout->ispin];
    hRGint_kernel->set_zero();
    const int nnrg = hRGint_kernel->get_nnr();

#pragma omp parallel
{
    // define HContainer here to reference.
    //Under the condition of gamma_only, hRGint will be instantiated.
    hamilt::HContainer<double> hRGint_thread(*hRGint_kernel);
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
            dpsir_ylm_z.get_ptr_2D()
        );
	
	    //calculating f_mu(r) = v(r)*psi_mu(r)*dv
	    const ModuleBase::Array_Pool<double> psir_vlbr3 = Gint_Tools::get_psir_vlbr3(
		    	this->bxyz, na_grid, LD_pool, block_index.data(), cal_flag.get_ptr_2D(), vldr3.data(), psir_ylm.get_ptr_2D());

	    //calculating df_mu(r) = vofk(r) * dpsi_mu(r) * dv
	    const ModuleBase::Array_Pool<double> dpsix_vlbr3 = Gint_Tools::get_psir_vlbr3(
			this->bxyz, na_grid, LD_pool, block_index.data(), cal_flag.get_ptr_2D(), vkdr3.data(), dpsir_ylm_x.get_ptr_2D());
	    const ModuleBase::Array_Pool<double> dpsiy_vlbr3 = Gint_Tools::get_psir_vlbr3(
			this->bxyz, na_grid, LD_pool, block_index.data(), cal_flag.get_ptr_2D(), vkdr3.data(), dpsir_ylm_y.get_ptr_2D());	
	    const ModuleBase::Array_Pool<double> dpsiz_vlbr3 = Gint_Tools::get_psir_vlbr3(
			this->bxyz, na_grid, LD_pool, block_index.data(), cal_flag.get_ptr_2D(), vkdr3.data(), dpsir_ylm_z.get_ptr_2D());


        //integrate (psi_mu*v(r)*dv) * psi_nu on grid
        //and accumulates to the corresponding element in Hamiltonian
        this->cal_meshball_vlocal(
            na_grid, LD_pool, block_iw.data(), block_size.data(), block_index.data(), grid_index, cal_flag.get_ptr_2D(),
            psir_ylm.get_ptr_2D(), psir_vlbr3.get_ptr_2D(), &hRGint_thread);
        //integrate (d/dx_i psi_mu*vk(r)*dv) * (d/dx_i psi_nu) on grid (x_i=x,y,z)
        //and accumulates to the corresponding element in Hamiltonian
        this->cal_meshball_vlocal(
            na_grid, LD_pool, block_iw.data(), block_size.data(), block_index.data(), grid_index, cal_flag.get_ptr_2D(),
            dpsir_ylm_x.get_ptr_2D(), dpsix_vlbr3.get_ptr_2D(), &hRGint_thread);
        this->cal_meshball_vlocal(
            na_grid, LD_pool, block_iw.data(), block_size.data(), block_index.data(), grid_index, cal_flag.get_ptr_2D(),
            dpsir_ylm_y.get_ptr_2D(), dpsiy_vlbr3.get_ptr_2D(), &hRGint_thread);
        this->cal_meshball_vlocal(
            na_grid, LD_pool, block_iw.data(), block_size.data(), block_index.data(), grid_index, cal_flag.get_ptr_2D(),
            dpsir_ylm_z.get_ptr_2D(), dpsiz_vlbr3.get_ptr_2D(), &hRGint_thread);
    }

#pragma omp critical
    {
        BlasConnector::axpy(nnrg,
                            1.0,
                            hRGint_thread.get_wrapper(),
                            1,
                            hRGint_kernel->get_wrapper(),
                            1);
    }
}

    ModuleBase::TITLE("Gint_interface", "cal_gint_vlocal_meta");
    ModuleBase::timer::tick("Gint_interface", "cal_gint_vlocal_meta");
}