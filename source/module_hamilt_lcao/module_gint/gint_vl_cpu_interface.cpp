#include "gint.h"
#include "module_base/memory.h"
#include "module_base/timer.h"

void Gint::cpu_vlocal_interface(Gint_inout* inout) {
    ModuleBase::TITLE("Gint_interface", "cal_gint_vlocal");
    ModuleBase::timer::tick("Gint_interface", "cal_gint_vlocal");
    const UnitCell& ucell = *this->ucell;
    const int max_size = this->gridt->max_atom;
    const int LD_pool = max_size * ucell.nwmax;
    const int lgd = this->gridt->lgd;
    const int nnrg = this->gridt->nnrg;
    const int ncyz = this->ny * this->nplane;
    const double dv = ucell.omega / this->ncxyz;
    const double delta_r = this->gridt->dr_uniform;
    if (!GlobalV::GAMMA_ONLY_LOCAL) {
        if (!pvpR_alloc_flag) {
            ModuleBase::WARNING_QUIT("Gint_interface::cal_gint",
                                     "pvpR has not been allocated yet!");
        } else {
            ModuleBase::GlobalFunc::ZEROS(this->pvpR_reduced[inout->ispin],
                                          nnrg);
        }
    }

    // define HContainer here to reference.
    hamilt::HContainer<double>* hRGint_thread =nullptr;
    //Under the condition of gamma_only, hRGint will be instantiated.
    if (GlobalV::GAMMA_ONLY_LOCAL){
        hRGint_thread = new hamilt::HContainer<double>(*this->hRGint);
    }
    //use vector instead of new-delete to avoid memory leak.
    std::vector<double> pvpR_thread = std::vector<double>(nnrg, 0.0);
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
        //prepare block information
	    int * block_iw, * block_index, * block_size;
	    bool** cal_flag;
    	Gint_Tools::get_block_info(*this->gridt, this->bxyz, na_grid, grid_index, block_iw, block_index, block_size, cal_flag);
	
	//evaluate psi and dpsi on grids
	    ModuleBase::Array_Pool<double> psir_ylm(this->bxyz, LD_pool);
	    Gint_Tools::cal_psir_ylm(*this->gridt, 
            this->bxyz, na_grid, grid_index, delta_r,
            block_index, block_size, 
            cal_flag,
            psir_ylm.get_ptr_2D());
        
	//calculating f_mu(r) = v(r)*psi_mu(r)*dv
	const ModuleBase::Array_Pool<double> psir_vlbr3 = Gint_Tools::get_psir_vlbr3(
			this->bxyz, na_grid, LD_pool, block_index, cal_flag, vldr3, psir_ylm.get_ptr_2D());

	//integrate (psi_mu*v(r)*dv) * psi_nu on grid
	//and accumulates to the corresponding element in Hamiltonian
    if(GlobalV::GAMMA_ONLY_LOCAL)
    {
		if(hRGint_thread == nullptr) { hRGint_thread = this->hRGint;}
		this->cal_meshball_vlocal_gamma(
			na_grid, LD_pool, block_iw, block_size, block_index, grid_index, cal_flag,
			psir_ylm.get_ptr_2D(), psir_vlbr3.get_ptr_2D(), hRGint_thread);
    }
    else
    {
        this->cal_meshball_vlocal_k(
            na_grid, LD_pool, grid_index, block_size, block_index, block_iw, cal_flag,
            psir_ylm.get_ptr_2D(), psir_vlbr3.get_ptr_2D(), this->pvpR_reduced[inout->ispin],ucell);
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
    
    if (GlobalV::GAMMA_ONLY_LOCAL) {
        {
            BlasConnector::axpy(this->hRGint->get_nnr(),
                                1.0,
                                hRGint_thread->get_wrapper(),
                                1,
                                this->hRGint->get_wrapper(),
                                1);
        }
    } else {
        {
            BlasConnector::axpy(nnrg,
                                1.0,
                                pvpR_thread.data(),
                                1,
                                pvpR_reduced[inout->ispin],
                                1);
        }
    ModuleBase::TITLE("Gint_interface", "cal_gint_vlocal");
    ModuleBase::timer::tick("Gint_interface", "cal_gint_vlocal");
}
}
void Gint::cpu_dvlocal_interface(Gint_inout* inout) {
    ModuleBase::TITLE("Gint_interface", "cal_gint_dvlocal");
    ModuleBase::timer::tick("Gint_interface", "cal_gint_dvlocal");
    const UnitCell& ucell = *this->ucell;
    const int max_size = this->gridt->max_atom;
    const int LD_pool = max_size * ucell.nwmax;
    const int lgd = this->gridt->lgd;
    const int nnrg = this->gridt->nnrg;
    const int ncyz = this->ny * this->nplane;
    const double dv = ucell.omega / this->ncxyz;
    const double delta_r = this->gridt->dr_uniform;

    if (GlobalV::GAMMA_ONLY_LOCAL) {
        ModuleBase::WARNING_QUIT("Gint_interface::cal_gint",
                                 "dvlocal only for k point!");
        ModuleBase::GlobalFunc::ZEROS(
                        this->pvdpRx_reduced[inout->ispin],
                        nnrg);
        ModuleBase::GlobalFunc::ZEROS(
                        this->pvdpRy_reduced[inout->ispin],
                        nnrg);
        ModuleBase::GlobalFunc::ZEROS(
                        this->pvdpRz_reduced[inout->ispin],
                        nnrg);
    }

#ifdef _OPENMP
    std::vector<double> pvdpRx_thread = std::vector<double>(nnrg, 0.0);
    std::vector<double> pvdpRy_thread = std::vector<double>(nnrg, 0.0);
    std::vector<double> pvdpRz_thread = std::vector<double>(nnrg, 0.0);
    hamilt::HContainer<double>* hRGint_thread = new hamilt::HContainer<double>(*this->hRGint);
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
#ifdef _OPENMP
        this->gint_kernel_dvlocal(na_grid,
                                  grid_index,
                                  delta_r,
                                  vldr3,
                                  LD_pool,
                                  pvdpRx_thread.data(),
                                  pvdpRy_thread.data(),
                                  pvdpRz_thread.data(),
                                  ucell);
#else
        this->gint_kernel_dvlocal(na_grid,
                                  grid_index,
                                  delta_r,
                                  vldr3,
                                  LD_pool,
                                  this->pvdpRx_reduced[inout->ispin],
                                  this->pvdpRy_reduced[inout->ispin],
                                  this->pvdpRz_reduced[inout->ispin],
                                  ucell);
#endif
        delete[] vldr3;
    }
#ifdef _OPENMP
    delete hRGint_thread;
#endif
    ModuleBase::TITLE("Gint_interface", "cal_gint_dvlocal");
    ModuleBase::timer::tick("Gint_interface", "cal_gint_dvlocal");
}

void Gint::cpu_vlocal_meta_interface(Gint_inout* inout) {
    ModuleBase::TITLE("Gint_interface", "cal_gint_vlocal_meta");
    ModuleBase::timer::tick("Gint_interface", "cal_gint_vlocal_meta");
    const UnitCell& ucell = *this->ucell;
    const int max_size = this->gridt->max_atom;
    const int LD_pool = max_size * ucell.nwmax;
    const int lgd = this->gridt->lgd;
    const int nnrg = this->gridt->nnrg;
    const int ncyz = this->ny * this->nplane;
    const double dv = ucell.omega / this->ncxyz;
    const double delta_r = this->gridt->dr_uniform;

#ifdef _OPENMP
    if (!GlobalV::GAMMA_ONLY_LOCAL) {
        if (!pvpR_alloc_flag) {
            ModuleBase::WARNING_QUIT("Gint_interface::cal_gint",
                                        "pvpR has not been allocated yet!");
        } else {
            ModuleBase::GlobalFunc::ZEROS(this->pvpR_reduced[inout->ispin],
                                            nnrg);
        }
    }
    // define HContainer here to reference.
    hamilt::HContainer<double>* hRGint_thread = nullptr;
    //Under the condition of gamma_only, hRGint will be instantiated.
    if (GlobalV::GAMMA_ONLY_LOCAL)
    {
        hRGint_thread =new hamilt::HContainer<double>(*this->hRGint);
    }
    //use vector instead of new-delete to avoid memory leak.
    std::vector<double>pvpR_thread = std::vector<double>(nnrg, 0.0);
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
            this->gint_kernel_vlocal_meta(na_grid,
                                          grid_index,
                                          delta_r,
                                          vldr3,
                                          vkdr3,
                                          LD_pool,
                                          pvpR_thread.data(),
                                          ucell,
                                          hRGint_thread);
#else
        if (GlobalV::GAMMA_ONLY_LOCAL) {
            this->gint_kernel_vlocal_meta(na_grid,
                                          grid_index,
                                          delta_r,
                                          vldr3,
                                          vkdr3,
                                          LD_pool,
                                          nullptr,
                                          ucell);
        } else {
            this->gint_kernel_vlocal_meta(na_grid,
                                          grid_index,
                                          delta_r,
                                          vldr3,
                                          vkdr3,
                                          LD_pool,
                                          this->pvpR_reduced[inout->ispin],
                                          ucell);
        }
#endif
        delete[] vldr3;
        delete[] vkdr3;
    }
#ifdef _OPENMP
    if (GlobalV::GAMMA_ONLY_LOCAL) {
#pragma omp critical(gint_gamma)
        {
            BlasConnector::axpy(this->hRGint->get_nnr(),
                                1.0,
                                hRGint_thread->get_wrapper(),
                                1,
                                this->hRGint->get_wrapper(),
                                1);
        }
    }
    else{
#pragma omp critical(gint_k)
        {
            BlasConnector::axpy(nnrg,
                                1.0,
                                pvpR_thread.data(),
                                1,
                                pvpR_reduced[inout->ispin],
                                1);
        }
    }
    delete hRGint_thread;
#endif
    ModuleBase::TITLE("Gint_interface", "cal_gint_vlocal_meta");
    ModuleBase::timer::tick("Gint_interface", "cal_gint_vlocal_meta");
}