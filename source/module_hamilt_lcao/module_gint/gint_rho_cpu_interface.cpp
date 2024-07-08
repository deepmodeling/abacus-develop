#include "gint.h"
#include "module_base/memory.h"
#include "module_base/timer.h"

void Gint::cpu_rho_interface(Gint_inout* inout) {
    ModuleBase::TITLE("Gint_interface", "cal_gint_rho");
    ModuleBase::timer::tick("Gint_interface", "cal_gint_rho");
    const UnitCell& ucell = *this->ucell;
    const int max_size = this->gridt->max_atom;
    const int LD_pool = max_size * ucell.nwmax;
    const int lgd = this->gridt->lgd;
    const int nnrg = this->gridt->nnrg;
    const int ncyz = this->ny * this->nplane;
    const double dv = ucell.omega / this->ncxyz;
    const double delta_r = this->gridt->dr_uniform;
    double* pvpR_thread = nullptr;
    hamilt::HContainer<double>* hRGint_thread = nullptr;
#ifdef _OPENMP
#pragma omp for
#endif
    for (int grid_index = 0; grid_index < this->nbxx; grid_index++) {
        const int na_grid = this->gridt->how_many_atoms[grid_index];
        if (na_grid == 0) {
            continue;
        }
        // int* vindex = Gint_Tools::get_vindex(ncyz, ibx, jby, kbz);
        int* vindex = Gint_Tools::get_vindex(this->bxyz,
                                             this->bx,
                                             this->by,
                                             this->bz,
                                             this->nplane,
                                             this->gridt->start_ind[grid_index],
                                             ncyz);
        this->gint_kernel_rho(na_grid,
                              grid_index,
                              delta_r,
                              vindex,
                              LD_pool,
                              ucell,
                              inout);
        delete[] vindex;
    }

    ModuleBase::TITLE("Gint_interface", "cal_gint_rho");
    ModuleBase::timer::tick("Gint_interface", "cal_gint_rho");
}

void Gint::cpu_tau_interface(Gint_inout* inout) {
    ModuleBase::TITLE("Gint_interface", "cal_gint_tau");
    ModuleBase::timer::tick("Gint_interface", "cal_gint_tau");
    const UnitCell& ucell = *this->ucell;
    const int max_size = this->gridt->max_atom;
    const int LD_pool = max_size * ucell.nwmax;
    const int lgd = this->gridt->lgd;
    const int nnrg = this->gridt->nnrg;
    const int ncyz = this->ny * this->nplane;
    const double dv = ucell.omega / this->ncxyz;
    const double delta_r = this->gridt->dr_uniform;

#ifdef _OPENMP
#pragma omp for
#endif
    for (int grid_index = 0; grid_index < this->nbxx; grid_index++) {
        const int na_grid = this->gridt->how_many_atoms[grid_index];
        if (na_grid == 0) {
            continue;
        }
        // int* vindex = Gint_Tools::get_vindex(ncyz, ibx, jby, kbz);
        int* vindex = Gint_Tools::get_vindex(this->bxyz,
                                             this->bx,
                                             this->by,
                                             this->bz,
                                             this->nplane,
                                             this->gridt->start_ind[grid_index],
                                             ncyz);
        this->gint_kernel_tau(na_grid,
                              grid_index,
                              delta_r,
                              vindex,
                              LD_pool,
                              inout,
                              ucell);
        delete[] vindex;
    }
    ModuleBase::TITLE("Gint_interface", "cal_gint_tau");
    ModuleBase::timer::tick("Gint_interface", "cal_gint_tau");
}
