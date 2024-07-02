#include "gint.h"
#include "module_base/memory.h"
#include "module_base/timer.h"

void Gint::cpu_force_interface(Gint_inout* inout) {
    ModuleBase::TITLE("Gint_interface", "cal_gint_force");
    ModuleBase::timer::tick("Gint_interface", "cal_gint_force");
    const UnitCell& ucell = *this->ucell;
    const int max_size = this->gridt->max_atom;
    const int LD_pool = max_size * ucell.nwmax;
    const int lgd = this->gridt->lgd;
    const int nnrg = this->gridt->nnrg;
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
#ifdef _OPENMP
        this->gint_kernel_force(na_grid,
                                grid_index,
                                delta_r,
                                vldr3,
                                LD_pool,
                                inout->ispin,
                                inout->isforce,
                                inout->isstress,
                                &fvl_dphi_thread,
                                &svl_dphi_thread,
                                ucell);
#else
        this->gint_kernel_force(na_grid,
                                grid_index,
                                delta_r,
                                vldr3,
                                LD_pool,
                                inout->ispin,
                                inout->isforce,
                                inout->isstress,
                                inout->fvl_dphi,
                                inout->svl_dphi,
                                ucell);
#endif
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
    const int lgd = this->gridt->lgd;
    const int nnrg = this->gridt->nnrg;
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
                                     LD_pool,
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
                                     LD_pool,
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
