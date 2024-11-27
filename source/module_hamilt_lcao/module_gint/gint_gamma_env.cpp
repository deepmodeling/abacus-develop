#include "gint_gamma.h"
#include "grid_technique.h"
#include "module_base/timer.h"
#include "module_base/ylm.h"
#include "module_base/array_pool.h"
#include "module_basis/module_ao/ORB_read.h"
#include "module_hamilt_pw/hamilt_pwdft/global.h"

void Gint_Gamma::cal_env(const double* wfc, double* rho, UnitCell& ucell)
{
    ModuleBase::TITLE("Grid_Integral", "cal_env");

    // it's a uniform grid to save orbital values, so the delta_r is a constant.
    const double delta_r = this->gridt->dr_uniform;
    const int max_size = this->gridt->max_atom;
    if (max_size <= 0){
        ModuleBase::WARNING_QUIT("Gint_Gamma::cal_env",
                                    "the max_size is less than 0!");
    }
    const int nbx = this->gridt->nbx;
    const int nby = this->gridt->nby;
    const int nbz = this->gridt->nbzp;
    const int ncyz = this->ny * this->nplane; // mohan add 2012-03-25
    const int bxyz = this->bxyz;

    #pragma omp parallel 
    {
        std::vector<int> block_iw(max_size, 0);
        std::vector<int> block_index(max_size+1, 0);
        std::vector<int> block_size(max_size, 0);
        std::vector<int> vindex(bxyz,0);
        #pragma omp for
        for (int grid_index = 0; grid_index < this->nbxx; grid_index++)
        {

            // get the value: how many atoms has orbital value on this grid.
            const int size = this->gridt->how_many_atoms[grid_index];
            if (size == 0)
                continue;

            // int *block_iw, *block_index, *block_size;
            ModuleBase::Array_Pool<bool> cal_flag(bxyz, size);
            Gint_Tools::get_block_info(*this->gridt,
                                       this->bxyz,
                                       size,
                                       grid_index,
                                       block_iw.data(),
                                       block_index.data(),
                                       block_size.data(),
                                       cal_flag.get_ptr_2D());
            const int LD_pool = block_index[size]; 

            // evaluate psi on grids
            ModuleBase::Array_Pool<double> psir_ylm(this->bxyz, LD_pool);
            Gint_Tools::cal_psir_ylm(*this->gridt,
                                     this->bxyz,
                                     size,
                                     grid_index,
                                     delta_r,
                                     block_index.data(),
                                     block_size.data(),
                                     cal_flag.get_ptr_2D(),
                                     psir_ylm.get_ptr_2D());

             Gint_Tools::get_vindex(this->bxyz,
                                    this->bx,
                                    this->by,
                                    this->bz,
                                    this->nplane,
                                    this->gridt->start_ind[grid_index],
                                    ncyz,
                                    vindex.data());

            for (int ia1 = 0; ia1 < size; ia1++)
            {
                const int mcell_index1 = this->gridt->bcell_start[grid_index] + ia1;
                const int iat = this->gridt->which_atom[mcell_index1];
                const int T1 = ucell.iat2it[iat];
                Atom* atom1 = &ucell.atoms[T1];
                const int I1 = ucell.iat2ia[iat];
                // get the start index of local orbitals.
                const int start1 = ucell.itiaiw2iwt(T1, I1, 0);
                for (int ib = 0; ib < this->bxyz; ib++)
                {
                    if (cal_flag[ib][ia1])
                    {
                        int iw1_lo = this->gridt->trace_lo[start1];
                        double* psi1 = &psir_ylm[ib][block_index[ia1]];
                        double tmp = 0.0;
                        for (int iw = 0; iw < atom1->nw; ++iw, ++iw1_lo)
                        {
                            tmp += psi1[iw] * wfc[iw1_lo];
                        } // iw
                        rho[vindex[ib]] += tmp;
                    } // cal_flag
                }     // ib
            }         // ia1
        }
    }
    return;
}
