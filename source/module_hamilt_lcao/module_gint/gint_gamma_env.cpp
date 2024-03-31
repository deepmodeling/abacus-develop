#include "gint_gamma.h"
#include "grid_technique.h"
#include "module_base/timer.h"
#include "module_base/ylm.h"
#include "module_basis/module_ao/ORB_read.h"
#include "module_hamilt_pw/hamilt_pwdft/global.h"

void Gint_Gamma::cal_env(const double* wfc, double* rho)
{
    ModuleBase::TITLE("Grid_Integral", "cal_env");

    // it's a uniform grid to save orbital values, so the delta_r is a constant.
    const double delta_r = GlobalC::ORB.dr_uniform;
    const int max_size = this->gridt->max_atom;
    const int LD_pool = max_size * GlobalC::ucell.nwmax;

    if (max_size != 0)
    {
        const int nbx = this->gridt->nbx;
        const int nby = this->gridt->nby;
        const int nbz_start = this->gridt->nbzp_start;
        const int nbz = this->gridt->nbzp;
        const int ncyz = this->ny * this->nplane; // mohan add 2012-03-25

        for (int grid_index = 0; grid_index < this->nbxx; grid_index++)
        {

            // get the value: how many atoms has orbital value on this grid.
            const int size = this->gridt->how_many_atoms[grid_index];
            if (size == 0)
                continue;

            int *block_iw, *block_index, *block_size;
            bool** cal_flag;
            Gint_Tools::get_block_info(*this->gridt,
                                       this->bxyz,
                                       size,
                                       grid_index,
                                       block_iw,
                                       block_index,
                                       block_size,
                                       cal_flag);

            // evaluate psi on grids
            Gint_Tools::Array_Pool<double> psir_ylm(this->bxyz, LD_pool);
            Gint_Tools::cal_psir_ylm(*this->gridt,
                                     this->bxyz,
                                     size,
                                     grid_index,
                                     delta_r,
                                     block_index,
                                     block_size,
                                     cal_flag,
                                     psir_ylm.ptr_2D);

            int* vindex
                = Gint_Tools::get_vindex(this->bxyz,
                                         this->bx,
                                         this->by,
                                         this->bz,
                                         this->nplane,
                                         this->gridt->start_ind[grid_index],
                                         ncyz);

            for (int ia1 = 0; ia1 < size; ia1++)
            {
                const int mcell_index1
                    = this->gridt->bcell_start[grid_index] + ia1;
                const int iat = this->gridt->which_atom[mcell_index1];
                const int T1 = GlobalC::ucell.iat2it[iat];
                Atom* atom1 = &GlobalC::ucell.atoms[T1];
                const int I1 = GlobalC::ucell.iat2ia[iat];
                // get the start index of local orbitals.
                const int start1 = GlobalC::ucell.itiaiw2iwt(T1, I1, 0);
                for (int ib = 0; ib < this->bxyz; ib++)
                {
                    if (cal_flag[ib][ia1])
                    {
                        int iw1_lo = this->gridt->trace_lo[start1];
                        double* psi1 = &psir_ylm.ptr_2D[ib][block_index[ia1]];
                        double tmp = 0.0;
                        for (int iw = 0; iw < atom1->nw; ++iw, ++iw1_lo)
                        {
                            tmp += psi1[iw] * wfc[iw1_lo];
                        } // iw
                        rho[vindex[ib]] += tmp;
                    } // cal_flag
                }     // ib
            }         // ia1

            delete[] vindex;
            delete[] block_iw;
            delete[] block_index;
            delete[] block_size;
            for (int ib = 0; ib < this->bxyz; ++ib)
            {
                delete[] cal_flag[ib];
            }
            delete[] cal_flag;
        }
    }

    return;
}
