#include "gint_k.h"
#include "../src_pw/global.h"
#include "global_fp.h" // mohan add 2021-01-30

#include "../module_base/ylm.h"
#include "../module_base/timer.h"

void Gint_k::cal_force_k(
	const bool isforce,
	const bool isstress,
	ModuleBase::matrix& fvl_dphi, 
	ModuleBase::matrix& svl_dphi, 
	const double *vl)
{
	ModuleBase::TITLE("Gint_k","cal_force_k");
	ModuleBase::timer::tick("Gint_k","cal_force_k");

    const int max_size = GlobalC::GridT.max_atom;

	if(max_size)
    {
        const int nbx = GlobalC::GridT.nbx;
        const int nby = GlobalC::GridT.nby;
        const int nbz_start = GlobalC::GridT.nbzp_start;
        const int nbz = GlobalC::GridT.nbzp;
        const double dv = GlobalC::ucell.omega/GlobalC::pw.ncxyz;
    
        const int ncyz = GlobalC::pw.ncy*GlobalC::pw.nczp; // mohan add 2012-03-25

        for (int i=0; i<nbx; i++)
        {
            const int ibx = i*GlobalC::pw.bx;
            for (int j=0; j<nby; j++)
            {
                const int jby = j*GlobalC::pw.by;
                for (int k=nbz_start; k<nbz_start+nbz; k++)
                {
                    const int kbz = k*GlobalC::pw.bz-GlobalC::pw.nczp_start;
    
                    const int grid_index = (k-nbz_start) + j * nbz + i * nby * nbz;
    
                    // get the value: how many atoms has orbital value on this grid.
                    const int na_grid = GlobalC::GridT.how_many_atoms[ grid_index ];
                    if(na_grid==0) continue;

                    // it's a uniform grid to save orbital values, so the delta_r is a constant.
                    const double delta_r = GlobalC::ORB.dr_uniform;						
                    
                    // here vindex refers to local potentials
                    int* vindex = Gint_Tools::get_vindex(ncyz, ibx, jby, kbz);	
                    
                    int * block_iw, * block_index, * block_size;
                    Gint_Tools::get_block_info(na_grid, grid_index, block_iw, block_index, block_size);

                    //------------------------------------------------------
                    // whether the atom-grid distance is larger than cutoff
                    //------------------------------------------------------
                    bool **cal_flag = Gint_Tools::get_cal_flag(na_grid, grid_index);

                    // set up band matrix psir_ylm and psir_DM
                    const int LD_pool = max_size*GlobalC::ucell.nwmax;

                    Gint_Tools::Array_Pool<double> psir_ylm(GlobalC::pw.bxyz, LD_pool);
                    Gint_Tools::Array_Pool<double> dpsir_ylm_x(GlobalC::pw.bxyz, LD_pool);
                    Gint_Tools::Array_Pool<double> dpsir_ylm_y(GlobalC::pw.bxyz, LD_pool);
                    Gint_Tools::Array_Pool<double> dpsir_ylm_z(GlobalC::pw.bxyz, LD_pool);

                    Gint_Tools::cal_dpsir_ylm(
                        na_grid, grid_index, delta_r,
                        block_index, block_size, 
                        cal_flag,
                        psir_ylm.ptr_2D,
                        dpsir_ylm_x.ptr_2D,
                        dpsir_ylm_y.ptr_2D,
                        dpsir_ylm_z.ptr_2D
                    );

                    double *vldr3 = Gint_Tools::get_vldr3(vl, ncyz, ibx, jby, kbz, dv);
                    const Gint_Tools::Array_Pool<double> psir_vlbr3    = Gint_Tools::get_psir_vlbr3(
						na_grid, LD_pool,
						block_index, cal_flag,
						vldr3, psir_ylm.ptr_2D);
                    
					Gint_Tools::Array_Pool<double> psir_vlbr3_DM(GlobalC::pw.bxyz, LD_pool);
                    if(isforce)
                    {
                    }
                    if(isstress)
                    {
                    }

                    free(vldr3);		vldr3=nullptr;
                    delete[] block_iw;
                    delete[] block_index;
                    delete[] block_size;

                    for(int ib=0; ib<GlobalC::pw.bxyz; ++ib)
                        free(cal_flag[ib]);
                    free(cal_flag);			cal_flag=nullptr;
                }//k
            }//j
        }//i
    }//max_size


	ModuleBase::timer::tick("Gint_k","cal_force_k");
	return;
}


