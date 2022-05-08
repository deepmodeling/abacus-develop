#include "gint_gamma.h"
#include "gint_tools.h"
#include "grid_technique.h"
#include "../module_orbital/ORB_read.h"
#include "../src_pw/global.h"
#include "../module_base/blas_connector.h"
#include "../src_parallel/parallel_reduce.h"
#include "../module_base/timer.h"

#include "global_fp.h" // mohan add 2021-01-30

void Gint_Gamma::cal_force_new(double** DM_in, const double*const vlocal, ModuleBase::matrix& force)
{
    ModuleBase::TITLE("Grid_Integral","cal_force_new");
    ModuleBase::timer::tick("Gint_Gamma","cal_force_new");
    this->save_atoms_on_grid(GlobalC::GridT);
    this->gamma_force_new(DM_in, vlocal, force);

    ModuleBase::timer::tick("Gint_Gamma","cal_force_new");
}

void Gint_Gamma::gamma_force_new(const double*const*const DM, const double*const vlocal, ModuleBase::matrix& force)
{
    ModuleBase::TITLE("Grid_Integral","gamma_force_new");
    ModuleBase::timer::tick("Gint_Gamma","gamma_force_new");

	if(max_size)
    {
        const int nbx = GlobalC::GridT.nbx;
        const int nby = GlobalC::GridT.nby;
        const int nbz_start = GlobalC::GridT.nbzp_start;
        const int nbz = GlobalC::GridT.nbzp;
    
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
                        na_grid, LD_pool, grid_index, delta_r,
                        block_index, block_size, 
                        cal_flag,
                        psir_ylm.ptr_2D, dpsir_ylm_x.ptr_2D, dpsir_ylm_y.ptr_2D, dpsir_ylm_z.ptr_2D
                    );
                    double *vldr3 = this->get_vldr3(vlocal, ncyz, ibx, jby, kbz);
                    const Gint_Tools::Array_Pool<double> psir_vlbr3 = Gint_Tools::get_psir_vlbr3(na_grid, LD_pool, block_index, cal_flag, vldr3, psir_ylm.ptr_2D);
                    this-> cal_meshball_force(grid_index, na_grid, LD_pool, block_iw, block_size, block_index,
                        cal_flag, psir_vlbr3.ptr_2D, dpsir_ylm_x.ptr_2D, dpsir_ylm_y.ptr_2D, dpsir_ylm_z.ptr_2D, DM, force);
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

    ModuleBase::timer::tick("Gint_Gamma","gamma_force_new");

}

void Gint_Gamma::cal_meshball_force(
    const int grid_index,
    const int na_grid,  					    // how many atoms on this (i,j,k) grid
	const int LD_pool,
	const int*const block_iw,				    // block_iw[na_grid],	index of wave functions for each block
	const int*const block_size, 			    // block_size[na_grid],	number of columns of a band
	const int*const block_index,		    	// block_index[na_grid+1], count total number of atomis orbitals
	const bool*const*const cal_flag,	    	// cal_flag[GlobalC::pw.bxyz][na_grid],	whether the atom-grid distance is larger than cutoff
	const double*const*const psir_vlbr3,	    // psir_vlbr3[GlobalC::pw.bxyz][LD_pool]
    const double*const*const dpsir_x,	    // psir_vlbr3[GlobalC::pw.bxyz][LD_pool]
    const double*const*const dpsir_y,	    // psir_vlbr3[GlobalC::pw.bxyz][LD_pool]
    const double*const*const dpsir_z,	    // psir_vlbr3[GlobalC::pw.bxyz][LD_pool]
    const double*const*const DM,
    ModuleBase::matrix &force)
{
    constexpr char side='L', uplo='U';
    constexpr char transa='N', transb='N';
    constexpr double alpha_symm=1, alpha_gemm=1, beta=1;    
    constexpr int inc=1;

    Gint_Tools::Array_Pool<double> psir_vlbr3_DM(GlobalC::pw.bxyz, LD_pool);
    ModuleBase::GlobalFunc::ZEROS(psir_vlbr3_DM.ptr_1D, GlobalC::pw.bxyz*LD_pool);

    for (int ia1=0; ia1<na_grid; ia1++)
    {
        const int iw1_lo=block_iw[ia1];
        for (int ia2=0; ia2<na_grid; ia2++)
        {
            int first_ib=0, last_ib=0;
            for(int ib=0; ib<GlobalC::pw.bxyz; ++ib)
            {
                if(cal_flag[ib][ia1] && cal_flag[ib][ia2])
                {
                    first_ib=ib;
                    break;
                }
            }
            for(int ib=GlobalC::pw.bxyz-1; ib>=0; --ib)
            {
                if(cal_flag[ib][ia1] && cal_flag[ib][ia2])
                {
                    last_ib=ib+1;
                    break;
                }
            }
            const int ib_length=last_ib-first_ib;
            if(ib_length<=0) continue;

            int cal_pair_num=0;
            for(int ib=first_ib; ib<last_ib; ++ib)
            {
                cal_pair_num += cal_flag[ib][ia1] && cal_flag[ib][ia2];
            }
            const int iw2_lo=block_iw[ia2];
            if(cal_pair_num>ib_length/4)
            {
                dgemm_(&transa, &transb, &block_size[ia2], &ib_length, &block_size[ia1], 
                    &alpha_gemm, &DM[iw1_lo][iw2_lo], &GlobalC::GridT.lgd, 
                    &psir_vlbr3[first_ib][block_index[ia1]], &LD_pool, 
                    &beta, &psir_vlbr3_DM.ptr_2D[first_ib][block_index[ia2]], &LD_pool);
            }
            else
            {
                for(int ib=first_ib; ib<last_ib; ++ib)
                {
                    if(cal_flag[ib][ia1] && cal_flag[ib][ia2])
                    {
                        dgemv_(&transa, &block_size[ia2], &block_size[ia1], 
                            &alpha_gemm, &DM[iw1_lo][iw2_lo], &GlobalC::GridT.lgd,
                            &psir_vlbr3[ib][block_index[ia1]], &inc,
                            &beta, &psir_vlbr3_DM.ptr_2D[ib][block_index[ia2]], &inc);
                    }
                }
            }
        }// ia2       
    } // ia1  

    for(int ia1=0;ia1<na_grid;ia1++)
    {
        const int mcell_index=GlobalC::GridT.bcell_start[grid_index] + ia1;
        const int iat=GlobalC::GridT.which_atom[mcell_index]; // index of atom

        for(int ib=0;ib<GlobalC::pw.bxyz;ib++)
        {
            const double rx = ddot_(&block_size[ia1], &psir_vlbr3_DM.ptr_2D[ib][block_index[ia1]], &inc, &dpsir_x[ib][block_index[ia1]], &inc);
            force(iat,0)+=rx*2.0;
            const double ry = ddot_(&block_size[ia1], &psir_vlbr3_DM.ptr_2D[ib][block_index[ia1]], &inc, &dpsir_y[ib][block_index[ia1]], &inc);
            force(iat,1)+=ry*2.0;
            const double rz = ddot_(&block_size[ia1], &psir_vlbr3_DM.ptr_2D[ib][block_index[ia1]], &inc, &dpsir_z[ib][block_index[ia1]], &inc);
            force(iat,2)+=rz*2.0;
        }
    }          

    return;
}