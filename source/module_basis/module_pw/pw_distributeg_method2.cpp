#include "pw_basis.h"
#include "module_base/mymath.h"
#include "module_base/global_function.h"


namespace ModulePW
{
///
/// Distribute planewaves in reciprocal space to cores.
/// Firstly, divide the sphere in reciprocal space into sticks, which are vertical to x-y plane.
/// Secondly, distribute these sticks to coreors.
/// 
/// Example
///                                ---- ixy increasing --->
/// index of sticks 0, 1, 2, ..., nst_per[0]-1, nst_per[0], ..., nst_per[1]-1, ...
///                |___________________________|______________________________|___
/// ip                           0                            1              ...
///                                 nst_per[i] approxiamte equal to nst_per[j]
/// Known: G, GT, GGT, fftny, fftnx, nz, poolnproc, poolrank, ggecut
/// output: ig2isz[ig], istot2ixy[is], is2fftixy[is], fftixy2ip[ixy], startnsz_per[ip], nst_per[ip], nst
///
void PW_Basis::distribution_method2()
{

    // initial the variables needed by all proc.
    int *st_bottom2D = new int[fftnxy];             // st_bottom2D[ixy], minimum z of stick on (x, y).
    int *st_length2D = new int[fftnxy];             // st_length2D[ixy], number of planewaves in stick on (x, y).
    delete[] this->nst_per; this->nst_per = new int[this->poolnproc]; // number of sticks on each core.
    delete[] this->npw_per;   this->npw_per = new int[this->poolnproc];  // number of planewaves on each core.
    delete[] this->fftixy2ip; this->fftixy2ip = new int[this->fftnxy];              // ip of core which contains the stick on (x, y).
    for (int ixy = 0; ixy < this->fftnxy; ++ixy)
        this->fftixy2ip[ixy] = -1;                 // meaning this stick has not been distributed or there is no stick on (x, y).
    if (poolrank == 0)
    {
        // (1) Count the total number of planewaves (tot_npw) and sticks (this->nstot).                  
        
        // Actually we will scan [(2 * ibox[0] + 1) * (2 * ibox[1] + 1)] points on x-y plane,
        // but we define st_length2D with (fftny * fftnx) points here, because we assume that the diameter
        // of the sphere is shorter than the sides of the cube.
        // calculate this->nstot and this->npwtot, liy, riy
        this->count_pw_st(st_length2D, st_bottom2D);
    }
#ifdef __MPI
    MPI_Bcast(&this->npwtot, 1, MPI_INT, 0, this->pool_world);
    MPI_Bcast(&this->nstot, 1, MPI_INT, 0, this->pool_world);
    MPI_Bcast(&liy, 1, MPI_INT, 0, this->pool_world);
    MPI_Bcast(&riy, 1, MPI_INT, 0, this->pool_world);
    MPI_Bcast(&lix, 1, MPI_INT, 0, this->pool_world);
    MPI_Bcast(&rix, 1, MPI_INT, 0, this->pool_world);
#endif
    delete[] this->istot2ixy; this->istot2ixy = new int[this->nstot];

    if(poolrank == 0)
    {
#ifdef __MPI
        // Parallel line
        // (2) Devide the sticks to each core, sticks are in the order of ixy increasing.
        // get nst_per and startnsz_per
        this->startnsz_per = new int[this->poolnproc];
        this->divide_sticks_2();

        // (3) Create the maps from ixy to ip, istot, and from istot to ixy
        // get istot2ixy, fftixy2ip, npw_per
        this->create_maps(st_length2D);
        //We do not need startnsz_per after it.
        delete[] this->startnsz_per;
        this->startnsz_per=nullptr;
#else
        // Serial line
        // get nst_per, npw_per, fftixy2ip, and istot2ixy
        this->nst_per[0] = this->nstot;
        this->npw_per[0] = this->npwtot;
        int st_move = 0;
        for (int ixy = 0; ixy < fftnxy; ++ixy)
        {
            if (st_length2D[ixy] > 0)
            {
                this->istot2ixy[st_move] = ixy / fftny * ny + ixy % fftny;
                this->fftixy2ip[ixy] = 0;
                st_move++;
            }
        }
#endif
    }
#ifdef __MPI
    MPI_Bcast(st_length2D, this->fftnxy, MPI_INT, 0, this->pool_world);
    MPI_Bcast(st_bottom2D, this->fftnxy, MPI_INT, 0, this->pool_world);
    MPI_Bcast(this->fftixy2ip, this->fftnxy, MPI_INT, 0, this->pool_world);
    MPI_Bcast(this->istot2ixy, this->nstot, MPI_INT, 0, this->pool_world);
    MPI_Bcast(this->nst_per, this->poolnproc, MPI_INT, 0 , this->pool_world);
    MPI_Bcast(this->npw_per, this->poolnproc, MPI_INT, 0 , this->pool_world);
#endif
    this->npw = this->npw_per[this->poolrank];
    this->nst = this->nst_per[this->poolrank];
    this->nstnz = this->nst * this->nz;

    // (5) Construct ig2isz and is2fftixy. 
    this->get_ig2isz_is2fftixy(st_bottom2D, st_length2D);

    delete[] st_bottom2D;
    delete[] st_length2D;
    return;
}

/// 
/// (2) Devide the sticks to each core according to the number of sticks
/// Sticks are in the order of ixy increasing.
/// known: this->nstot, this->poolnproc
/// output: nst_per, this->startnsz_per
/// 
void PW_Basis::divide_sticks_2()
{
    ModuleBase::GlobalFunc::ZEROS(nst_per, this->poolnproc);
    
    int average_nst = this->nstot / this->poolnproc;
    int mods = this->nstot % this->poolnproc;
    
    this->startnsz_per[0] = 0;
    for (int ip = 0; ip < this->poolnproc; ++ip)
    {
        nst_per[ip] = average_nst;
        if (ip < mods) nst_per[ip]++;
        if (ip >= 1) this->startnsz_per[ip] = this->startnsz_per[ip-1] + this->nst_per[ip-1] * this->nz; 
    }
}

// 
// (3) Create the maps from ixy to ip, istot, and from istot to ixy, and construt npw_per simultaneously.
// known: st_length2D
// output: this->fftixy2ip, this->istot2ixy, npw_per
// 
void PW_Basis::create_maps(
    int* st_length2D  // the number of planewaves that belong to the stick located on (x, y), stored in 2d x-y plane.
)
{
    ModuleBase::GlobalFunc::ZEROS(this->istot2ixy, this->nstot);
    ModuleBase::GlobalFunc::ZEROS(this->npw_per, poolnproc);
    int ip = 0;
    int st_move = 0; // the number of sticks that have been found.
    for (int ixy = 0; ixy < this->fftnxy; ++ixy)
    {
        if (st_length2D[ixy] > 0)
        {
            this->istot2ixy[st_move] = ixy / fftny * ny + ixy % fftny;
            this->fftixy2ip[ixy] = ip;
            this->npw_per[ip] += st_length2D[ixy];
            st_move++;
            if (ip < this->poolnproc - 1)
            {
                // all of sticks on current core are found, skip to next core
                if (st_move * this->nz >= this->startnsz_per[ip + 1]) ip++; 
            }
        }
    }
}
}