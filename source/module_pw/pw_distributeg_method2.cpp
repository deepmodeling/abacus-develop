#include "pw_basis.h"
#include "../module_base/mymath.h"
#include "../src_parallel/parallel_global.h"
#include "../module_base/global_function.h"
#include "iostream"
#include "../module_base/timer.h"


namespace ModulePW
{
//
// Distribute planewaves in reciprocal space to coreors.
// Firstly, divide the sphere in reciprocal space into sticks, which are vertical to x-y plane.
// Secondly, distribute these sticks to coreors.
// 
// Example
//                                ---- ixy increasing --->
// index of sticks 0, 1, 2, ..., nst_per[0]-1, nst_per[0], ..., nst_per[1], ...
//                |___________________________|____________________________|___
// ip                           0                            1              ...
// 
//Known: G, GT, GGT, ny, nx, nz, poolnproc, poolrank, ggecut
//output: ig2isz[ig], istot2bigixy[is], ixy2istot[nxy], is2ixy[is], ixy2ip[ixy], startnsz_per[ip], nst_per[ip], nst
//
void PW_Basis::distribution_method2()
{
    ModuleBase::timer::tick("PW_Basis", "distributeg_method2");

    // initial the variables needed by all proc.
    int tot_npw = 0;                     // total number of planewaves.
    this->nstot = 0;                     // total number of sticks.
   // int st_start = 0;                    // index of the first stick on current proc.
    int *st_bottom2D = NULL;             // st_bottom2D[ixy], minimum z of stick on (x, y).
    int *st_length2D = NULL;             // st_length2D[ixy], number of planewaves in stick on (x, y).

    this->nst_per = new int[this->poolnproc]; // number of sticks on each core.
    if (poolrank == 0)
    {
        // (1) Count the total number of planewaves (tot_npw) and sticks (this->nstot).                  
        
        // Actually we will scan [(2 * ibox[0] + 1) * (2 * ibox[1] + 1)] points on x-y plane,
        // but we define st_length2D with (ny * nx) points here, because we assume that the diameter
        // of the sphere is shorter than the sides of the cube.
        st_length2D = new int[nxy];                    // the number of planewaves that belong to the stick located on (x, y).
        st_bottom2D = new int[nxy];                    // the z-coordinate of the bottom of stick on (x, y).
        ModuleBase::GlobalFunc::ZEROS(st_length2D, this->nxy);
        ModuleBase::GlobalFunc::ZEROS(st_bottom2D, this->nxy);

        this->count_pw_st(tot_npw, st_length2D, st_bottom2D);
        // for test --------------------------------------------------
        std::cout << "The first step done\n";
        std::cout << "tot_npw   " << tot_npw << '\n';
        std::cout << "this->nstot   " << this->nstot << '\n';
        for (int ix = 0; ix < nx; ++ix)
        {
            for (int iy = 0; iy < ny; ++iy)
            {
                std::cout << st_length2D[ix * ny + iy] << std::setw(4);
            }
            std::cout << '\n';
        }
        // ------------------------------------------------------------

        // (2) Devide the sticks to each core, sticks are in the order of ixy increasing.

        ModuleBase::GlobalFunc::ZEROS(nst_per, this->poolnproc);
        this->divide_sticks2();
        // for test ----------------------------------------------------------------------------
        std::cout << "nst_per  ";
        for (int ip = 0; ip < this->poolnproc; ++ip) std::cout << nst_per[ip] << std::setw(4);
        std::cout << "\n";
        //-------------------------------------------------------------------------------------- 

        // (3) Create the maps from ixy to ip, istot, and from istot to ixy

        int *npw_per = new int[this->poolnproc];  // number of planewaves on each core.
        ModuleBase::GlobalFunc::ZEROS(npw_per, this->poolnproc);
        this->create_maps(st_length2D, npw_per);
        // for test ----------------------------------------------------------------------------
        std::cout << "npw_per  ";
        for (int ip = 0; ip < this->poolnproc; ++ip) std::cout << npw_per[ip] << std::setw(4);
        std::cout << "\n";
        //-------------------------------------------------------------------------------------- 

        // (4) Send npw_per, nst[poolrank], st_* to all cores.
        this->npw = npw_per[0];
        this->nst = nst_per[0];

#ifdef __MPI
        for (int ip = 1; ip < this->poolnproc; ++ip)
        {
            MPI_Send(&npw_per[ip], 1, MPI_INT, ip, 0, POOL_WORLD);
            MPI_Send(&nst_per[ip], 1, MPI_INT, ip, 0, POOL_WORLD);
        }
#endif
        delete[] npw_per;
    }
    else
    {
#ifdef __MPI
        MPI_Status ierror;
        MPI_Recv(&npw, 1, MPI_INT, 0, 0, POOL_WORLD, &ierror);  // number of planewaves in current proc.
        MPI_Recv(&nst, 1, MPI_INT, 0, 0, POOL_WORLD, &ierror);
        std::cout << this->poolrank << " recive done.\n";
#endif
    }
#ifdef __MPI
    MPI_Bcast(&tot_npw, 1, MPI_INT, 0, POOL_WORLD);
    MPI_Bcast(&this->nstot, 1, MPI_INT, 0, POOL_WORLD);
    MPI_Bcast(&liy, 1, MPI_INT, 0, POOL_WORLD);
    MPI_Bcast(&riy, 1, MPI_INT, 0, POOL_WORLD);
    MPI_Bcast(this->nst_per, this->poolnproc, MPI_INT, 0 , POOL_WORLD);
    if (this->poolrank != 0)
    {
        st_bottom2D = new int[this->nxy];                      // minimum z of stick.
        st_length2D = new int[this->nxy];                      // number of planewaves in stick.
        this->ixy2ip = new int[this->nxy];              // ip of core which contains stick on (x, y).
        this->ixy2istot = new int[this->nxy];
        this->istot2bigixy = new int[this->nstot];
    }

    MPI_Bcast(st_length2D, this->nxy, MPI_INT, 0, POOL_WORLD);
    MPI_Bcast(st_bottom2D, this->nxy, MPI_INT, 0, POOL_WORLD);
    MPI_Bcast(this->ixy2ip, this->nxy, MPI_INT, 0, POOL_WORLD);
    MPI_Bcast(this->istot2bigixy, this->nstot, MPI_INT, 0, POOL_WORLD);
    MPI_Bcast(this->ixy2istot, this->nxy, MPI_INT, 0, POOL_WORLD);

    std::cout << "Bcast done\n";
#endif
    this->nstnz = this->nst * this->nz;

    // (5) Construct ig2isz and is2ixy. 
    this->get_ig2isz_is2ixy(st_bottom2D, st_length2D);

    if (st_bottom2D != NULL) delete[] st_bottom2D;
    if (st_length2D != NULL) delete[] st_length2D;
    // for test ----------------------------------------------
    if (poolrank==0) std::cout << "The fifth step done\n";
    // -------------------------------------------------------
    ModuleBase::timer::tick("PW_Basis", "distributeg_method1");
    return;
}

// 
// (2) Devide the sticks to each core, sticks are in the order of ixy increasing.
// known: this->nstot, this->poolnproc
// output: nst_per, this->nstnz_per, this->startnsz_per
// 
void PW_Basis::divide_sticks2()
{
    this->nstnz_per = new int[this->poolnproc]; // nz * nst(number of sticks) on each core.
    this->startnsz_per = new int[this->poolnproc];
    ModuleBase::GlobalFunc::ZEROS(this->nstnz_per, this->poolnproc);
    ModuleBase::GlobalFunc::ZEROS(this->startnsz_per, this->poolnproc);
    
    int average_nst = this->nstot / this->poolnproc;
    int mods = this->nstot % this->poolnproc;
    for (int ip = 0; ip < this->poolnproc; ++ip)
    {
        nst_per[ip] = average_nst;
        if (ip < mods) nst_per[ip]++;
        this->nstnz_per[ip] = nst_per[ip] * this->nz;
        if (ip >= 1) this->startnsz_per[ip] = this->startnsz_per[ip-1] + this->nstnz_per[ip-1]; 
    }
}

// 
// (3) Create the maps from ixy to ip, istot, and from istot to ixy, and construt npw_per simultaneously.
// known: st_length2D
// output: this->ixy2ip, this->ixy2istot, this->istot2bigixy, npw_per
// 
void PW_Basis::create_maps(
    int* st_length2D,  // the number of planewaves that belong to the stick located on (x, y), stored in 2d x-y plane.
    int* npw_per       // number of planewaves on each core.
)
{
    this->ixy2ip = new int[this->nxy];
    this->ixy2istot = new int[this->nxy];
    this->istot2bigixy = new int[this->nstot];

    ModuleBase::GlobalFunc::ZEROS(this->istot2bigixy, this->nstot);
    int ip = 0;
    int st_move = 0; // the number of sticks that have been found.
    for (int ixy = 0; ixy < this->nxy; ++ixy)
    {
        if (st_length2D[ixy] > 0)
        {
            this->ixy2istot[ixy] = st_move;
            this->istot2bigixy[st_move] = ixy / ny * bigny + ixy % ny;
            this->ixy2ip[ixy] = ip;
            npw_per[ip] += st_length2D[ixy];
            st_move++;
            if (ip < this->poolnproc - 1)
            {
                // all of sticks on current core are found, skip to next core
                if (st_move * this->nz >= this->startnsz_per[ip + 1]) ip++; 
            }
        }
        else
        {
            this->ixy2istot[ixy] = -1;
            this->ixy2ip[ixy] = -1;
        }
    }
}
}