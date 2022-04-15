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
// 
// Example
//                |  ---- ixy increasing ---> |  ---- ixy increasing --->  |...
// index of sticks 0, 1, 2, ..., nst_per[0]-1, nst_per[0], ..., nst_per[1], ...
//                |___________________________|____________________________|___
// ip                           0                            1              ...
//                             npw    approximate equal to  npw   approximate equal to...
// 
// Secondly, distribute these sticks to coreors.
//Known: G, GT, GGT, ny, nx, nz, poolnproc, poolrank, ggecut
//output: ig2isz[ig], istot2bigixy[is], ixy2istot[nxy], is2ixy[is], ixy2ip[ixy], startnsz_per[ip], nst_per[ip], nst
//
void PW_Basis::distribution_method1()
{
    ModuleBase::timer::tick("PW_Basis", "distributeg_method1");

    // initial the variables needed by all proc.
    int tot_npw = 0;                     // total number of planewaves.
    this->nstot = 0;                     // total number of sticks.
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
        // for (int ix = 0; ix < nx; ++ix)
        // {
        //     for (int iy = 0; iy < ny; ++iy)
        //     {
        //         std::cout << st_length2D[ix * ny + iy] << std::setw(4);
        //     }
        //     std::cout << '\n';
        // }
        // ------------------------------------------------------------
#ifdef __MPI
        // Parallel line

        // (2) Collect the x, y indexs, length, bottom of the sticks.
        int* st_i = new int[this->nstot];                           // x or x + nx (if x < 0) of stick.
        int* st_j = new int[this->nstot];                           // y or y + ny (if y < 0) of stick.
        int* st_length = new int[this->nstot];                      // number of planewaves in stick.
        
        this->collect_st(st_length2D, st_bottom2D, st_i, st_j, st_length);

        // ------------------------------------------------------------
        std::cout << "\nThe second step done\n";
        // ------------------------------------------------------------

        // (3) Distribute sticks to cores.
        int *npw_per = new int[this->poolnproc];  // number of planewaves on each core.
        this->nstnz_per = new int[this->poolnproc]; // nz * nst(number of sticks) on each core.
        this->startnsz_per = new int[this->poolnproc];
        ModuleBase::GlobalFunc::ZEROS(npw_per, poolnproc);
        ModuleBase::GlobalFunc::ZEROS(this->nst_per, poolnproc);
        ModuleBase::GlobalFunc::ZEROS(this->nstnz_per, poolnproc);
        ModuleBase::GlobalFunc::ZEROS(startnsz_per, poolnproc);
        
        this->ixy2ip = new int[this->nxy];              // ip of core which contains stick on (x, y).
        for (int ixy = 0; ixy < this->nxy; ++ixy)
        {
            this->ixy2ip[ixy] = -1;                 // meaning this stick has not been distributed or there is no stick on (x, y).
        }
        this->divide_sticks(st_i, st_j, st_length, npw_per);
        delete[] st_length;

         // for test -----------------------------------------------------------------------------
        std::cout << "The 3-1 step done\n";
        // std::cout << "st_i    ";
        // for (int is = 0; is < this->nstot; ++is) cout << st_i[is] << setw(4) ;
        // std::cout << std::endl;
        // --------------------------------------------------------------------------------------

        this->get_istot2bigixy(st_i, st_j);
        delete[] st_i;
        delete[] st_j;
        // for test -----------------------------------------------------------------------------
        std::cout << "The 3-2 step done\n";
        // --------------------------------------------------------------------------------------

        // (4) Send npw_per, nst[poolrank], st_* to all cores.
        this->npw = npw_per[0];
        this->nst = nst_per[0];

        for (int ip = 1; ip < this->poolnproc; ++ip)
        {
            MPI_Send(&npw_per[ip], 1, MPI_INT, ip, 0, POOL_WORLD);
            MPI_Send(&nst_per[ip], 1, MPI_INT, ip, 0, POOL_WORLD);
        }
        delete[] npw_per;

#else
        // Serial line
        this->npw = tot_npw;
        this->nst = this->nstot;

        this->nstnz_per = new int[1];
        this->nstnz_per[0] = this->nst * this->nz;
        this->startnsz_per = new int[1];
        this->startnsz_per[0] = 0;

        this->ixy2istot = new int[nxy];
        this->istot2bigixy = new int[this->nstot];
        this->ixy2ip = new int[nxy];              // ip of core which contains stick on (x, y).
        int st_move = 0;
        for (int ixy = 0; ixy < nxy; ++ixy)
        {
            if (st_length2D[ixy] > 0)
            {
                this->ixy2istot[ixy] = st_move;
                this->istot2bigixy[st_move] = ixy / ny * bigny + ixy % ny;
                this->ixy2ip[ixy] = 0;
                st_move++;
            }
            else
            {
            ixy2istot[ixy] = -1;
            ixy2ip[ixy] = -1;
            }
        }
#endif
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
    MPI_Bcast(this->nst_per, this->poolnproc, MPI_INT, 0 , POOL_WORLD);

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
// (2) Collect the x, y indexs, length of the sticks.
// Firstly, we scan the area and construct temp_st_*.
// Then, as we will distribute the longest sticks preferentially in Step(3),
// we will sort temp_st_length from largest to smallest, and reaarange st_* to the same order.
// known: tot_npw, this->nstot, st_length2D, st_bottom2D
// output: st_i, st_j, st_length
//
void PW_Basis::collect_st(
    int* st_length2D,                               // the number of planewaves that belong to the stick located on (x, y), stored in 2d x-y plane.
    int* st_bottom2D,                               // the z-coordinate of the bottom of stick on (x, y), stored in 2d x-y plane.
    int* st_i,                                      // x or x + nx (if x < 0) of stick.
    int* st_j,                                      // y or y + ny (if y < 0) of stick.
    int* st_length                                 // number of planewaves in stick, stored in 1d array with this->nstot elements.
)
{
    int *temp_st_i = new int[this->nstot];                      // x or x + nx (if x < 0) of stick.
    int *temp_st_j = new int[this->nstot];                      // y or y + ny (if y < 0) of stick.
    double *temp_st_length = new double[this->nstot];           // length of sticks.
    ModuleBase::GlobalFunc::ZEROS(temp_st_length, this->nstot);

    int ibox[3] = {0, 0, 0};                            // an auxiliary vector, determine the boundary of the scanning area.
    ibox[0] = int(this->nx / 2) + 1;                    // scan x from -ibox[0] to ibox[0].
    ibox[1] = int(this->ny / 2) + 1;                    // scan y from -ibox[1] to ibox[1].
    ibox[2] = int(this->nz / 2) + 1;                    // scan z from -ibox[2] to ibox[2].

    ModuleBase::Vector3<double> f;
    int is = 0; // index of stick.

    int iy_start = -ibox[1]; // determine the scaning area along x-direct, if gamma-only, only positive axis is used.
    int iy_end = ibox[1];
    if (this->gamma_only)
    {
        iy_start = 0;
        iy_end = this->ny - 1;
    }
    for (int ix = -ibox[0]; ix <= ibox[0]; ++ix)
    {
        for (int iy = iy_start; iy <= iy_end; ++iy)
        {
            // we have shifted all sticks to the first quadrant in x-y plane before.
            // (ix, iy, iz) is the direct coordinates of planewaves.
            // x and y is the coordinates of shifted sticks in x-y plane.
            // for example, if nx = ny = 10, we will shift the stick on (-1, 2) to (9, 2),
            // so that its index in st_length and st_bottom is 9 * 10 + 2 = 92.
            int x = ix;
            int y = iy;
            if (x < 0) x += nx;
            if (y < 0) y += ny;
            int index = x * this->ny + y;
            if (st_length2D[index] > 0) // meaning there is a stick on (x, y) point.
            {
                bool find_stick = false;
                for (int iz = st_bottom2D[index]; iz < st_bottom2D[index] + st_length2D[index]; ++iz)
                {
                    f.x = ix;
                    f.y = iy;
                    f.z = iz;
                    double modulus = f * (GGT * f);
                    if (modulus <= ggecut)
                    {
                        find_stick = true;
                        break;
                    }                  
                }
                if (find_stick)
                {
                    temp_st_i[is] = x;
                    temp_st_j[is] = y;
                    temp_st_length[is] = static_cast<double>(st_length2D[index]);
                    ++is;
                }   
            }
        }
    }
    assert(is == this->nstot);
    std::cout<<"collect sticks done\n";

    // As we will distribute the longest sticks preferentially in Step(3), we rearrange st_* in the order of length decreasing.

    int *st_sorted_index = new int[this->nstot]; // indexs in the order of length increasing.
    st_sorted_index[0] = 0;
    ModuleBase::heapsort(this->nstot, temp_st_length, st_sorted_index); // sort st_* in the order of length decreasing.

    for (int istot = 0; istot < this->nstot; ++istot)
    {
        st_length[istot] = static_cast<int>(temp_st_length[istot]);
        st_i[istot] = temp_st_i[st_sorted_index[istot]];
        st_j[istot] = temp_st_j[st_sorted_index[istot]];
    }
    // std::cout << "st_length    ";
    // for (int is = 0; is < this->nstot; ++is) std::cout << st_length[is] << std::setw(4);
    // std::cout << "\n";

    delete[] temp_st_i;
    delete[] temp_st_j;
    delete[] temp_st_length;
    delete[] st_sorted_index;
    return;
}

//
// (3-1) Distribute sticks to cores.
// We have rearranged sticks in the order of length decreasing, so that we will distribute the longest stick preferentially here.
// For each stick, we find the core that contains the least planewaves firstly, and distribute the stick to it,
// then update npw_per, this->nstnz_per, this->ixy2ip, and this->startnsz_per.
// known: tot_npw, this->nstot, st_i, st_j, st_length
// output: npw_per, nst_per, this->nstnz_per, this->ixy2ip, this->startnsz_per
//
void PW_Basis::divide_sticks(
    int* st_i,          // x or x + nx (if x < 0) of stick.
    int* st_j,          // y or y + ny (if y < 0) of stick.
    int* st_length,     // the stick on (x, y) consists of st_length[x*ny+y] planewaves.
    int* npw_per       // number of planewaves on each core.
)
{
    int ipmin = 0; // The ip of core containing least number of planewaves.
    for (int is = 0; is < this->nstot; ++is)
    {
        // find the ip of core containing the least planewaves.
        for (int ip = 0; ip < this->poolnproc; ++ip)
        {
            //const int non_zero_grid = nst_per[ip] * this->nz;       // number of reciprocal planewaves on this core.
            const int npwmin = npw_per[ipmin];
            const int npw_ip = npw_per[ip];
            const int nstmin = nst_per[ipmin];
            const int nst_ip = nst_per[ip];

            if (npw_ip == 0)
            {
                // if (non_zero_grid + nz < this->nrxx) // assert reciprocal planewaves is less than real space planewaves.
                // {
                    ipmin = ip;
                    break;
                // }
            }
            else if (npw_ip < npwmin)
            {
                // if (non_zero_grid + nz < this->nrxx) // assert reciprocal planewaves is less than real space planewaves.
                // {
                    ipmin = ip;
                // }
            }
            else if (npw_ip == npwmin && nst_ip < nstmin)
            {
                // if (non_zero_grid + nz < this->nrxx)
                // {
                    ipmin = ip;
                // }
            }
        }
        nst_per[ipmin]++;
        this->nstnz_per[ipmin] += this->nz;
        npw_per[ipmin] += st_length[is];
        this->ixy2ip[st_i[is] * this->ny + st_j[is]] = ipmin;
    }
    // for test --------------------------------------------------------------------------
    // for(int i = 0; i < this->poolnproc; ++i) std::cout<<this->nstnz_per[i]<<std::setw(4);
    // std::cout<<"\n";
    // for(int i = 0; i < this->poolnproc; ++i) std::cout<<nst_per[i]<<std::setw(4);
    // std::cout<<"\n";
    // -----------------------------------------------------------------------------------
    for (int ip = 1; ip < poolnproc; ++ip)
    {
        this->startnsz_per[ip] = this->startnsz_per[ip - 1] + this->nstnz_per[ip - 1];
    }
    return;
}

//
// (3-2) Rearrange sticks in the order of the ip of core increasing, in each core, sticks are sorted in the order of ixy increasing.
// (st_start + st_move) is the new index of sticks.
// Then get istot2bigixy (istot2bigixy[is]: iy + ix * ny of is^th stick among all sticks) on the first core
// and ixy2istot (ixy2istot[iy + ix * ny]: is of stick on (iy, ix) among all sticks).
// known: this->nstot, st_i, st_j, this->startnsz_per
// output: istot2bigixy, ixy2istot
//
void PW_Basis::get_istot2bigixy(
    int* st_i,          // x or x + nx (if x < 0) of stick.
    int* st_j          // y or y + ny (if y < 0) of stick.
)
{
    assert(this->poolrank == 0);
    this->ixy2istot = new int[this->nx * this->ny];
    this->istot2bigixy = new int[this->nstot];
    int* st_move = new int[this->poolnproc]; // st_move[ip]: this is the st_move^th stick on ip^th core.
    for (int ixy = 0; ixy < this->nx * this->ny; ++ixy)
    {
        this->ixy2istot[ixy] = -1;
    }
    ModuleBase::GlobalFunc::ZEROS(this->istot2bigixy, this->nstot);
    ModuleBase::GlobalFunc::ZEROS(st_move, this->poolnproc);

    for (int ixy = 0; ixy < this->nxy; ++ixy)
    {
        int ip = this->ixy2ip[ixy];
        if (ip != -1)
        {        
            this->istot2bigixy[this->startnsz_per[ip] / this->nz + st_move[ip]] = (ixy / ny)*bigny + ixy % ny;
            this->ixy2istot[ixy] = this->startnsz_per[ip] / this->nz + st_move[ip];
            st_move[ip]++;
        }
    }
    delete[] st_move;
    return;
}

}
