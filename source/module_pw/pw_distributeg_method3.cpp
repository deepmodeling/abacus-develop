#include "pw_basis.h"
#include "../module_base/mymath.h"
#include "../src_parallel/parallel_global.h"
#include "../module_base/global_function.h"
#include "iostream"



namespace ModulePW
{
//
// Distribute planewaves in reciprocal space to coreors.
// Firstly, divide the sphere in reciprocal space into sticks, which are vertical to x-y plane.
// Secondly, distribute these sticks to coreors.
//Known: G, GT, GGT, nx, ny, nz, poolnproc, poolrank, ggecut
//output: ig2isz[ig], istot2ixy[is], ixy2istot[this->nxy], is2ixy[is], ixy2ip[ixy], startnsz_per[ip], nst_per[ip], gg[ig], gcar[ig], gdirect[ig], nst
//
void PW_Basis::distribution_method1()
{

    // initial the variables needed by all proc.
    int tot_npw = 0;                     // total number of planewaves.
    this->nstot = 0;
    int st_start = 0;                    // index of the first stick on current proc.
    int *st_i = NULL;                    // x or x + nx (if x < 0) of stick.
    int *st_j = NULL;                    // y or y + ny (if y < 0) of stick.
    int *st_bottom = NULL;               // minimum z of stick.
    int *st_length = NULL;               // number of planewaves in stick.
    int *istot2ip = NULL;                   // ip of core that contains is^th stick, map is to ip.

    if (poolrank == 0)
    {
        // (1) Count the total number of planewaves (tot_npw) and sticks (this->nstot).                  
        
        // Actually we will scan [(2 * ibox[0] + 1) * (2 * ibox[1] + 1)] points on x-y plane,
        // but we define st_length2D with (nx * ny) points here, because we assume that the diameter
        // of the sphere is shorter than the sides of the cube.
        int *st_length2D = new int[this->nxy];                    // the number of planewaves that belong to the stick located on (x, y).
        int *st_bottom2D = new int[this->nxy];                    // the z-coordinate of the bottom of stick on (x, y).
        ModuleBase::GlobalFunc::ZEROS(st_length2D, this->nxy);
        ModuleBase::GlobalFunc::ZEROS(st_bottom2D, this->nxy);

        this->count_pw_st(tot_npw, st_length2D, st_bottom2D);
        // for test --------------------------------------------------
        std::cout << "The first step done\n";
        std::cout << "tot_npw   " << tot_npw << '\n';
        std::cout << "this->nstot   " << this->nstot << '\n';
        for (int iy = 0; iy < ny; ++iy)
        {
            for (int ix = 0; ix < nx; ++ix)
            {
                std::cout << st_length2D[iy * nx + ix] << std::setw(4);
            }
            std::cout << '\n';
        }
        // ------------------------------------------------------------

        // (2) Collect the x, y indexs, length, bottom of the sticks.
        st_i = new int[this->nstot];                           // x or x + nx (if x < 0) of stick.
        st_j = new int[this->nstot];                           // y or y + ny (if y < 0) of stick.
        st_bottom = new int[this->nstot];                      // minimum z of stick.
        st_length = new int[this->nstot];                      // number of planewaves in stick.
        
        this->collect_st(st_length2D, st_bottom2D, st_i, st_j, st_length, st_bottom);

        delete[] st_length2D;
        delete[] st_bottom2D;
        // for test ---------------------------------------------------
        std::cout << "st_i    ";
        for (int is = 0; is < this->nstot; ++is) cout << st_i[is] << setw(4) ;
        std::cout << "\nThe second step done\n";
        // ------------------------------------------------------------

#ifdef __MPI
        // Parallel line
        // (3) Distribute sticks to cores.
        int *npw_per = new int[this->poolnproc];  // number of planewaves on each core.
        int *nst_per = new int[this->poolnproc]; // number of sticks on each core.
        this->nstnz_per = new int[this->poolnproc]; // nz * nst(number of sticks) on each core.
        this->startnsz_per = new int[poolnproc];
        ModuleBase::GlobalFunc::ZEROS(npw_per, this->poolnproc);
        ModuleBase::GlobalFunc::ZEROS(nst_per, this->poolnproc);
        ModuleBase::GlobalFunc::ZEROS(this->nstnz_per, this->poolnproc);
        ModuleBase::GlobalFunc::ZEROS(startnsz_per, this->poolnproc);

        istot2ip = new int[this->nstot];                 // ip of core that contains is^th stick, map is to ip.
        for (int istot = 0; istot < this->nstot; ++istot)
        {
            istot2ip[istot] = -1;                        // meaning this stick has not been distributed.
        }
        this->ixy2ip = new int[this->nxy];              // ip of core which contains stick on (x, y).
        for (int ixy = 0; ixy < this->nxy; ++ixy)
        {
            this->ixy2ip[ixy] = -1;                 // meaning this stick has not been distributed or there is no stick on (x, y).
        }
        this->divide_sticks(tot_npw, st_i, st_j, st_length, npw_per, nst_per, istot2ip);
        // for test -----------------------------------------------------------------------------
        std::cout << "The 3-1 step done\n";
        // --------------------------------------------------------------------------------------

        this->get_istot2ixy(st_i, st_j, istot2ip);
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
        delete[] nst_per;
#else
        // Serial line
        this->npw = tot_npw;
        this->nst = this->nstot;

        this->nstnz_per = new int[1];
        this->nstnz_per[0] = this->nst * this->nz;
        this->startnsz_per = new int[1];
        this->startnsz_per[0] = 0;

        this->ixy2istot = new int[this->nxy];
        this->istot2ixy = new int[this->nstot];
        this->ixy2ip = new int[this->nxy];              // ip of core which contains stick on (x, y).
        for (int ixy = 0; ixy < this->nxy; ++ixy)
        {
            ixy2istot[ixy] = -1;
            ixy2ip[ixy] = -1;
        }
        for (int is = 0; is < this->nstot; ++is)
        {
            int ixy = st_i[is] + st_j[is] * nx;
            int ixy = st_i[is] + st_j[is] * nx;
            ixy2istot[ixy] = is;
            istot2ixy[is] = ixy;
            ixy2ip[ixy] = 0;
        }

        istot2ip = new int[this->nstot];
        ModuleBase::GlobalFunc::ZEROS(istot2ip, this->nstot);
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
    if (this->poolrank != 0)
    {
        st_i = new int[this->nstot];                           // x or x + nx (if x < 0) of stick.
        st_j = new int[this->nstot];                           // y or y + ny (if y < 0) of stick.
        st_bottom = new int[this->nstot];                      // minimum z of stick.
        st_length = new int[this->nstot];                      // number of planewaves in stick.
        istot2ip = new int[this->nstot];                 // ip of core that contains is^th stick, map is to ip.
        this->ixy2ip = new int[this->nxy];              // ip of core which contains stick on (x, y).
        this->ixy2istot = new int[this->nxy];
        this->istot2ixy = new int[this->nstot];
    }

    MPI_Bcast(st_i, this->nstot, MPI_INT, 0, POOL_WORLD);
    MPI_Bcast(st_j, this->nstot, MPI_INT, 0, POOL_WORLD);
    MPI_Bcast(st_length, this->nstot, MPI_INT, 0, POOL_WORLD);
    MPI_Bcast(st_bottom, this->nstot, MPI_INT, 0, POOL_WORLD);
    MPI_Bcast(istot2ip, this->nstot, MPI_INT, 0, POOL_WORLD);
    MPI_Bcast(ixy2ip, this->nxy, MPI_INT, 0, POOL_WORLD);
    MPI_Bcast(istot2ixy, this->nstot, MPI_INT, 0, POOL_WORLD);
    MPI_Bcast(ixy2istot, this->nxy, MPI_INT, 0, POOL_WORLD);
    
    std::cout << "Bcast done\n";
#endif
    this->nstnz = this->nst * this->nz;
    // for test ----------------------------------------------
    if (poolrank==0) std::cout << "The fifth step done\n";
    // -------------------------------------------------------

    // (5) Construct ig2isz, and is2ixy. 
    this->get_ig2isz_is2ixy(st_i, st_j, st_bottom, st_length, istot2ip);

    if (st_i != NULL) delete[] st_i;
    if (st_j != NULL) delete[] st_j;
    if (st_bottom != NULL) delete[] st_bottom;
    if (st_length != NULL) delete[] st_length;
    if (istot2ip != NULL) delete[] istot2ip;

    return;
}

//
// (1) We count the total number of planewaves (tot_npw) and sticks (this->nstot) here.
// Meanwhile, we record the number of planewaves on (x, y) in st_length2D, and store the smallest z-coordinate of each stick in st_bottom2D,
// so that we can scan a much smaller area in step(2).
// known: nx, ny, nz, ggecut, GGT
// output: tot_npw, this->nstot, st_length2D, st_bottom2D
//
void PW_Basis::count_pw_st(
        int &tot_npw,     // total number of planewaves.
        int* st_length2D, // the number of planewaves that belong to the stick located on (x, y).
        int* st_bottom2D  // the z-coordinate of the bottom of stick on (x, y).
)
{
    int ibox[3] = {0, 0, 0};                            // an auxiliary vector, determine the boundary of the scanning area.
    ibox[0] = int(this->nx / 2) + 1;                    // scan x from -ibox[0] to ibox[0].
    ibox[1] = int(this->ny / 2) + 1;                    // scan y from -ibox[1] to ibox[1].
    ibox[2] = int(this->nz / 2) + 1;                    // scan z from -ibox[2] to ibox[2].

    ModuleBase::Vector3<double> f;

    int ix_start = -ibox[0]; // determine the scaning area along x-direct, if gamma-only, only positive axis is used.
    int ix_end = ibox[0];
    if (this->gamma_only)
    {
        ix_start = 0;
        ix_end = this->nx;
    }
    for (int iy = -ibox[1]; iy <= ibox[1]; ++iy)
    {
        for (int ix = ix_start; ix <= ix_end; ++ix)
        {
            // we shift all sticks to the first quadrant in x-y plane here.
            // (ix, iy, iz) is the direct coordinates of planewaves.
            // x and y is the coordinates of shifted sticks in x-y plane.
            // for example, if nx = ny = 10, we will shift the stick on (-1, 2) to (9, 2),
            // so that its index in st_length and st_bottom is 9 + 10 * 2 = 29.
            int x = ix;
            int y = iy;
            if (x < 0) x += this->nx;
            if (y < 0) y += this->ny;
            int index = y * this->nx + x;

            int length = 0; // number of planewave in stick (x, y).
            for (int iz = -ibox[2]; iz <= ibox[2]; ++iz)
            {
                f.x = ix;
                f.y = iy;
                f.z = iz;
                double modulus = f * (this->GGT * f);
                if (modulus <= this->ggecut)
                {
                    if (length == 0) st_bottom2D[index] = iz; // length == 0 means this point is the bottom of stick (x, y).
                    ++tot_npw;
                    ++length;
                }
            }
            if (length > 0)
            {
                st_length2D[index] = length;
                ++this->nstot;
            }
        }
    }
    return;
}

//        
// (2) Collect the x, y indexs, length, z_start of the sticks.
// Firstly, we scan the area and construct temp_st_*.
// Then, as we will distribute the longest sticks preferentially in Step(3),
// we will sort temp_st_length from lagest to smallest, and reaarange st_* to the same order.
// known: tot_npw, this->nstot, st_length2D, st_bottom2D
// output: gg_global, gdirect_global, st_i, st_j, st_length, st_bottom
//
void PW_Basis::collect_st(
    int* st_length2D,                               // the number of planewaves that belong to the stick located on (x, y), stored in 2d x-y plane.
    int* st_bottom2D,                               // the z-coordinate of the bottom of stick on (x, y), stored in 2d x-y plane.
    int* st_i,                                      // x or x + nx (if x < 0) of stick.
    int* st_j,                                      // y or y + ny (if y < 0) of stick.
    int* st_length,                                 // number of planewaves in stick, stored in 1d array with this->nstot elements.
    int* st_bottom                                  // minimum z of stick, stored in 1d array with this->nstot elements.
)
{
    int *temp_st_i = new int[this->nstot];                      // x or x + nx (if x < 0) of stick.
    int *temp_st_j = new int[this->nstot];                      // y or y + ny (if y < 0) of stick.
    int *temp_st_bottom = new int[this->nstot];                 // minimum z of stick.
    double *temp_st_length = new double[this->nstot];           // length of sticks.
    ModuleBase::GlobalFunc::ZEROS(temp_st_length, this->nstot);

    int ibox[3] = {0, 0, 0};                            // an auxiliary vector, determine the boundary of the scanning area.
    ibox[0] = int(this->nx / 2) + 1;                    // scan x from -ibox[0] to ibox[0].
    ibox[1] = int(this->ny / 2) + 1;                    // scan y from -ibox[1] to ibox[1].
    ibox[2] = int(this->nz / 2) + 1;                    // scan z from -ibox[2] to ibox[2].

    ModuleBase::Vector3<double> f;
    int ig = 0; // index of planewave.
    int is = 0; // index of stick.

    int ix_start = -ibox[0]; // determine the scaning area along x-direct, if gamma-only, only positive axis is used.
    int ix_end = ibox[0];
    if (this->gamma_only)
    {
        ix_start = 0;
        ix_end = this->nx;
    }
    for (int iy = -ibox[1]; iy <= ibox[1]; ++iy)
    {
        for (int ix = ix_start; ix <= ix_end; ++ix)
        {
            // we have shifted all sticks to the first quadrant in x-y plane before.
            // (ix, iy, iz) is the direct coordinates of planewaves.
            // x and y is the coordinates of shifted sticks in x-y plane.
            // for example, if nx = ny = 10, we will shift the stick on (-1, 2) to (9, 2),
            // so that its index in st_length and st_bottom is 9 + 10 * 2 = 29.
            int x = ix;
            int y = iy;
            if (x < 0) x += this->nx;
            if (y < 0) y += this->ny;
            int index = y * nx + x;
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
                    temp_st_bottom[is] = st_bottom2D[index];
                    ++is;
                    std::cout << "is   " << is << '\n'; 
                }                  
            }
        }
    }
    assert(is == this->nstot);

    std::cout<<"collect sticks done\n";

    // As we will distribute the longest sticks preferentially in Step(3), we rearrange st_* in the order of length decreasing.

    int *st_sorted_index = new int[this->nstot]; // indexs in the order of length increasing.
    st_sorted_index[0] = 0;
    ModuleBase::heapsort(this->nstot, temp_st_length, st_sorted_index); // sort st_* in the order of length increasing.

    int index = 0;  // indexs in the order of length decreasing.
    for (int istot = 0; istot < this->nstot; ++istot)
    {
        index = (this->nstot - 1) - istot;
        st_length[index] = static_cast<int>(temp_st_length[istot]);
        st_i[index] = temp_st_i[st_sorted_index[istot]];
        st_j[index] = temp_st_j[st_sorted_index[istot]];
        st_bottom[index] = temp_st_bottom[st_sorted_index[istot]];
    }
    std::cout << "st_length    ";
    for (int is = 0; is < this->nstot; ++is) std::cout << st_length[is] << std::setw(4);
    std::cout << "\n";

    delete[] temp_st_i;
    delete[] temp_st_j;
    delete[] temp_st_bottom;
    delete[] st_sorted_index;
    return;
}

//
// (3-1) Distribute sticks to cores.
// We have rearranged sticks in the order of length decreasing, so that we will distribute the longest stick preferentially here.
// For each stick, we find the core that contains the least planewaves firstly, and distribute the stick to it,
// then update npw_per, this->nstnz_per, istot2ip this->ixy2ip and this->startnsz_per.
// known: tot_npw, this->nstot, st_i, st_j, st_length
// output: npw_per, nst_per, this->nstnz_per, istot2ip, this->ixy2ip, this->startnsz_per
//
void PW_Basis::divide_sticks(
    const int tot_npw,  // total number of planewaves.
    int* st_i,          // x or x + nx (if x < 0) of stick.
    int* st_j,          // y or y + ny (if y < 0) of stick.
    int* st_length,     // the stick on (x, y) consists of st_length[x*ny+y] planewaves.
    int* npw_per,       // number of planewaves on each core.
    int* nst_per,       // number of sticks on each core.
    int* istot2ip          // ip of core containing is^th stick, map is to ip.
)
{
    int ipmin = 0; // The ip of core containing least number of planewaves.
    for (int is = 0; is < this->nstot; ++is)
    {
        // find the ip of core containing the least planewaves.
        for (int ip = 0; ip < this->poolnproc; ++ip)
        {
            const int non_zero_grid = nst_per[ip] * this->nz;       // number of reciprocal planewaves on this core.
            const int npwmin = npw_per[ipmin];
            const int npw_ip = npw_per[ip];
            const int nstmin = nst_per[ipmin];
            const int nst_ip = nst_per[ip];

            if (npw_ip == 0)
            {
                if (non_zero_grid + nz < this->nrxx) // assert reciprocal planewaves is less than real space planewaves.
                {
                    ipmin = ip;
                    break;
                }
            }
            else if (npw_ip < npwmin)
            {
                if (non_zero_grid + nz < this->nrxx) // assert reciprocal planewaves is less than real space planewaves.
                {
                    ipmin = ip;
                }
            }
            else if (npw_ip == npwmin && nst_ip < nstmin)
            {
                if (non_zero_grid + nz < this->nrxx)
                {
                    ipmin = ip;
                }
            }
        }
        nst_per[ipmin]++;
        this->nstnz_per[ipmin] += this->nz;
        npw_per[ipmin] += st_length[is];
        istot2ip[is] = ipmin;
        this->ixy2ip[st_j[is] * this->nx + st_i[is]] = ipmin;
    }
    // for test --------------------------------------------------------------------------
    for(int i = 0; i < this->poolnproc; ++i) std::cout<<this->nstnz_per[i]<<std::setw(4);
    std::cout<<"\n";
    for(int i = 0; i < this->poolnproc; ++i) std::cout<<nst_per[i]<<std::setw(4);
    std::cout<<"\n";
    // -----------------------------------------------------------------------------------
    for (int ip = 1; ip < poolnproc; ++ip)
    {
        this->startnsz_per[ip] = this->startnsz_per[ip - 1] + this->nstnz_per[ip - 1];
    }
    return;
}

//
// (3-2) Rearrange sticks in the order of the ip of core increasing, (st_start + st_move) is the new index of sticks.
// Then get istot2ixy (istot2ixy[is]: ix + iy * nx of is^th stick among all sticks) on the first core
// and ixy2istot (ixy2istot[ix + iy * nx]: is of stick on (ix, iy) among all sticks).
// known: this->nstot, st_i, st_j, istot2ip, this->startnsz_per
// output: istot2ixy, ixy2istot
//
void PW_Basis::get_istot2ixy(
    int* st_i,          // x or x + nx (if x < 0) of stick.
    int* st_j,          // y or y + ny (if y < 0) of stick.
    int* istot2ip         // ip of core containing is^th stick, map is to ip.
)
{
    assert(this->poolrank == 0);
    this->ixy2istot = new int[this->nx * this->ny];
    this->istot2ixy = new int[this->nstot];
    int* st_move = new int[this->poolnproc]; // st_move[ip]: this is the st_move^th stick on ip^th core.
    for (int ixy = 0; ixy < this->nx * this->ny; ++ixy)
    {
        this->ixy2istot[ixy] = -1;
    }
    ModuleBase::GlobalFunc::ZEROS(this->istot2ixy, this->nstot);
    ModuleBase::GlobalFunc::ZEROS(st_move, this->poolnproc);


    for (int is = 0; is < this->nstot; ++is)
    {
        int ip = istot2ip[is];
        this->istot2ixy[this->startnsz_per[ip] / this->nz + st_move[ip]] = st_j[is] * nx + st_i[is];
        this->ixy2istot[st_j[is] * nx + st_i[is]] = this->startnsz_per[ip] / this->nz + st_move[ip];
        st_move[ip]++;
    }
    delete[] st_move;
    return;
}

//
// (5) Construct ig2isz, and is2ixy.
// is2ixy contains the x-coordinate and y-coordinate of sticks on current core.
// ig2isz contains the z-coordinate of planewaves on current core.
// We will scan all the sticks and find the planewaves on them, then store the information into ig2isz and is2ixy.
// known: this->nstot, st_i, st_j, st_bottom, st_length, istot2ip
// output: ig2isz, is2ixy
// 
void PW_Basis::get_ig2isz_is2ixy(
    int* st_i,          // x or x + nx (if x < 0) of stick.
    int* st_j,          // y or y + ny (if y < 0) of stick.
    int* st_bottom,     // minimum z of stick, stored in 1d array with this->nstot elements.
    int* st_length,     // the stick on (x, y) consists of st_length[x*ny+y] planewaves.
    int* istot2ip          // ip of core containing is^th stick, map is to ip.
)
{
    this->ig2isz = new int[this->npw]; // map ig to the z coordinate of this planewave.
    ModuleBase::GlobalFunc::ZEROS(this->ig2isz, this->npw);
    this->is2ixy = new int[this->nst]; // map is (index of sticks) to ixy (ix + iy * nx).
    for (int is = 0; is < this->nst; ++is) 
    {
        is2ixy[is] = -1;
    }

    ModuleBase::Vector3<double> f;
    int st_move = 0; // this is the st_move^th stick on current core.
    int pw_filled = 0; // how many current core's planewaves have been found.
    for (int is = 0; is < this->nstot; ++is)
    {
        if (istot2ip[is] == this->poolrank)
        {
            int zstart = st_bottom[is];
            for (int iz = zstart; iz < zstart + st_length[is]; ++iz)
            {
                int z = iz;
                if (z < 0) z += this->nz;
                this->ig2isz[pw_filled] = st_move * this->nz + z;
                pw_filled++;
            }
            this->is2ixy[st_move] = st_j[is] * this->nx + st_i[is];
            st_move++;
        }
        if (st_move == this->nst && pw_filled == this->npw) break;
    }
    return;
}


}
