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
///                |  ---- ixy increasing ---> |   ---- ixy increasing --->   |...
/// index of sticks 0, 1, 2, ..., nst_per[0]-1, nst_per[0], ..., nst_per[1]-1, ...
///                |___________________________|______________________________|___
/// ip                           0                            1              ...
///                             npw    approximate equal to  npw   approximate equal to...
/// 
/// Known: G, GT, GGT, fftny, fftnx, nz, poolnproc, poolrank, ggecut
/// output: ig2isz[ig], istot2ixy[is], is2fftixy[is], fftixy2ip[ixy], startnsz_per[ip], nst_per[ip], nst
///
void PW_Basis::distribution_method1()
{
    // initial the variables needed by all process
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
        // but we define st_length2D with (fftny * fftnx) points here, because the diameter
        // of the sphere should be shorter than the sides of the cube.
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
        // (2) Collect the x, y indexs, and length of the sticks.
        int* st_i = new int[this->nstot];                           // x or x + fftnx (if x < 0) of stick.
        int* st_j = new int[this->nstot];                           // y or y + fftny (if y < 0) of stick.
        int* st_length = new int[this->nstot];                      // number of planewaves in stick.  
        this->collect_st(st_length2D, st_bottom2D, st_i, st_j, st_length);

        // (3) Distribute the sticks to cores.
        //get nst_per, npw_per, fftixy2ip, and startnsz_per
        this->startnsz_per = new int[this->poolnproc];
        this->divide_sticks_1(st_i, st_j, st_length);
        delete[] st_length;

        // (4) Get map from istot to (iy, ix)
        this->get_istot2ixy(st_i, st_j);
        delete[] st_i;
        delete[] st_j;
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
/// (2) Collect the x, y indexs, length of the sticks.
/// Firstly, we scan the area and construct temp_st_*.
/// Then, as we will distribute the longest sticks preferentially in Step(3),
/// we will sort temp_st_length from largest to smallest, and reaarange st_* to the same order.
/// known: tot_npw, this->nstot, st_length2D, st_bottom2D
/// output: st_i, st_j, st_length
///
void PW_Basis::collect_st(
    int* st_length2D,                               // the number of planewaves that belong to the stick located on (x, y), stored in 2d x-y plane.
    int* st_bottom2D,                               // the z-coordinate of the bottom of stick on (x, y), stored in 2d x-y plane.
    int* st_i,                                      // x or x + fftnx (if x < 0) of stick.
    int* st_j,                                      // y or y + fftny (if y < 0) of stick.
    int* st_length                                  // number of planewaves in stick, stored in 1d array with this->nstot elements.
)
{
    int *temp_st_i = new int[this->nstot];                      // x or x + fftnx (if x < 0) of stick.
    int *temp_st_j = new int[this->nstot];                      // y or y + fftny (if y < 0) of stick.
    double *temp_st_length = new double[this->nstot];           // length of sticks.
    ModuleBase::GlobalFunc::ZEROS(temp_st_length, this->nstot);

    ModuleBase::Vector3<double> f;
    int is = 0; // index of stick.

    int ix_end = int(this->nx / 2) + 1;
    int ix_start = -ix_end; 
    int iy_end = int(this->ny / 2) + 1;
    int iy_start = -iy_end; 
    if (this->full_pw)
    {
        ix_end = int(this->nx / 2);
        ix_start = ix_end - this->nx + 1; 
        iy_end = int(this->ny / 2);
        iy_start = iy_end - this->ny + 1; 
    }

    if (this->gamma_only)
    {
        if(this->xprime)
        {
            ix_start = 0;
            ix_end = this->fftnx - 1;
        }
        else
        {
            iy_start = 0;
            iy_end = this->fftny - 1;
        }
    }

    for (int ix = ix_start; ix <= ix_end; ++ix)
    {
        for (int iy = iy_start; iy <= iy_end; ++iy)
        {
            // we have shifted all sticks to the first quadrant in x-y plane before.
            // (ix, iy, iz) is the direct coordinates of planewaves.
            // x and y is the coordinates of shifted sticks in x-y plane.
            // for example, if fftnx = fftny = 10, we will shift the stick on (-1, 2) to (9, 2),
            // so that its index in st_length and st_bottom is 9 * 10 + 2 = 92.
            int x = ix;
            int y = iy;
            if (x < 0) x += nx;
            if (y < 0) y += ny;
            int index = x * this->fftny + y;
            if (st_length2D[index] > 0) // meaning there is a stick on (x, y) point.
            {
                bool find_stick = false;
                if (!this->full_pw)
                {
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
                }
                else
                {
                    find_stick = true;
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
    // std::cout<<"collect sticks done\n";

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

///
/// (3-1) Distribute sticks to cores according to the number of plane waves.
/// We have rearranged sticks in the order of length decreasing, so that we will distribute the longest stick preferentially here.
/// For each stick, we find the core that contains the least planewaves firstly, and distribute the stick to it,
/// then update npw_per, this->fftixy2ip, and this->startnsz_per.
/// known: tot_npw, this->nstot, st_i, st_j, st_length
/// output: npw_per, nst_per, this->fftixy2ip, this->startnsz_per
///
void PW_Basis::divide_sticks_1(
    int* st_i,          // x or x + fftnx (if x < 0) of stick.
    int* st_j,          // y or y + fftny (if y < 0) of stick.
    int* st_length     // the stick on (x, y) consists of st_length[x*fftny+y] planewaves.
)
{
    ModuleBase::GlobalFunc::ZEROS(this->nst_per, poolnproc);
    ModuleBase::GlobalFunc::ZEROS(this->npw_per, poolnproc);
    int ipmin = 0; // The ip of core containing least number of planewaves.
    for (int is = 0; is < this->nstot; ++is)
    {
        // find the ip of core containing the least planewaves.
        for (int ip = 0; ip < this->poolnproc; ++ip)
        {
            const int npwmin = this->npw_per[ipmin];
            const int npw_ip = this->npw_per[ip];
            const int nstmin = nst_per[ipmin];
            const int nst_ip = nst_per[ip];

            if (npw_ip == 0)
            {
                ipmin = ip;
                break;
            }
            else if (npw_ip < npwmin)
            {
                ipmin = ip;
            }
            else if (npw_ip == npwmin && nst_ip < nstmin)
            {
                ipmin = ip;
            }
        }
        this->nst_per[ipmin]++;
        this->npw_per[ipmin] += st_length[is];
        this->fftixy2ip[st_i[is] * this->fftny + st_j[is]] = ipmin;
    }

     this->startnsz_per[0] = 0;
    for (int ip = 1; ip < poolnproc; ++ip)
    {
        this->startnsz_per[ip] = this->startnsz_per[ip - 1] + this->nst_per[ip - 1] * this->nz;
    }
    return;
}

///
/// (3-2) Rearrange sticks in the order of the ip of core increasing, in each core, sticks are sorted in the order of ixy increasing.
/// (st_start + st_move) is the new index of sticks.
/// Then get istot2ixy (istot2ixy[is]: iy + ix * fftny of is^th stick among all sticks) on the first core
/// known: this->nstot, st_i, st_j, this->startnsz_per
/// output: istot2ixy
///
void PW_Basis::get_istot2ixy(
    int* st_i,          // x or x + fftnx (if x < 0) of stick.
    int* st_j          // y or y + fftny (if y < 0) of stick.
)
{
    assert(this->poolrank == 0);
    int* st_move = new int[this->poolnproc]; // st_move[ip]: this is the st_move^th stick on ip^th core.
    ModuleBase::GlobalFunc::ZEROS(this->istot2ixy, this->nstot);
    ModuleBase::GlobalFunc::ZEROS(st_move, this->poolnproc);

    for (int ixy = 0; ixy < this->fftnxy; ++ixy)
    {
        int ip = this->fftixy2ip[ixy];
        if (ip != -1)
        {        
            this->istot2ixy[this->startnsz_per[ip] / this->nz + st_move[ip]] = (ixy / fftny)*ny + ixy % fftny;
            st_move[ip]++;
        }
    }
    delete[] st_move;
    return;
}

}
