#include "pw_basis.h"
#include "../module_base/tool_quit.h"
#include "../module_base/global_function.h"
#include <iostream>

namespace ModulePW
{
//
//distribute plane waves to different cores
//Known: G, GT, GGT, nx, ny, nz, poolnproc, poolrank, ggecut
//output: ig2isz[ig], istot2ixy[is], ixy2istot[ixy], is2ixy[is], ixy2ip[ixy], startnsz_per[ip], nstnz_per[ip], gg[ig], gcar[ig], gdirect[ig], nst, nstot
//
void PW_Basis::distribute_g()
{
    if(this->distribution_type == 1)
    {
        this->distribution_method1();
    }
    else if(this->distribution_type == 2)
    {
        this->distribution_method2();
    }
    else
    {
        ModuleBase::WARNING_QUIT("divide", "No such division type.");
    }
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

    int iy_start = -ibox[1]; // determine the scaning area along x-direct, if gamma-only, only positive axis is used.
    int iy_end = ibox[1];
    if (this->gamma_only)
    {
        iy_start = 0;
        iy_end = this->ny - 1;
    }
    this->liy = this->riy = 0;
    for (int ix = -ibox[0]; ix <= ibox[0]; ++ix)
    {
        for (int iy = iy_start; iy <= iy_end; ++iy)
        {
            // we shift all sticks to the first quadrant in x-y plane here.
            // (ix, iy, iz) is the direct coordinates of planewaves.
            // x and y is the coordinates of shifted sticks in x-y plane.
            // for example, if ny = nx = 10, we will shift the stick on (-1, 2) to (9, 2),
            // so that its index in st_length and st_bottom is 9 * 10 + 2 = 92.
            int x = ix;
            int y = iy;
            if (x < 0) x += this->nx;
            if (y < 0) y += this->ny;
            int index = x * this->ny + y;

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
                    if(iy < this->riy) this->riy = iy;
                    if(iy > this->liy) this->liy = iy;
                }
            }
            if (length > 0)
            {
                st_length2D[index] = length;
                ++this->nstot;
            }
        }
    }
    riy += this->ny;
    std::cout<<"liy "<<liy<<" ; riy "<<riy<<std::endl;
    return;
}

//
// (5) Construct ig2isz, and is2ixy.
// is2ixy contains the x-coordinate and y-coordinate of sticks on current core.
// ig2isz contains the z-coordinate of planewaves on current core.
// We will scan all the sticks and find the planewaves on them, then store the information into ig2isz and is2ixy.
// known: this->nstot, st_bottom2D, st_length2D
// output: ig2isz, is2ixy
// 
void PW_Basis::get_ig2isz_is2ixy(
    int* st_bottom2D,     // minimum z of stick, stored in 1d array with this->nstot elements.
    int* st_length2D     // the stick on (x, y) consists of st_length[x*ny+y] planewaves.
)
{
    if (this->npw == 0)
    {
        this->ig2isz = new int[1]; // map ig to the z coordinate of this planewave.
        this->ig2isz[0] = 0;
        this->is2ixy = new int[1]; // map is (index of sticks) to ixy (ix + iy * nx).
        this->is2ixy[0] = -1;
        return;
    }

    this->ig2isz = new int[this->npw]; // map ig to the z coordinate of this planewave.
    ModuleBase::GlobalFunc::ZEROS(this->ig2isz, this->npw);
    this->is2ixy = new int[this->nst]; // map is (index of sticks) to ixy (ix + iy * nx).
    for (int is = 0; is < this->nst; ++is) 
    {
        this->is2ixy[is] = -1;
    }

    int st_move = 0; // this is the st_move^th stick on current core.
    int pw_filled = 0; // how many current core's planewaves have been found.
    for (int ixy = 0; ixy < this->nxy; ++ixy)
    {
        if (this->ixy2ip[ixy] == this->poolrank)
        {
            int zstart = st_bottom2D[ixy];
            for (int iz = zstart; iz < zstart + st_length2D[ixy]; ++iz)
            {
                int z = iz;
                if (z < 0) z += this->nz;
                this->ig2isz[pw_filled] = st_move * this->nz + z;
                pw_filled++;
                // assert(pw_filled <= this->npw);
            }
            this->is2ixy[st_move] = ixy;
            st_move++;
            // assert(st_move <= this->nst);
        }
        if (st_move == this->nst && pw_filled == this->npw) break;
    }
    return;
}
}
