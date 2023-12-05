#include "pw_basis.h"
#include "module_base/tool_quit.h"
#include "module_base/global_function.h"
#include "module_base/timer.h"

namespace ModulePW
{
///
/// distribute plane waves to different cores
/// Known: G, GT, GGT, fftnx, fftny, nz, poolnproc, poolrank, ggecut
/// output: ig2isz[ig], fftixyz2ig[isz], istot2ixy[is], is2fftixy[is], fftixy2ip[ixy], gg[ig], gcar[ig], gdirect[ig],
/// nst, nstot
///
void PW_Basis::distribute_g()
{
    ModuleBase::timer::tick(this->classname, "distributeg");
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
    ModuleBase::CHECK_WARNING_QUIT((this->npw == 0), "pw_distributeg.cpp", "Current core has no plane waves! Please reduce the cores.");
    ModuleBase::timer::tick(this->classname, "distributeg");
    return;
}

///
/// (1) We count the total number of planewaves (tot_npw) and sticks (this->nstot) here.
/// Meanwhile, we record the number of planewaves on (x, y) in st_length2D, and store the smallest z-coordinate of each stick in st_bottom2D,
/// so that we can scan a much smaller area in step(2).
/// known: fftnx, fftny, nz, ggecut, GGT
/// output: tot_npw, this->nstot, st_length2D, st_bottom2D, this->riy, this->liy
///
void PW_Basis::count_pw_st(
        int* st_length2D, // the number of planewaves that belong to the stick located on (x, y).
        int* st_bottom2D  // the z-coordinate of the bottom of stick on (x, y).
)
{
    ModuleBase::GlobalFunc::ZEROS(st_length2D, this->fftnxy);
    ModuleBase::GlobalFunc::ZEROS(st_bottom2D, this->fftnxy);
    ModuleBase::Vector3<double> f;

    // determine the scaning area along x-direct, if gamma-only && xprime, only positive axis is used.
    int ix_end = int(this->nx / 2) + 1;
    int ix_start = -ix_end; 
    // determine the scaning area along y-direct, if gamma-only && !xprime, only positive axis is used.
    int iy_end = int(this->ny / 2) + 1;
    int iy_start = -iy_end; 

    int iz_end = int(this->nz / 2) + 1;
    int iz_start = -iz_end;

    if (this->full_pw)
    {
        ix_end = int(this->nx / 2);
        ix_start = ix_end - this->nx + 1; 

        iy_end = int(this->ny / 2);
        iy_start = iy_end - this->ny + 1; 

        iz_end = int(this->nz / 2);
        iz_start = iz_end - this->nz + 1;
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

    this->liy = this->riy = 0;
    this->lix = this->rix = 0;
    this->npwtot = 0;
    this->nstot = 0;
    for (int ix = ix_start; ix <= ix_end; ++ix)
    {
        for (int iy = iy_start; iy <= iy_end; ++iy)
        {
            // we shift all sticks to the first quadrant in x-y plane here.
            // (ix, iy, iz) is the direct coordinates of planewaves.
            // x and y is the coordinates of shifted sticks in x-y plane.
            // for example, if fftny = fftnx = 10, we will shift the stick on (-1, 2) to (9, 2),
            // so that its index in st_length and st_bottom is 9 * 10 + 2 = 92.
            int x = ix;
            int y = iy;
            if (x < 0) x += this->nx;
            if (y < 0) y += this->ny;
            int index = x * this->fftny + y;

            int length = 0; // number of planewave on stick (x, y).
            for (int iz = iz_start; iz <= iz_end; ++iz)
            {
                f.x = ix;
                f.y = iy;
                f.z = iz;
                double modulus = f * (this->GGT * f);
                if (modulus <= this->ggecut || this->full_pw)
                {
                    if (length == 0) st_bottom2D[index] = iz; // length == 0 means this point is the bottom of stick (x, y).
                    ++this->npwtot;
                    ++length;
                    if(iy < this->riy) this->riy = iy;
                    if(iy > this->liy) this->liy = iy;
                    if(ix < this->rix) this->rix = ix;
                    if(ix > this->lix) this->lix = ix;
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
    rix += this->nx;
    return;
}

///
/// (5) Construct ig2isz, fftixyz2ig, and is2fftixy.
/// is2fftixy contains the x-coordinate and y-coordinate of sticks on current core.
/// ig2isz contains the z-coordinate of planewaves on current core.
/// fftixyz2ig maps ixyz to local index of planewaves.
/// scan all the sticks and find the planewaves on them, then store the information into ig2isz, fftixyz2ig and
/// is2fftixy. known: this->nstot, st_bottom2D, st_length2D output: ig2isz, fftixyz2ig, is2fftixy
///
void PW_Basis::get_ig2isz_is2fftixy(
    int* st_bottom2D,     // minimum z of stick, stored in 1d array with this->nstot elements.
    int* st_length2D     // the stick on (x, y) consists of st_length[x*fftny+y] planewaves.
)
{
    if (this->npw == 0)
    {
        delete[] this->ig2isz; this->ig2isz = nullptr; // map ig to the z coordinate of this planewave.
        delete[] this->fftixyz2ig;
        this->fftixyz2ig = nullptr;                          // map isz to ig.
        delete[] this->is2fftixy; this->is2fftixy = nullptr; // map is (index of sticks) to ixy (iy + ix * fftny).
#if defined(__CUDA) || defined(__ROCM)
        if (this->device == "gpu") {
            delmem_int_op()(gpu_ctx, this->d_is2fftixy);
            d_is2fftixy = nullptr;
        }
#endif
        return;
    }

    delete[] this->ig2isz; this->ig2isz = new int[this->npw]; // map ig to the z coordinate of this planewave.
    ModuleBase::GlobalFunc::ZEROS(this->ig2isz, this->npw);
    delete[] this->fftixyz2ig;
    this->fftixyz2ig = new int[this->nxyz]; // map isz to ig.
    for (int isz = 0; isz < this->nxyz; ++isz)
    {
        this->fftixyz2ig[isz] = -1;
    }
    delete[] this->is2fftixy; this->is2fftixy = new int[this->nst]; // map is (index of sticks) to ixy (iy + ix * fftny).
    for (int is = 0; is < this->nst; ++is) 
    {
        this->is2fftixy[is] = -1;
    }

    int st_move = 0; // this is the st_move^th stick on current core.
    int pw_filled = 0; // how many current core's planewaves have been found.
    for (int ixy = 0; ixy < this->fftnxy; ++ixy)
    {
        if (this->fftixy2ip[ixy] == this->poolrank)
        {
            int zstart = st_bottom2D[ixy];
            for (int iz = zstart; iz < zstart + st_length2D[ixy]; ++iz)
            {
                int z = iz;
                if (z < 0) z += this->nz;
                this->ig2isz[pw_filled] = st_move * this->nz + z;
                this->fftixyz2ig[ixy * this->nz + z] = pw_filled;
                pw_filled++;
            }
            this->is2fftixy[st_move] = ixy;
            st_move++;
            if(xprime && ixy/fftny == 0) ng_xeq0 = pw_filled;
        }
        if (st_move == this->nst && pw_filled == this->npw) break;
    }
#if defined(__CUDA) || defined(__ROCM)
    if (this->device == "gpu") {
        resmem_int_op()(gpu_ctx, d_is2fftixy, this->nst);
        syncmem_int_h2d_op()(gpu_ctx, cpu_ctx, this->d_is2fftixy, this->is2fftixy, this->nst);
    }
#endif
    return;
}
}
