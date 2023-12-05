#include "module_base/timer.h"
#include "pw_basis.h"

namespace ModulePW
{

PW_Basis_Sup::~PW_Basis_Sup()
{
}

///
/// distribute plane wave basis and real-space grids to different processors
/// set up maps for fft and create arrays for MPI_Alltoall
/// set up ffts
///
void PW_Basis_Sup::setuptransform(const ModulePW::PW_Basis* pw_rho)
{
    ModuleBase::timer::tick(this->classname, "setuptransform");
    this->distribute_r();
    this->distribute_g(pw_rho);
    this->getstartgr();
    this->ft.clear();
    if (this->xprime)
        this->ft.initfft(this->nx,
                         this->ny,
                         this->nz,
                         this->lix,
                         this->rix,
                         this->nst,
                         this->nplane,
                         this->poolnproc,
                         this->gamma_only,
                         this->xprime);
    else
        this->ft.initfft(this->nx,
                         this->ny,
                         this->nz,
                         this->liy,
                         this->riy,
                         this->nst,
                         this->nplane,
                         this->poolnproc,
                         this->gamma_only,
                         this->xprime);
    this->ft.setupFFT();
    ModuleBase::timer::tick(this->classname, "setuptransform");
}

///
/// distribute plane waves to different cores
/// Known: G, GT, GGT, fftnx, fftny, nz, poolnproc, poolrank, ggecut
/// output: ig2isz[ig], fftixyz2ig[isz], istot2ixy[is], is2fftixy[is], fftixy2ip[ixy], gg[ig], gcar[ig], gdirect[ig],
/// nst, nstot
///
void PW_Basis_Sup::distribute_g(const ModulePW::PW_Basis* pw_rho)
{
    ModuleBase::timer::tick(this->classname, "distributeg");
    this->distribution_method3(pw_rho);
    ModuleBase::CHECK_WARNING_QUIT((this->npw == 0),
                                   "pw_distributeg.cpp",
                                   "Current core has no plane waves! Please reduce the cores.");
    ModuleBase::timer::tick(this->classname, "distributeg");
    return;
}

///
/// Distribute planewaves in reciprocal space to cores.
/// Firstly, divide the sphere in reciprocal space into sticks, which are vertical to x-y plane.
/// Secondly, distribute these sticks to cores.
///
/// Example
///                |  ---- ixy increasing ---> |   ---- ixy increasing --->   |...
/// index of sticks 0, 1, 2, ..., nst_per[0]-1, nst_per[0], ..., nst_per[1]-1, ...
///                |___________________________|______________________________|___
/// ip                           0                            1              ...
///                             npw    approximate equal to  npw   approximate equal to...
///
/// Note: This method are ONLY used for dense grids in uspp, and it is not suitable for other cases.
///       The smooth grids is constructed by distribution_method1().
///       Then, in order to conserve the consistence of planewaves between dense and smooth grids,
///       we divide sticks corresponding to smooth girds first, and then the left ones are divided
///       to ensure the approximate equality of planewaves on each core.
///
/// Known: smooth grids
/// Known: G, GT, GGT, fftny, fftnx, nz, poolnproc, poolrank, ggecut
/// output: ig2isz[ig], fftixyz2ig[isz], istot2ixy[is], is2fftixy[is], fftixy2ip[ixy], startnsz_per[ip], nst_per[ip],
/// nst
///
void PW_Basis_Sup::distribution_method3(const ModulePW::PW_Basis* pw_rho)
{
    // initial the variables needed by all process
    int* st_bottom2D = new int[fftnxy]; // st_bottom2D[ixy], minimum z of stick on (x, y).
    int* st_length2D = new int[fftnxy]; // st_length2D[ixy], number of planewaves in stick on (x, y).
    delete[] this->nst_per;
    this->nst_per = new int[this->poolnproc]; // number of sticks on each core.
    delete[] this->npw_per;
    this->npw_per = new int[this->poolnproc]; // number of planewaves on each core.
    delete[] this->fftixy2ip;
    this->fftixy2ip = new int[this->fftnxy]; // ip of core which contains the stick on (x, y).
    for (int ixy = 0; ixy < this->fftnxy; ++ixy)
        this->fftixy2ip[ixy] = -1; // meaning this stick has not been distributed or there is no stick on (x, y).
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
    delete[] this->istot2ixy;
    this->istot2ixy = new int[this->nstot];

    if (poolrank == 0)
    {
#ifdef __MPI
        // Parallel line
        // (2) Collect the x, y indexs, and length of the sticks.
        int* st_i = new int[this->nstot];      // x or x + fftnx (if x < 0) of stick.
        int* st_j = new int[this->nstot];      // y or y + fftny (if y < 0) of stick.
        int* st_length = new int[this->nstot]; // number of planewaves in stick.
        this->collect_st(st_length2D, st_bottom2D, st_i, st_j, st_length);

        // (3) Distribute the sticks to cores.
        // get nst_per, npw_per, fftixy2ip, and startnsz_per
        this->startnsz_per = new int[this->poolnproc];
        this->divide_sticks_3(st_length2D, st_i, st_j, st_length, pw_rho->fftixy2ip, pw_rho->nx, pw_rho->ny);
        delete[] st_length;

        // (4) Get map from istot to (iy, ix)
        this->get_istot2ixy(st_i, st_j);
        delete[] st_i;
        delete[] st_j;
        // We do not need startnsz_per after it.
        delete[] this->startnsz_per;
        this->startnsz_per = nullptr;
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
    MPI_Bcast(this->nst_per, this->poolnproc, MPI_INT, 0, this->pool_world);
    MPI_Bcast(this->npw_per, this->poolnproc, MPI_INT, 0, this->pool_world);
#endif
    this->npw = this->npw_per[this->poolrank];
    this->nst = this->nst_per[this->poolrank];
    this->nstnz = this->nst * this->nz;

    // (5) Construct ig2isz, fftixyz2ig and is2fftixy.
    this->get_ig2isz_is2fftixy(st_bottom2D, st_length2D, pw_rho);

    delete[] st_bottom2D;
    delete[] st_length2D;
    return;
}

///
/// (3-1) Distribute sticks to cores.
/// The smooth grids is constructed by distribution_method1().
/// Then, in order to conserve the consistence of planewaves between dense and smooth grids,
/// we divide sticks corresponding to smooth girds first.
/// We have rearranged sticks in the order of length decreasing, so that we will distribute the longest in the lefted
/// stick preferentially here. For each stick, we find the core that contains the least planewaves firstly, and
/// distribute the stick to it, then update npw_per, this->fftixy2ip, and this->startnsz_per.
/// known: fftixy2ip[ixy], fftnxy, fftny, nx, ny of smooth grids
/// known: tot_npw, this->nstot, st_i, st_j, st_length
/// output: npw_per, nst_per, this->fftixy2ip, this->startnsz_per
///
void PW_Basis_Sup::divide_sticks_3(
    const int* st_length2D, // st_length2D[ixy], number of planewaves in stick on (x, y).
    const int* st_i,        // x or x + fftnx (if x < 0) of stick.
    const int* st_j,        // y or y + fftny (if y < 0) of stick.
    const int* st_length,   // the stick on (x, y) consists of st_length[x*fftny+y] planewaves.
    const int* fftixy2ip_s, // fftixy2ip of smooth grids
    const int& nx_s,        // nx of smooth grids
    const int& ny_s)        // ny of smooth grids
{
    ModuleBase::GlobalFunc::ZEROS(this->nst_per, poolnproc);
    ModuleBase::GlobalFunc::ZEROS(this->npw_per, poolnproc);

    int fftny_s = ny_s;
    int fftnx_s = nx_s;
    if (this->gamma_only)
    {
        if (this->xprime)
            fftnx_s = int(nx_s / 2) + 1;
        else
            fftny_s = int(ny_s / 2) + 1;
    }

    int fftnxy_s = fftnx_s * fftny_s;

    // (1) Distribute sticks corresponding to smooth grids first.
    for (int ixy = 0; ixy < fftnxy_s; ++ixy)
    {
        int ix = ixy / fftny_s;
        int iy = ixy % fftny_s;
        if (ix >= int(nx_s / 2) + 1)
            ix -= nx_s;
        if (iy >= int(ny_s / 2) + 1)
            iy -= ny_s;

        if (ix < 0)
            ix += nx;
        if (iy < 0)
            iy += ny;
        int index = ix * this->fftny + iy;
        int ip = fftixy2ip_s[ixy];
        if (ip >= 0)
        {
            this->fftixy2ip[index] = ip;
            this->nst_per[ip]++;
            this->npw_per[ip] += st_length2D[index];
        }
    }

    // distribute the longest in the lefted stick preferentially.
    int ipmin = 0; // The ip of core containing least number of planewaves.
    for (int is = 0; is < this->nstot; ++is)
    {
        // skip sticks corresponding to smooth grids.
        if (this->fftixy2ip[st_i[is] * this->fftny + st_j[is]] >= 0)
        {
            continue;
        }

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
/// (5) Construct ig2isz, fftixyz2ig and is2fftixy.
/// is2fftixy contains the x-coordinate and y-coordinate of sticks on current core.
/// ig2isz contains the z-coordinate of planewaves on current core.
/// fftixyz2ig maps ixyz to local index of planewaves.
/// We will scan all the sticks and find the planewaves on them, then store the information into ig2isz and is2fftixy.
/// known: smooth grids
/// known: this->nstot, st_bottom2D, st_length2D
/// output: ig2isz, fftixyz2ig, is2fftixy
///
void PW_Basis_Sup::get_ig2isz_is2fftixy(
    int* st_bottom2D, // minimum z of stick, stored in 1d array with this->nstot elements.
    int* st_length2D, // the stick on (x, y) consists of st_length[x*fftny+y] planewaves.
    const ModulePW::PW_Basis* pw_rho)
{
    if (this->npw == 0)
    {
        delete[] this->ig2isz;
        this->ig2isz = nullptr; // map ig to the z coordinate of this planewave.
        delete[] this->fftixyz2ig;
        this->fftixyz2ig = nullptr; // map isz to ig.
        delete[] this->is2fftixy;
        this->is2fftixy = nullptr; // map is (index of sticks) to ixy (iy + ix * fftny).
#if defined(__CUDA) || defined(__ROCM)
        if (this->device == "gpu")
        {
            delmem_int_op()(gpu_ctx, this->d_is2fftixy);
            d_is2fftixy = nullptr;
        }
#endif
        return;
    }

    delete[] this->ig2isz;
    this->ig2isz = new int[this->npw]; // map ig to the z coordinate of this planewave.
    ModuleBase::GlobalFunc::ZEROS(this->ig2isz, this->npw);
    delete[] this->fftixyz2ig;
    this->fftixyz2ig = new int[this->nxyz]; // map isz to ig.
    for (int isz = 0; isz < this->nxyz; ++isz)
    {
        this->fftixyz2ig[isz] = -1;
    }
    delete[] this->is2fftixy;
    this->is2fftixy = new int[this->nst]; // map is (index of sticks) to ixy (iy + ix * fftny).
    for (int is = 0; is < this->nst; ++is)
    {
        this->is2fftixy[is] = -1;
    }
    int* fftixy2is = new int[this->fftnxy]; // map ixy to is.
    for (int ixy = 0; ixy < this->fftnxy; ++ixy)
    {
        fftixy2is[ixy] = -1;
    }
    bool* found = new bool[this->fftnxyz]; // whether the planewave on (x, y, z) has been found on the smooth grid.
    for (int i = 0; i < this->fftnxyz; ++i)
    {
        found[i] = false;
    }

    // get is2fftixy
    int st_move = 0; // this is the st_move^th stick on current core.
    for (int ixy = 0; ixy < this->fftnxy; ++ixy)
    {
        if (this->fftixy2ip[ixy] == this->poolrank)
        {
            this->is2fftixy[st_move] = ixy;
            fftixy2is[ixy] = st_move;
            st_move++;
        }
        if (st_move == this->nst)
            break;
    }

    // distribute planewaves in the same order as smooth grids first.
    int pw_filled = 0; // how many current core's planewaves have been found.
    for (int ig = 0; ig < pw_rho->npw; ig++)
    {
        int isz = pw_rho->ig2isz[ig];
        int iz = isz % pw_rho->nz;
        int is = isz / pw_rho->nz;
        int ixy = pw_rho->is2fftixy[is];
        int ix = ixy / pw_rho->fftny;
        int iy = ixy % pw_rho->fftny;
        if (ix >= int(pw_rho->nx / 2) + 1)
            ix -= pw_rho->nx;
        if (iy >= int(pw_rho->ny / 2) + 1)
            iy -= pw_rho->ny;
        if (iz >= int(pw_rho->nz / 2) + 1)
            iz -= pw_rho->nz;

        if (ix < 0)
            ix += this->nx;
        if (iy < 0)
            iy += this->ny;
        if (iz < 0)
            iz += this->nz;
        int ixy_now = ix * this->fftny + iy;
        int index = ixy_now * this->nz + iz;
        int is_now = fftixy2is[ixy_now];
        int isz_now = is_now * this->nz + iz;
        this->ig2isz[ig] = isz_now;
        this->fftixyz2ig[ixy_now * this->nz + iz] = ig;
        pw_filled++;
        found[index] = true;
        if (xprime && ix == 0)
            ng_xeq0++;
    }
    assert(pw_filled == pw_rho->npw);

    // distribute the lefted planewaves.
    for (int ixy = 0; ixy < this->fftnxy; ++ixy)
    {
        if (this->fftixy2ip[ixy] == this->poolrank)
        {
            int zstart = st_bottom2D[ixy];
            for (int iz = zstart; iz < zstart + st_length2D[ixy]; ++iz)
            {
                int z = iz;
                if (z < 0)
                    z += this->nz;
                if (!found[ixy * this->nz + z])
                {
                    found[ixy * this->nz + z] = true;
                    int is = fftixy2is[ixy];
                    this->ig2isz[pw_filled] = is * this->nz + z;
                    this->fftixyz2ig[ixy * this->nz + z] = pw_filled;
                    pw_filled++;
                    if (xprime && ixy / fftny == 0)
                        ng_xeq0++;
                }
            }
        }
        if (pw_filled == this->npw)
            break;
    }

    delete[] fftixy2is;
    delete[] found;

#if defined(__CUDA) || defined(__ROCM)
    if (this->device == "gpu")
    {
        resmem_int_op()(gpu_ctx, d_is2fftixy, this->nst);
        syncmem_int_h2d_op()(gpu_ctx, cpu_ctx, this->d_is2fftixy, this->is2fftixy, this->nst);
    }
#endif
    return;
}

} // namespace ModulePW