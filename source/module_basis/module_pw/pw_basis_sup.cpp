#include "pw_basis_sup.h"

#include "module_base/timer.h"

namespace ModulePW
{

PW_Basis_Sup::~PW_Basis_Sup()
{
    delete[] igs2igd;
    delete[] igd2igs;
}

///
/// distribute plane wave basis and real-space grids to different processors
/// set up maps for fft and create arrays for MPI_Alltoall
/// set up ffts
///
void PW_Basis_Sup::setuptransform(const int* fftixy2ip_s, // fftixy2ip of smooth grids
                                  const int& nx_s,        // nx of smooth grids
                                  const int& ny_s         // ny of smooth grids
)
{
    ModuleBase::timer::tick(this->classname, "setuptransform");
    this->distribute_r();
    this->distribute_g(fftixy2ip_s, nx_s, ny_s);
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
/// output: ig2isz[ig], istot2ixy[is], is2fftixy[is], fftixy2ip[ixy], gg[ig], gcar[ig], gdirect[ig], nst, nstot
///
void PW_Basis_Sup::distribute_g(const int* fftixy2ip_s, // fftixy2ip of smooth grids
                                const int& nx_s,        // nx of smooth grids
                                const int& ny_s         // ny of smooth grids
)
{
    ModuleBase::timer::tick(this->classname, "distributeg");
    this->distribution_method3(fftixy2ip_s, nx_s, ny_s);
    ModuleBase::CHECK_WARNING_QUIT((this->npw == 0),
                                   "pw_distributeg.cpp",
                                   "Current core has no plane waves! Please reduce the cores.");
    ModuleBase::timer::tick(this->classname, "distributeg");
    return;
}

void PW_Basis_Sup::link_igs_igd(const ModuleBase::Vector3<double>* gcar_s, // G vectors  of smooth grids
                                const int& npw_s                           // npw of smooth grids
)
{
    delete[] igs2igd;
    delete[] igd2igs;
    igs2igd = new int[npw_s];
    igd2igs = new int[this->npw];
    for (int i = 0; i < npw_s; ++i)
    {
        igs2igd[i] = -1;
    }
    for (int i = 0; i < this->npw; ++i)
    {
        igd2igs[i] = -1;
    }
    for (int igs = 0; igs < npw_s; igs++)
    {
        for (int igd = 0; igd < this->npw; igd++)
        {
            if (gcar_s[igs] == this->gcar[igd])
            {
                igs2igd[igs] = igd;
                igd2igs[igd] = igs;
            }
        }
    }
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
/// Known: fftixy2ip[ixy], nx, ny of smooth grids
/// Known: G, GT, GGT, fftny, fftnx, nz, poolnproc, poolrank, ggecut
/// output: ig2isz[ig], istot2ixy[is], is2fftixy[is], fftixy2ip[ixy], startnsz_per[ip], nst_per[ip], nst
///
void PW_Basis_Sup::distribution_method3(const int* fftixy2ip_s, const int& nx_s, const int& ny_s)
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
        this->divide_sticks_3(st_length2D, st_i, st_j, st_length, fftixy2ip_s, nx_s, ny_s);
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

    // (5) Construct ig2isz and is2fftixy.
    this->get_ig2isz_is2fftixy(st_bottom2D, st_length2D);

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

} // namespace ModulePW