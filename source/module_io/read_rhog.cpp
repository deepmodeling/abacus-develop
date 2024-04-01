#include "binstream.h"
#include "module_base/global_function.h"
#include "module_base/global_variable.h"
#include "module_base/parallel_global.h"
#include "module_base/timer.h"
#include "module_base/vector3.h"
#include "rhog_io.h"

bool ModuleIO::read_rhog(const std::string& filename, const ModulePW::PW_Basis* pw_rhod, std::complex<double>** rhog)
{
    ModuleBase::TITLE("ModuleIO", "read_rhog");
    ModuleBase::timer::tick("ModuleIO", "read_rhog");

    const int nx = pw_rhod->nx;
    const int ny = pw_rhod->ny;
    const int nz = pw_rhod->nz;

    Binstream ifs;
    bool error = false;
    int gamma_only_in, npwtot_in, nspin_in, size;

    double b1[3], b2[3], b3[3];

    if (GlobalV::RANK_IN_POOL == 0)
    {
        ifs.open(filename, "r");
        if (!ifs)
        {
            error = true;
        }
    }

#ifdef __MPI
    MPI_Bcast(&error, 1, MPI_C_BOOL, 0, POOL_WORLD);
#endif

    if (error)
    {
        ModuleBase::WARNING("ModuleIO::read_rhog", "Can't open file " + filename);
        return false;
    }

    if (GlobalV::RANK_IN_POOL == 0)
    {
        ifs >> size >> gamma_only_in >> npwtot_in >> nspin_in >> size;
        ifs >> size >> b1[0] >> b1[1] >> b1[2] >> b2[0] >> b2[1] >> b2[2] >> b3[0] >> b3[1] >> b3[2] >> size;

        if (gamma_only_in != pw_rhod->gamma_only)
        {
            // there is a treatment that can transform between gamma_only and non-gamma_only
            // however, it is not implemented here
            error = true;
            ifs.close();
        }
        if (npwtot_in > pw_rhod->npwtot)
        {
            ModuleBase::WARNING("ModuleIO::read_rhog", "some planewaves in file are not used");
        }
        else if (npwtot_in < pw_rhod->npwtot)
        {
            ModuleBase::WARNING("ModuleIO::read_rhog", "some planewaves in file are missing");
        }
        if (nspin_in < GlobalV::NSPIN)
        {
            ModuleBase::WARNING("ModuleIO::read_rhog", "some spin channels in file are missing");
        }
    }

#ifdef __MPI
    MPI_Bcast(&error, 1, MPI_C_BOOL, 0, POOL_WORLD);
#endif

    if (error)
    {
        ModuleBase::WARNING("ModuleIO::read_rhog", "gamma_only read from file is inconsistent with INPUT");
        return false;
    }

#ifdef __MPI
    MPI_Bcast(&gamma_only_in, 1, MPI_INT, 0, POOL_WORLD);
    MPI_Bcast(&npwtot_in, 1, MPI_INT, 0, POOL_WORLD);
    MPI_Bcast(&nspin_in, 1, MPI_INT, 0, POOL_WORLD);
    MPI_Bcast(b1, 3, MPI_DOUBLE, 0, POOL_WORLD);
    MPI_Bcast(b2, 3, MPI_DOUBLE, 0, POOL_WORLD);
    MPI_Bcast(b3, 3, MPI_DOUBLE, 0, POOL_WORLD);
#endif

    ModuleBase::Vector3<int>* miller = new ModuleBase::Vector3<int>[npwtot_in];
    if (GlobalV::RANK_IN_POOL == 0)
    {
        ifs >> size;
        for (int i = 0; i < npwtot_in; ++i)
        {
            ifs >> miller[i].x >> miller[i].y >> miller[i].z;
        }
        ifs >> size;
    }
#ifdef __MPI
    MPI_Bcast(miller, 3 * npwtot_in, MPI_INT, 0, POOL_WORLD);
#endif

    // set to zero
    for (int is = 0; is < GlobalV::NSPIN; ++is)
    {
        ModuleBase::GlobalFunc::ZEROS(rhog[is], pw_rhod->npw);
    }

    // maps ixyz tp ig
    int* fftixyz2ig = new int[pw_rhod->nxyz]; // map isz to ig.
    for (int i = 0; i < pw_rhod->nxyz; ++i)
    {
        fftixyz2ig[i] = -1;
    }
    for (int ig = 0; ig < pw_rhod->npw; ++ig)
    {
        int isz = pw_rhod->ig2isz[ig];
        int iz = isz % nz;
        int is = isz / nz;
        int ixy = pw_rhod->is2fftixy[is];
        int ixyz = iz + nz * ixy;
        fftixyz2ig[ixyz] = ig;
    }

    std::complex<double>* rhog_in = new std::complex<double>[npwtot_in];
    for (int is = 0; is < nspin_in; ++is)
    {
        if (GlobalV::RANK_IN_POOL == 0)
        {
            ifs >> size;
            for (int i = 0; i < npwtot_in; ++i)
            {
                ifs >> rhog_in[i];
            }
            ifs >> size;
        }
#ifdef __MPI
        MPI_Bcast(rhog_in, npwtot_in, MPI_DOUBLE_COMPLEX, 0, POOL_WORLD);
#endif

        for (int i = 0; i < npwtot_in; ++i)
        {
            int ix = miller[i].x;
            int iy = miller[i].y;
            int iz = miller[i].z;

            if (ix <= -int((nx + 1) / 2) || ix >= int(nx / 2) + 1 || iy <= -int((ny + 1) / 2) || iy >= int(ny / 2) + 1
                || iz <= -int((nz + 1) / 2) || iz >= int(nz / 2) + 1)
            {
                // these planewaves are not used
                continue;
            }

            if (ix < 0)
                ix += nx;
            if (iy < 0)
                iy += ny;
            if (iz < 0)
                iz += nz;
            int fftixy = iy + pw_rhod->fftny * ix;
            if (GlobalV::RANK_IN_POOL == pw_rhod->fftixy2ip[fftixy])
            {
                int fftixyz = iz + nz * fftixy;
                int ig = fftixyz2ig[fftixyz];
                rhog[is][ig] = rhog_in[i];
            }
        }

        if (nspin_in == 2 && GlobalV::NSPIN == 4 && is == 1)
        {
            for (int ig = 0; ig < pw_rhod->npw; ++ig)
            {
                rhog[3][ig] = rhog[1][ig];
            }
            ModuleBase::GlobalFunc::ZEROS(rhog[1], pw_rhod->npw);
            ModuleBase::GlobalFunc::ZEROS(rhog[2], pw_rhod->npw);
        }
    }

    delete[] fftixyz2ig;
    delete[] miller;
    delete[] rhog_in;

    if (GlobalV::RANK_IN_POOL == 0)
    {
        ifs.close();
    }

    ModuleBase::timer::tick("ModuleIO", "read_rhog");
    return true;
}
