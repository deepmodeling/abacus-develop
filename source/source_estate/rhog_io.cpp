#include "source_base/module_out/binstream.h"
#include "source_base/vector3.h"
#include "source_base/module_parallel/para_mpi_func.h"
#include "rhog_io.h"
#include <algorithm>
#include <numeric>
#include <unistd.h>

namespace
{
inline void warn(std::ostream* os,
                 const Parallel::ParaWorld& pw_world,
                 const std::string& file,
                 const std::string& desc)
{
    if (pw_world.rank() == 0 && os != nullptr)
    {
        *os << " " << file << "  warning : " << desc << std::endl;
    }
}
} // namespace

bool elecstate::read_rhog(const std::string& filename,
                         const ModulePW::PW_Basis* pw_rhod,
                         const int nspin,
                         std::complex<double>** rhog,
                         const Parallel::ParaWorld& pw_world,
                         std::ostream* os_warning)
{
    if (pw_rhod == nullptr)
    {
        warn(os_warning, pw_world, "elecstate::read_rhog", "pw_rhod is null");
        return false;
    }
    if (rhog == nullptr)
    {
        warn(os_warning, pw_world, "elecstate::read_rhog", "rhog is null");
        return false;
    }
    if (nspin != 1 && nspin != 2 && nspin != 4)
    {
        warn(os_warning, pw_world, "elecstate::read_rhog", "nspin must be 1, 2, or 4");
        return false;
    }
    if (pw_rhod->nx <= 0 || pw_rhod->ny <= 0 || pw_rhod->nz <= 0)
    {
        warn(os_warning, pw_world, "elecstate::read_rhog", "PW_Basis grid dimensions must be positive");
        return false;
    }

    const int nx = pw_rhod->nx;
    const int ny = pw_rhod->ny;
    const int nz = pw_rhod->nz;

    Binstream ifs;
    bool error = false;
    int gamma_only_in = 0;
    int npwtot_in = 0;
    int nspin_in = 0;
    int size = 0;
    double b1[3], b2[3], b3[3];

    if (pw_world.rank() == 0)
    {
        ifs.open(filename, "r");
        if (!ifs)
        {
            error = true;
        }
    }

    Parallel::bcast_bool(error, pw_world);

    if (error)
    {
        warn(os_warning, pw_world, "elecstate::read_rhog", "Can't open file " + filename);
        return false;
    }

    if (pw_world.rank() == 0)
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
            warn(os_warning, pw_world, "elecstate::read_rhog", "some planewaves in file are not used");
        }
        else if (npwtot_in < pw_rhod->npwtot)
        {
            warn(os_warning, pw_world, "elecstate::read_rhog", "some planewaves in file are missing");
        }
        if (nspin_in < nspin)
        {
            warn(os_warning, pw_world, "elecstate::read_rhog", "some spin channels in file are missing");
        }
    }

    Parallel::bcast_bool(error, pw_world);

    if (error)
    {
        warn(os_warning, pw_world, "elecstate::read_rhog", "gamma_only read from file is inconsistent with INPUT");
        return false;
    }

    Parallel::bcast_int(gamma_only_in, pw_world);
    Parallel::bcast_int(npwtot_in, pw_world);
    Parallel::bcast_int(nspin_in, pw_world);
    Parallel::bcast_double(b1, 3, pw_world);
    Parallel::bcast_double(b2, 3, pw_world);
    Parallel::bcast_double(b3, 3, pw_world);

    std::vector<int> miller(npwtot_in * 3);
    // once use ModuleBase::Vector3, it is highly bug-prone to assume the memory layout of the class.
    // The x, y and z of Vector3 will not always to be contiguous.
    // Instead, a relatively safe choice is to use std::vector, the memory layout is assumed
    // to be npwtot_in rows and 3 columns.
    if (pw_world.rank() == 0)
    {
        ifs >> size;
        for (int i = 0; i < npwtot_in; ++i) // loop over rows...
        {
            ifs >> miller[i*3] >> miller[i*3+1] >> miller[i*3+2];
        }
        ifs >> size;
    }
    Parallel::bcast_int(miller.data(), miller.size(), pw_world);

    // set to zero
    for (int is = 0; is < nspin; ++is)
    {
        std::fill(rhog[is], rhog[is] + pw_rhod->npw, std::complex<double>(0.0, 0.0));
    }
    // maps ixyz tp ig
    std::vector<int> fftixyz2ig(pw_rhod->nxyz, -1); // map isz to ig.
    for (int ig = 0; ig < pw_rhod->npw; ++ig)
    {
        int isz = pw_rhod->ig2isz[ig];
        int iz = isz % nz;
        int is = isz / nz;
        int ixy = pw_rhod->is2fftixy[is];
        int ixyz = iz + nz * ixy;
        fftixyz2ig[ixyz] = ig;
    }
    std::vector<std::complex<double>> rhog_in(npwtot_in);
    for (int is = 0; is < nspin_in; ++is)
    {
        if (pw_world.rank() == 0)
        {
            ifs >> size;
            for (int i = 0; i < npwtot_in; ++i)
            {
                ifs >> rhog_in[i];
            }
            ifs >> size;
        }
        Parallel::bcast_complex(rhog_in.data(), rhog_in.size(), pw_world);

        for (int i = 0; i < npwtot_in; ++i)
        {
            int ix = miller[i * 3];
            int iy = miller[i * 3 + 1];
            int iz = miller[i * 3 + 2];

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
            if (pw_world.rank() == pw_rhod->fftixy2ip[fftixy])
            {
                int fftixyz = iz + nz * fftixy;
                int ig = fftixyz2ig[fftixyz];
                rhog[is][ig] = rhog_in[i];
            }
        }

        if (nspin_in == 2 && nspin == 4 && is == 1)
        {
            for (int ig = 0; ig < pw_rhod->npw; ++ig)
            {
                rhog[3][ig] = rhog[1][ig];
            }
            std::fill(rhog[1], rhog[1] + pw_rhod->npw, std::complex<double>(0.0, 0.0));
            std::fill(rhog[2], rhog[2] + pw_rhod->npw, std::complex<double>(0.0, 0.0));
        }
    }

    if (pw_world.rank() == 0)
    {
        ifs.close();
    }
    return true;
}

bool elecstate::write_rhog(const std::string& fchg,
                          const bool gamma_only,
                          const ModulePW::PW_Basis* pw_rho,
                          const int nspin,
                          const ModuleBase::Matrix3& GT,
                          std::complex<double>** rhog,
                          const Parallel::ParaWorld& pw_world,
                          std::ostream* os_warning)
{
    if (pw_rho == nullptr)
    {
        warn(os_warning, pw_world, "elecstate::write_rhog", "pw_rho is null");
        return false;
    }
    if (rhog == nullptr)
    {
        warn(os_warning, pw_world, "elecstate::write_rhog", "rhog is null");
        return false;
    }
    if (nspin != 1 && nspin != 2 && nspin != 4)
    {
        warn(os_warning, pw_world, "elecstate::write_rhog", "nspin must be 1, 2, or 4");
        return false;
    }

    // only rank 0 in the domain writes the header; all ranks cooperate
    // on sequential writes synchronized by barriers.
    const int irank = pw_world.rank();
    const int nrank = pw_world.size();

    // write the header (by rank 0): gamma_only, ngm_g, nspin
    int size = 3;
    int ngm_g = pw_rho->npwtot;
    int gam = gamma_only;
    int nsp = nspin;

    std::ofstream ofs;
    Parallel::barrier(pw_world);

    if (irank == 0)
    {
        ofs.open(fchg, std::ios::binary);
        if (!ofs)
        {
            warn(os_warning, pw_world, "elecstate::write_rhog", "File I/O failure: cannot open file " + fchg);
            return false;
        }
        ofs.write(reinterpret_cast<char*>(&size), sizeof(size));
        ofs.write(reinterpret_cast<char*>(&gam), sizeof(gam));
        ofs.write(reinterpret_cast<char*>(&ngm_g), sizeof(ngm_g));
        ofs.write(reinterpret_cast<char*>(&nsp), sizeof(nsp));
        ofs.write(reinterpret_cast<char*>(&size), sizeof(size));
        // write the lattice vectors
        std::vector<double> b = {GT.e11, GT.e12, GT.e13, GT.e21, GT.e22, GT.e23, GT.e31, GT.e32, GT.e33};
        size = 9;
        ofs.write(reinterpret_cast<char*>(&size), sizeof(size));
        for (int i = 0; i < 9; ++i)
        {
            ofs.write(reinterpret_cast<char*>(&b[i]), sizeof(b[i]));
        }
        ofs.write(reinterpret_cast<char*>(&size), sizeof(size));
        ofs.close();
    }
    Parallel::barrier(pw_world);
    Parallel::barrier(pw_world);

    // write the G-vectors in Miller indices
    size = 3 * ngm_g;
    if (irank == 0)
    {
        ofs.open(fchg, std::ios::binary | std::ios::app);
        ofs.write(reinterpret_cast<char*>(&size), sizeof(size));
        ofs.close();
    }
    Parallel::barrier(pw_world);

    for (int i = 0; i < nrank; ++i)
    {
        if (i == irank)
        {
            ofs.open(fchg, std::ios::binary | std::ios::app);
            for (int ig = 0; ig < pw_rho->npw; ++ig)
            {
                const ModuleBase::Vector3<double> g = pw_rho->gdirect[ig];
                std::vector<int> miller = {int(g.x), int(g.y), int(g.z)};
                ofs.write(reinterpret_cast<char*>(&miller[0]), sizeof(miller[0]));
                ofs.write(reinterpret_cast<char*>(&miller[1]), sizeof(miller[1]));
                ofs.write(reinterpret_cast<char*>(&miller[2]), sizeof(miller[2]));
            }
            ofs.close();
        }
        Parallel::barrier(pw_world);
    }

    if (irank == 0)
    {
        ofs.open(fchg, std::ios::binary | std::ios::app);
        ofs.write(reinterpret_cast<char*>(&size), sizeof(size));
        ofs.close();
    }
    Parallel::barrier(pw_world);

    // write the rho(G) values
    std::complex<double> sum_check;
    size = ngm_g;
    for (int ispin = 0; ispin < nspin; ++ispin)
    {
        if (irank == 0)
        {
            ofs.open(fchg, std::ios::binary | std::ios::app);
            ofs.write(reinterpret_cast<char*>(&size), sizeof(size));
            ofs.close();
        }
        Parallel::barrier(pw_world);

        for (int i = 0; i < nrank; ++i)
        {
            if (i == irank)
            {
                ofs.open(fchg, std::ios::binary | std::ios::app);
                sum_check = 0.0;
                for (int ig = 0; ig < pw_rho->npw; ++ig)
                {
                    sum_check += rhog[ispin][ig];
                    ofs.write(reinterpret_cast<char*>(&rhog[ispin][ig]), sizeof(rhog[ispin][ig]));
                }
                ofs.close();
            }
            Parallel::barrier(pw_world);
        }

        if (irank == 0)
        {
            ofs.open(fchg, std::ios::binary | std::ios::app);
            ofs.write(reinterpret_cast<char*>(&size), sizeof(size));
            ofs.close();
        }
        Parallel::barrier(pw_world);
    }
    return true;
}

// self-consistency test with the following python code
// import numpy as np

// with open("rhog_read.txt") as f:
//     read = f.readlines()

// with open("rhog_write.txt") as f:
//     write = f.readlines()

// # convert c++ stype complex number (a,b) to python complex
// def to_complex(s):
//     a, b = s.replace("(", "").replace(")", "").split(",")
//     return complex(float(a), float(b))

// read = [[to_complex(rhog) for rhog in spin.strip().split()] for spin in read]
// write = [[to_complex(rhog) for rhog in spin.strip().split()] for spin in write]

// diff = np.array(read) - np.array(write)
// print(np.max(np.abs(diff)))
// test system: integrated test 118_PW_CHG_BINARY
// yielding error 5.290000000000175e-11