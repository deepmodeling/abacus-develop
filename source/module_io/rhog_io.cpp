#include "binstream.h"
#include "module_base/global_function.h"
#include "module_base/global_variable.h"
#include "module_base/parallel_global.h"
#include "module_base/timer.h"
#include "module_base/vector3.h"
#include "rhog_io.h"
#include <numeric>
#include <unistd.h>

bool ModuleIO::read_rhog(const std::string& filename, const ModulePW::PW_Basis* pw_rhod, std::complex<double>** rhog)
{
    ModuleBase::TITLE("ModuleIO", "read_rhog");
    ModuleBase::timer::tick("ModuleIO", "read_rhog");

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
    std::vector<int> miller(npwtot_in * 3); 
    // once use ModuleBase::Vector3, it is highly bug-prone to assume the memory layout of the class. 
    // The x, y and z of Vector3 will not always to be contiguous.
    // Instead, a relatively safe choice is to use std::vector, the memory layout is assumed
    // to be npwtot_in rows and 3 columns.
    if (GlobalV::RANK_IN_POOL == 0)
    {
        ifs >> size;
        for (int i = 0; i < npwtot_in; ++i) // loop over rows...
        {
            ifs >> miller[i*3] >> miller[i*3+1] >> miller[i*3+2];
        }
        ifs >> size;
    }
#ifdef __MPI
    MPI_Bcast(miller.data(), miller.size(), MPI_INT, 0, POOL_WORLD);
#endif
    // set to zero
    for (int is = 0; is < GlobalV::NSPIN; ++is)
    {
        ModuleBase::GlobalFunc::ZEROS(rhog[is], pw_rhod->npw);
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
        MPI_Bcast(rhog_in.data(), rhog_in.size(), MPI_DOUBLE_COMPLEX, 0, POOL_WORLD);
#endif

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

    if (GlobalV::RANK_IN_POOL == 0)
    {
        ifs.close();
    }
    // for debug, write the rhog to a file (not binary)
    // if (GlobalV::RANK_IN_POOL == 0)
    // {
    //     std::ofstream ofs("rhog_read.txt");
    //     for (int i = 0; i < nspin_in; ++i)
    //     {
    //         for (int ig = 0; ig < pw_rhod->npw; ++ig)
    //         {
    //             ofs << rhog[i][ig] << " ";
    //         }
    //         ofs << std::endl;
    //     }
    //     ofs.close();
    // }
    ModuleBase::timer::tick("ModuleIO", "read_rhog");
    return true;
}

bool ModuleIO::write_rhog(const std::string& fchg,
                          const bool gamma_only, // from INPUT
                          const ModulePW::PW_Basis* pw_rho, // pw_rho in runtime
                          const int nspin, // GlobalV
                          const ModuleBase::Matrix3& GT, // from UnitCell, useful for calculating the miller
                          std::complex<double>** rhog,
                          const int ipool,
                          const int irank,
                          const int nrank)
{
    ModuleBase::TITLE("ModuleIO", "write_rhog");
    ModuleBase::timer::tick("ModuleIO", "write_rhog");
    if (ipool != 0) { return true; }
    // only one pool writes the rhog, because rhog in all pools are identical.

    // for large-scale data, it is not wise to collect all distributed components to the
    // master process and then write the data to the file. Instead, we can write the data
    // processer by processer.

    // Quantum ESPRESSO requires the G-vectors collected should be in the order like as if
    // there is only 1 process, this order is recorded in

    // fftixy2ip will be useful for the order of the G-vectors
    // each time we iterate on ig, then find the rho_g over all the processes
    // for ig in npwtot, then find the local index of ig on processor, ig -> fftixy2ip -> igl


    // write the header (by rank 0): gamma_only, ngm_g, nspin
    int size = 3;
    // because "reinterpret_cast" cannot drop the "const", so use intermediate variable
    int ngm_g = pw_rho->npwtot;
    int gam = gamma_only; // IMPLICIT DATA TYPE CONVERSION!
    int nsp = nspin;

    std::ofstream ofs;
#ifdef __MPI
    MPI_Barrier(POOL_WORLD); 
    // this is still a global variable... should be moved into param
    // list as `const MPI_Comm& comm`
    if (irank == 0)
    {
        // printf(" CHGDEN >>> Writing header by rank %d...\n", irank);
#endif
    ofs.open(fchg, std::ios::binary); // open the file by all processors
    if (!ofs)
    {
        ModuleBase::WARNING_QUIT("ModuleIO::write_rhog", "File I/O failure: cannot open file " + fchg);
        return false;
    }
    ofs.write(reinterpret_cast<char*>(&size), sizeof(size));
    ofs.write(reinterpret_cast<char*>(&gam), sizeof(gam));
    ofs.write(reinterpret_cast<char*>(&ngm_g), sizeof(ngm_g));
    ofs.write(reinterpret_cast<char*>(&nsp), sizeof(nsp));
    ofs.write(reinterpret_cast<char*>(&size), sizeof(size));
    // write the lattice vectors, GT is the reciprocal lattice vectors, need 2pi?
    std::vector<double> b = {GT.e11, GT.e12, GT.e13, GT.e21, GT.e22, GT.e23, GT.e31, GT.e32, GT.e33};
    size = 9;
    ofs.write(reinterpret_cast<char*>(&size), sizeof(size));
    for (int i = 0; i < 9; ++i)
    {
        ofs.write(reinterpret_cast<char*>(&b[i]), sizeof(b[i]));
    }
    ofs.write(reinterpret_cast<char*>(&size), sizeof(size));
    ofs.close();
#ifdef __MPI
    // printf(" CHGDEN >>> Complete header writting by rank %d\n", irank);
    }
    MPI_Barrier(POOL_WORLD); // wait for rank 0 to finish writing the header
    // printf(" CHGDEN >>> rank %d ready for continue writing...\n", irank);
    MPI_Barrier(POOL_WORLD);
#endif

    // write the G-vectors in Miller indices, the Miller indices can be calculated by
    // the dot product of the G-vectors and the reciprocal lattice vectors
    // parallelization needed considered here. Because the sequence of the G-vectors
    // is not important, we can write the G-vectors processer by processer
    size = 3*ngm_g;
#ifdef __MPI
    if(irank == 0)
    {
        // printf(" CHGDEN >>> Writing header of Miller indices by rank %d...\n", irank);
#endif
        ofs.open(fchg, std::ios::binary | std::ios::app); // open the file by rank 0
        ofs.write(reinterpret_cast<char*>(&size), sizeof(size));
        ofs.close();
#ifdef __MPI
        // printf(" CHGDEN >>> Complete header of Miller indices writting by rank %d\n", irank);
    }
    MPI_Barrier(POOL_WORLD); // wait for rank 0 to finish writing the header of miller indices
#endif
#ifdef __MPI
    for(int i = 0; i < nrank; ++i) // write the miller indices processer by processer
    {
        if(i == irank)
        {
            // printf(" CHGDEN >>> Writing Miller indices by rank %d...\n", irank);
#endif
            ofs.open(fchg, std::ios::binary | std::ios::app); // open the file by processer i
            for(int ig = 0; ig < pw_rho->npw; ++ig)
            {
                const ModuleBase::Vector3<double> g = pw_rho->gdirect[ig]; // g direct is (ix, iy, iz), miller index (integer), centered at (0, 0, 0)
                std::vector<int> miller = {int(g.x), int(g.y), int(g.z)};
                ofs.write(reinterpret_cast<char*>(&miller[0]), sizeof(miller[0]));
                ofs.write(reinterpret_cast<char*>(&miller[1]), sizeof(miller[1]));
                ofs.write(reinterpret_cast<char*>(&miller[2]), sizeof(miller[2]));
            }
            ofs.close();
#ifdef __MPI
            // printf(" CHGDEN >>> Complete Miller indices writting by rank %d\n", irank);
        }
        MPI_Barrier(POOL_WORLD); // wait for the current rank to finish writing the miller indices
    }
#endif
#ifdef __MPI
    if(irank == 0)
    {
#endif
        ofs.open(fchg, std::ios::binary | std::ios::app); // open the file by rank 0
        ofs.write(reinterpret_cast<char*>(&size), sizeof(size));
        ofs.close();
#ifdef __MPI
    }
    MPI_Barrier(POOL_WORLD); // wait for rank 0 to finish writing the miller indices
#endif

    // write the rho(G) values
    std::complex<double> sum_check;
    size = ngm_g;
    for(int ispin = 0; ispin < nspin; ++ispin)
    {
#ifdef __MPI
        if(irank == 0)
        {
            // printf(" CHGDEN >>> Writing header of rho(G) values by rank %d...\n", irank);
#endif
            ofs.open(fchg, std::ios::binary | std::ios::app); // open the file by rank 0
            ofs.write(reinterpret_cast<char*>(&size), sizeof(size));
            ofs.close();
#ifdef __MPI
            // printf(" CHGDEN >>> Complete header of rho(G) values writting by rank %d\n", irank);
        }
        MPI_Barrier(POOL_WORLD); // wait for rank 0 to finish writing the header of rho(G)
#endif
#ifdef __MPI
        for(int i = 0; i < nrank; ++i) // write the rho(G) values processer by processer
        {
            if(i == irank)
            {
                // printf(" CHGDEN >>> Writing rho(G) values by rank %d...\n", irank);
#endif
                ofs.open(fchg, std::ios::binary | std::ios::app); // open the file by processer i
                sum_check = 0.0;
                for(int ig = 0; ig < pw_rho->npw; ++ig)
                {
                    sum_check += rhog[ispin][ig];
                    ofs.write(reinterpret_cast<char*>(&rhog[ispin][ig]), sizeof(rhog[ispin][ig]));
                }
                // assert(std::abs(sum_check) > 1.0e-10); // check if the sum of rho(G) is valid
                ofs.close();
#ifdef __MPI
                // printf(" CHGDEN >>> Complete rho(G) values writting by rank %d\n", irank);
            }
            MPI_Barrier(POOL_WORLD); // wait for the current rank to finish writing the rho(G) values
        }
#endif

#ifdef __MPI
        if(irank == 0)
        {
#endif
            ofs.open(fchg, std::ios::binary | std::ios::app); // open the file by rank 0
            ofs.write(reinterpret_cast<char*>(&size), sizeof(size));
            ofs.close();
#ifdef __MPI
        }
        MPI_Barrier(POOL_WORLD); // wait for rank 0 to finish writing the rho(G) values
#endif
    }
    // for debug, write the rhog to a file (not binary)
    // if (irank == 0)
    // {
    //     std::ofstream ofs("rhog_write.txt");
    //     for (int i = 0; i < nspin; ++i)
    //     {
    //         for (int ig = 0; ig < pw_rho->npw; ++ig)
    //         {
    //             ofs << rhog[i][ig] << " ";
    //         }
    //         ofs << std::endl;
    //     }
    //     ofs.close();
    // }
    ModuleBase::timer::tick("ModuleIO", "write_rhog");
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