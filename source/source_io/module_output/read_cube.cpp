#include "source_io/module_output/cube_io.h"
#include <cstdint>
#include <limits>
#include "source_base/parallel_grid.h"
#include "source_io/module_parameter/parameter.h"
#include <cstring>  // use std::memcpy

bool ModuleIO::read_vdata_palgrid(
    const Parallel_Grid& pgrid,
    const int my_rank,
    std::ostream& ofs_running,
    const std::string& fn,
    double* const data,
    const int natom)
{
    ModuleBase::TITLE("ModuleIO", "read_vdata_palgrid");

    // Only the root rank parses the file. On failure it must abort the whole
    // run instead of returning: the other ranks enter pgrid.bcast() below and
    // would block in MPI_Recv waiting for data that never comes.
    std::ifstream ifs(fn.c_str());
    if (my_rank == 0)
    {
        if (!ifs)
        {
            ofs_running << " !!! Couldn't find the file: " << fn << std::endl;
            ModuleBase::WARNING_QUIT("ModuleIO::read_vdata_palgrid",
                                     "couldn't find the cube file: " + fn);
        }
        ofs_running << " Find the file " << fn << " , try to read it." << std::endl;
    }

    // read the full grid data
    const int nx = pgrid.get_nx();
    const int ny = pgrid.get_ny();
    const int nz = pgrid.get_nz();
    const int& nxyz = nx * ny * nz;
    std::vector<double> data_xyz_full(nxyz, 0.0);
    if (my_rank == 0)
    {
        std::vector<std::string> comment;
        int natom = 0;
        std::vector<double> origin;
        std::vector<int> nvoxel;
        int nx_read = 0;
        int ny_read = 0;
        int nz_read = 0;
        std::vector<double> dx(3);
        std::vector<double> dy(3);
        std::vector<double> dz(3);
        std::vector<std::vector<double>> axis_vecs;
        std::vector<int> atom_type;
        std::vector<double> atom_charge;
        std::vector<std::vector<double>> atom_pos;
        std::vector<double> data_read;

        // validate the cube content before copying or interpolating the data
        if (!ModuleIO::read_cube(fn, comment, natom, origin, nx_read, ny_read, nz_read,
                                 dx, dy, dz, atom_type, atom_charge, atom_pos, data_read))
        {
            ofs_running << " !!! Failed to parse the cube file: " << fn << std::endl;
            ModuleBase::WARNING_QUIT("ModuleIO::read_vdata_palgrid",
                                     "failed to parse the cube file: " + fn);
        }

        // if mismatch, trilinear interpolate
        if (nx == nx_read && ny == ny_read && nz == nz_read)
        {
            std::memcpy(data_xyz_full.data(), data_read.data(), nxyz * sizeof(double));
        }
        else
        {
            trilinear_interpolate(data_read.data(), nx_read, ny_read, nz_read, nx, ny, nz, data_xyz_full.data());
        }
    }

    // distribute
#ifdef __MPI 
    pgrid.bcast(data_xyz_full.data(), data, my_rank, PARAM.inp.esolver_type == "sdft");
#else
    std::memcpy(data, data_xyz_full.data(), nxyz * sizeof(double));
#endif
    return true;
}

void ModuleIO::trilinear_interpolate(
    const double* const data_in,
    const int& nx_read,
    const int& ny_read,
    const int& nz_read,
    const int& nx,
    const int& ny,
    const int& nz,
    double* data_out)
{
    ModuleBase::TITLE("ModuleIO", "trilinear_interpolate");

    double** read_rho = new double*[nz_read];
    for (int iz = 0; iz < nz_read; iz++)
    {
        read_rho[iz] = new double[nx_read * ny_read];
    }
    for (int ix = 0; ix < nx_read; ix++)
    {
        for (int iy = 0; iy < ny_read; iy++)
        {
            for (int iz = 0; iz < nz_read; iz++)
            {
                read_rho[iz][ix * ny_read + iy] = data_in[(ix * ny_read + iy) * nz_read + iz];
            }
        }
    }

    for (int ix = 0; ix < nx; ix++)
    {
        double fracx = 0.5 * (static_cast<double>(nx_read) / nx * (1.0 + 2.0 * ix) - 1.0);
        fracx = std::fmod(fracx, nx_read);
        int lowx = static_cast<int>(fracx);
        double dx = fracx - lowx;
        int highx = (lowx == nx_read - 1) ? 0 : lowx + 1; // the point nz_read is the same as 0
        for (int iy = 0; iy < ny; iy++)
        {
            double fracy = 0.5 * (static_cast<double>(ny_read) / ny * (1.0 + 2.0 * iy) - 1.0);
            fracy = std::fmod(fracy, ny_read);
            int lowy = static_cast<int>(fracy);
            double dy = fracy - lowy;
            int highy = (lowy == ny_read - 1) ? 0 : lowy + 1;
            for (int iz = 0; iz < nz; iz++)
            {
                double fracz = 0.5 * (static_cast<double>(nz_read) / nz * (1.0 + 2.0 * iz) - 1.0);
                fracz = std::fmod(fracz, nz_read);
                int lowz = static_cast<int>(fracz);
                double dz = fracz - lowz;
                int highz = (lowz == nz_read - 1) ? 0 : lowz + 1;

                double result = read_rho[lowz][lowx * ny_read + lowy] * (1 - dx) * (1 - dy) * (1 - dz)
                                + read_rho[lowz][highx * ny_read + lowy] * dx * (1 - dy) * (1 - dz)
                                + read_rho[lowz][lowx * ny_read + highy] * (1 - dx) * dy * (1 - dz)
                                + read_rho[highz][lowx * ny_read + lowy] * (1 - dx) * (1 - dy) * dz
                                + read_rho[lowz][highx * ny_read + highy] * dx * dy * (1 - dz)
                                + read_rho[highz][highx * ny_read + lowy] * dx * (1 - dy) * dz
                                + read_rho[highz][lowx * ny_read + highy] * (1 - dx) * dy * dz
                                + read_rho[highz][highx * ny_read + highy] * dx * dy * dz;

                data_out[(ix * ny + iy) * nz + iz] = result;    // x > y > z order, consistent with the cube file
            }
        }
    }

    for (int iz = 0; iz < nz_read; iz++)
    {
        delete[] read_rho[iz];
    }
    delete[] read_rho;
}

bool ModuleIO::read_cube(const std::string& file,
    std::vector<std::string>& comment,
    int& natom,
    std::vector<double>& origin,
    int& nx,
    int& ny,
    int& nz,
    std::vector<double>& dx,
    std::vector<double>& dy,
    std::vector<double>& dz,
    std::vector<int>& atom_type,
    std::vector<double>& atom_charge,
    std::vector<std::vector<double>>& atom_pos,
    std::vector<double>& data)
{
    std::ifstream ifs(file);

    if (!ifs) 
    { 
	    return false; 
    }

    comment.resize(2);
    for (auto& c : comment) 
    { 
	    std::getline(ifs, c); 
    }

    ifs >> natom;
    if (ifs.fail() || natom < 0)
    {
        return false;
    }

    origin.resize(3);
    for (auto& cp : origin)
    {
        ifs >> cp;
    }
    if (ifs.fail())
    {
        return false;
    }

    dx.resize(3);
    dy.resize(3);
    dz.resize(3);
    ifs >> nx >> dx[0] >> dx[1] >> dx[2];
    ifs >> ny >> dy[0] >> dy[1] >> dy[2];
    ifs >> nz >> dz[0] >> dz[1] >> dz[2];
    if (ifs.fail() || nx <= 0 || ny <= 0 || nz <= 0)
    {
        return false;
    }

    atom_type.resize(natom);
    atom_charge.resize(natom);
    atom_pos.resize(natom, std::vector<double>(3));
    for (int i = 0; i < natom; ++i)
    {
        ifs >> atom_type[i] >> atom_charge[i] >> atom_pos[i][0] >> atom_pos[i][1] >> atom_pos[i][2];
    }
    if (ifs.fail())
    {
        return false;
    }

    // guard against int overflow before allocating the data buffer
    const std::int64_t nxyz_64 = static_cast<std::int64_t>(nx) * ny * nz;
    if (nxyz_64 > std::numeric_limits<int>::max())
    {
        return false;
    }
    const int nxyz = static_cast<int>(nxyz_64);
    data.resize(nxyz);
    for (int i = 0; i < nxyz; ++i)
    {
        ifs >> data[i];
    }
    if (ifs.fail())
    {
        return false;
    }

    ifs.close();
    return true;
}
