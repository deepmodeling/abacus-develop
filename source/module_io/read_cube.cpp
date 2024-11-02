#include "module_io/cube_io.h"
#include <limits>
#include "module_hamilt_pw/hamilt_pwdft/parallel_grid.h"
// #include "module_base/global_variable.h" // GlobalV reference removed

bool ModuleIO::read_vdata_palgrid(
    const Parallel_Grid& pgrid,
    const int my_rank,
    std::ofstream& ofs_running,
    const std::string& fn,
    double* const data,
    const int natom)
{
    ModuleBase::TITLE("ModuleIO", "read_vdata_palgrid");

    // check if the file exists
    std::ifstream ifs(fn.c_str());
    if (!ifs)
    {
        std::string tmp_warning_info = "!!! Couldn't find the file: " + fn;
        ofs_running << tmp_warning_info << std::endl;
        return false;
    }
    else
    {
        ofs_running << " Find the file " << fn << " , try to read it." << std::endl;
    }

    // read the full grid data
    const int& nx = pgrid.nx;
    const int& ny = pgrid.ny;
    const int& nz = pgrid.nz;
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

        // we've already checked the file existence, so we don't need the returned value here
        ModuleIO::read_cube(fn, comment, natom, origin, nx_read, ny_read, nz_read, dx, dy, dz, atom_type, atom_charge, atom_pos, data_read);

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
    pgrid.bcast(data_xyz_full.data(), data, my_rank);
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

    if (!ifs) { return false; }

    comment.resize(2);
    for (auto& c : comment) { std::getline(ifs, c); }

    ifs >> natom;
    origin.resize(3);
    for (auto& cp : origin) { ifs >> cp; }

    dx.resize(3);
    dy.resize(3);
    dz.resize(3);
    ifs >> nx >> dx[0] >> dx[1] >> dx[2];
    ifs >> ny >> dy[0] >> dy[1] >> dy[2];
    ifs >> nz >> dz[0] >> dz[1] >> dz[2];

    atom_type.resize(natom);
    atom_charge.resize(natom);
    atom_pos.resize(natom, std::vector<double>(3));
    for (int i = 0;i < natom;++i)
    {
        ifs >> atom_type[i] >> atom_charge[i] >> atom_pos[i][0] >> atom_pos[i][1] >> atom_pos[i][2];
    }

    const int nxyz = nx * ny * nz;
    data.resize(nxyz);
    for (int i = 0;i < nxyz;++i) { ifs >> data[i]; }

    ifs.close();
    return true;
}