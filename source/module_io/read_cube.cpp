#include "module_io/cube_io.h"
// #include "module_base/global_variable.h" // GlobalV reference removed

bool ModuleIO::read_cube(
#ifdef __MPI
    const Parallel_Grid*const Pgrid,
#endif
    const int my_rank,
    const std::string esolver_type,
    const int rank_in_stogroup,
    const int is,
    std::ofstream& ofs_running,
    const int nspin,
    const std::string& fn,
    double*const data,
    const int nx,
    const int ny,
    const int nz,
    double& ef,
    const UnitCell*const ucell,
    int& prenspin,
    const bool warning_flag)
{
    ModuleBase::TITLE("ModuleIO","read_cube");
    std::ifstream ifs(fn.c_str());
    if (!ifs)
    {
        std::string tmp_warning_info = "!!! Couldn't find the charge file of ";
        tmp_warning_info += fn;
        ofs_running << tmp_warning_info << std::endl;
        return false;
    }
    else
    {
        ofs_running << " Find the file, try to read charge from file." << std::endl;
    }

    bool quit=false;

    ifs.ignore(300, '\n'); // skip the header

    if(nspin != 4)
    {
        int v_in;
        ifs >> v_in;
        if (v_in != nspin)
        {
            std::cout << " WARNING: nspin mismatch:  " << nspin << " in INPUT parameters but " << v_in << " in " << fn << std::endl;
            return false;
        }
    }
    else
    {
        ifs >> prenspin;
    }
    ifs.ignore(150, ')');

    ifs >> ef;
    ofs_running << " read in fermi energy = " << ef << std::endl;

    ifs.ignore(150, '\n');

    ModuleBase::CHECK_INT(ifs, ucell->nat);
    ifs.ignore(150, '\n');

    int nx_read = 0;
    int ny_read = 0;
    int nz_read = 0;
    double fac = ucell->lat0;
    std::string temp;
    if (warning_flag)
    {
        ifs >> nx_read;
        ModuleBase::CHECK_DOUBLE(ifs, fac * ucell->latvec.e11 / double(nx), quit);
        ModuleBase::CHECK_DOUBLE(ifs, fac * ucell->latvec.e12 / double(nx), quit);
        ModuleBase::CHECK_DOUBLE(ifs, fac * ucell->latvec.e13 / double(nx), quit);
        ifs >> ny_read;
        ModuleBase::CHECK_DOUBLE(ifs, fac * ucell->latvec.e21 / double(ny), quit);
        ModuleBase::CHECK_DOUBLE(ifs, fac * ucell->latvec.e22 / double(ny), quit);
        ModuleBase::CHECK_DOUBLE(ifs, fac * ucell->latvec.e23 / double(ny), quit);
        ifs >> nz_read;
        ModuleBase::CHECK_DOUBLE(ifs, fac * ucell->latvec.e31 / double(nz), quit);
        ModuleBase::CHECK_DOUBLE(ifs, fac * ucell->latvec.e32 / double(nz), quit);
        ModuleBase::CHECK_DOUBLE(ifs, fac * ucell->latvec.e33 / double(nz), quit);
    }
    else
    {
        ifs >> nx_read;
        ifs >> temp >> temp >> temp;
        ifs >> ny_read;
        ifs >> temp >> temp >> temp;
        ifs >> nz_read;
        ifs >> temp >> temp >> temp;
    }

    for (int it = 0; it < ucell->ntype; it++)
    {
        for (int ia = 0; ia < ucell->atoms[it].na; ia++)
        {
            ifs >> temp; // skip atomic number
            ifs >> temp; // skip Z valance
            if (warning_flag)
            {
                ModuleBase::CHECK_DOUBLE(ifs, fac * ucell->atoms[it].tau[ia].x, quit);
                ModuleBase::CHECK_DOUBLE(ifs, fac * ucell->atoms[it].tau[ia].y, quit);
                ModuleBase::CHECK_DOUBLE(ifs, fac * ucell->atoms[it].tau[ia].z, quit);
            }
            else
            {
                ifs >> temp >> temp >> temp;
            }
        }
    }

    const bool flag_read_rank = (my_rank == 0 || (esolver_type == "sdft" && rank_in_stogroup == 0));
#ifdef __MPI
    if(nx == nx_read && ny == ny_read && nz == nz_read)
        ModuleIO::read_cube_core_match(ifs, Pgrid, flag_read_rank, data, nx*ny, nz);
    else
        ModuleIO::read_cube_core_mismatch(ifs, Pgrid, flag_read_rank, data, nx, ny, nz, nx_read, ny_read, nz_read);
#else
    ofs_running << " Read SPIN = " << is + 1 << " charge now." << std::endl;
    if(nx == nx_read && ny == ny_read && nz == nz_read)
        ModuleIO::read_cube_core_match(ifs, data, nx*ny, nz);
    else
        ModuleIO::read_cube_core_mismatch(ifs, data, nx, ny, nz, nx_read, ny_read, nz_read);
#endif

    return true;
}

void ModuleIO::read_cube_core_match(
    std::ifstream &ifs,
#ifdef __MPI
    const Parallel_Grid*const Pgrid,
    const bool flag_read_rank,
#endif
    double*const data,
    const int nxy,
    const int nz)
{
#ifdef __MPI
    if (flag_read_rank)
    {
        std::vector<std::vector<double>> read_rho(nz, std::vector<double>(nxy));
        for (int ixy = 0; ixy < nxy; ixy++)
            for (int iz = 0; iz < nz; iz++)
                ifs >> read_rho[iz][ixy];
        for (int iz = 0; iz < nz; iz++)
            Pgrid->zpiece_to_all(read_rho[iz].data(), iz, data);
    }
    else
    {
        std::vector<double> zpiece(nxy);
        for (int iz = 0; iz < nz; iz++)
            Pgrid->zpiece_to_all(zpiece.data(), iz, data);
    }
#else
    for (int ixy = 0; ixy < nxy; ixy++)
        for (int iz = 0; iz < nz; iz++)
            ifs >> data[iz * nxy + ixy];
#endif
}

void ModuleIO::read_cube_core_mismatch(
    std::ifstream &ifs,
#ifdef __MPI
    const Parallel_Grid*const Pgrid,
    const bool flag_read_rank,
#endif
    double*const data,
    const int nx,
    const int ny,
    const int nz,
    const int nx_read,
    const int ny_read,
    const int nz_read)
{
#ifdef __MPI
    const int nxy = nx * ny;
    if (flag_read_rank)
    {
        std::vector<std::vector<double>> read_rho(nz, std::vector<double>(nxy));
        ModuleIO::trilinear_interpolate(ifs, nx_read, ny_read, nz_read, nx, ny, nz, read_rho);
        for (int iz = 0; iz < nz; iz++)
            Pgrid->zpiece_to_all(read_rho[iz].data(), iz, data);
    }
    else
    {
        std::vector<double> zpiece(nxy);
        for (int iz = 0; iz < nz; iz++)
            Pgrid->zpiece_to_all(zpiece.data(), iz, data);
    }
#else
    ModuleIO::trilinear_interpolate(ifs, nx_read, ny_read, nz_read, nx, ny, nz, data);
#endif
}

void ModuleIO::trilinear_interpolate(std::ifstream& ifs,
                                     const int& nx_read,
                                     const int& ny_read,
                                     const int& nz_read,
                                     const int& nx,
                                     const int& ny,
                                     const int& nz,
#ifdef __MPI
                                     std::vector<std::vector<double>> &data
#else
                                     double* data
#endif
)
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
                ifs >> read_rho[iz][ix * ny_read + iy];
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

#ifdef __MPI
                data[iz][ix * ny + iy] = result;
#else
                data[iz * nx * ny + ix * ny + iy] = result;
#endif
            }
        }
    }

    for (int iz = 0; iz < nz_read; iz++)
    {
        delete[] read_rho[iz];
    }
    delete[] read_rho;
}