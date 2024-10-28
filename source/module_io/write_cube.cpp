#include "module_base/element_name.h"
#include "module_io/cube_io.h"
#include "module_parameter/parameter.h"
#include<vector>
#include "module_hamilt_pw/hamilt_pwdft/global.h"   // tmp, for Pgrid

void ModuleIO::write_cube(
    const double*const data,
    const int is,
    const int nspin,
    const int iter,
    const std::string& fn,
    const int nx,
    const int ny,
    const int nz,
    const double ef,
    const UnitCell*const ucell,
    const int precision,
    const int out_fermi)
{
    ModuleBase::TITLE("ModuleIO", "write_cube");

    const int my_rank = GlobalV::MY_RANK;

    time_t start;
    time_t end;
    std::ofstream ofs_cube;

    if (my_rank == 0)
    {
        start = time(NULL);

        if (iter == 0)
        {
            ofs_cube.open(fn.c_str());
        }
        else
        {
            ofs_cube.open(fn.c_str(), std::ios::app);
        }

        if (!ofs_cube)
        {
            ModuleBase::WARNING("ModuleIO::write_cube", "Can't create Output File!");
        }

        /// output header for cube file
        ofs_cube << "STEP: " << iter << "  Cubefile created from ABACUS. Inner loop is z, followed by y and x" << std::endl;
        ofs_cube << nspin << " (nspin) ";

        ofs_cube << std::fixed;
        ofs_cube << std::setprecision(6);
        if (out_fermi == 1)
        {
            if (PARAM.globalv.two_fermi)
            {
                if (is == 0)
                {
                    ofs_cube << ef << " (fermi energy for spin=1, in Ry)" << std::endl;
                }
                else if (is == 1)
                {
                    ofs_cube << ef << " (fermi energy for spin=2, in Ry)" << std::endl;
                }
            }
            else
            {
                ofs_cube << ef << " (fermi energy, in Ry)" << std::endl;
            }
        }
        else
        {
            ofs_cube << std::endl;
        }

        ofs_cube << ucell->nat << " 0.0 0.0 0.0 " << std::endl;
        double fac = ucell->lat0;
        ofs_cube << nx << " " << fac * ucell->latvec.e11 / double(nx) << " " << fac * ucell->latvec.e12 / double(nx)
                 << " " << fac * ucell->latvec.e13 / double(nx) << std::endl;
        ofs_cube << ny << " " << fac * ucell->latvec.e21 / double(ny) << " " << fac * ucell->latvec.e22 / double(ny)
                 << " " << fac * ucell->latvec.e23 / double(ny) << std::endl;
        ofs_cube << nz << " " << fac * ucell->latvec.e31 / double(nz) << " " << fac * ucell->latvec.e32 / double(nz)
                 << " " << fac * ucell->latvec.e33 / double(nz) << std::endl;

        std::string element = "";
        for (int it = 0; it < ucell->ntype; it++)
        {
            // erase the number in label, such as Fe1.
            element = ucell->atoms[it].label;
            std::string::iterator temp = element.begin();
            while (temp != element.end())
            {
                if ((*temp >= '1') && (*temp <= '9'))
                {
                    temp = element.erase(temp);
                }
                else
                {
                    temp++;
                }
            }

            for (int ia = 0; ia < ucell->atoms[it].na; ia++)
            {
                // convert from label to atomic number
                int z = 0;
                for (int j = 0; j != ModuleBase::element_name.size(); j++)
                {
                    if (element == ModuleBase::element_name[j])
                    {
                        z = j + 1;
                        break;
                    }
                }
                ofs_cube << " " << z << " " << ucell->atoms[it].ncpp.zv << " " << fac * ucell->atoms[it].tau[ia].x
                         << " " << fac * ucell->atoms[it].tau[ia].y << " " << fac * ucell->atoms[it].tau[ia].z
                         << std::endl;
            }
        }
        ofs_cube.unsetf(std::ostream::fixed);
        ofs_cube << std::setprecision(precision);
        ofs_cube << std::scientific;
    }

    ModuleIO::write_cube_core(ofs_cube, data, nx * ny, nz, 6);

    if (my_rank == 0)
    {
        end = time(NULL);
        ModuleBase::GlobalFunc::OUT_TIME("write_cube", start, end);

        /// for cube file
        ofs_cube.close();
    }

    return;
}


void ModuleIO::write_cube_core(
    std::ofstream& ofs_cube,
    const double*const data,
    const int nxy,
    const int nz,
    const int n_data_newline)
{
    ModuleBase::TITLE("ModuleIO", "write_cube_core");

#ifdef __MPI

    const int my_rank = GlobalV::MY_RANK;
    const int my_pool = GlobalV::MY_POOL;

    // only do in the first pool.
    if (my_pool == 0)
    {
        /// for cube file
        const int nxyz = nxy * nz;
        std::vector<double> data_cube(nxyz, 0.0);

        GlobalC::Pgrid.reduce_to_fullrho(data_cube.data(), data);

        // for cube file
        if (my_rank == 0)
        {
            for (int ixy = 0; ixy < nxy; ixy++)
            {
                for (int iz = 0; iz < nz; iz++)
                {
                    ofs_cube << " " << data_cube[ixy * nz + iz];
                    if ((iz % n_data_newline == n_data_newline-1) && (iz != nz - 1))
                    {
                        ofs_cube << "\n";
                    }
                }
                ofs_cube << "\n";
            }
        }
        /// for cube file
    }
    MPI_Barrier(MPI_COMM_WORLD);
#else
    for (int ixy = 0; ixy < nxy; ixy++)
    {
        for (int iz = 0; iz < nz; iz++)
        {
            ofs_cube << " " << data[iz * nxy + ixy];
            // ++count_cube;
            if ((iz % n_data_newline == n_data_newline-1) && (iz != nz - 1))
            {
                ofs_cube << "\n";
            }
        }
        ofs_cube << "\n";
    }
#endif
}
