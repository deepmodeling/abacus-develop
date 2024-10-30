#include "module_base/element_name.h"
#include "module_io/cube_io.h"
#include "module_parameter/parameter.h"
#include<vector>
#include "module_hamilt_pw/hamilt_pwdft/parallel_grid.h"

void ModuleIO::write_grid(
    const Parallel_Grid& pgrid,
    const double* const data,
    const int is,
    const int nspin,
    const int iter,
    const std::string& fn,
    const double ef,
    const UnitCell*const ucell,
    const int precision,
    const int out_fermi)
{
    ModuleBase::TITLE("ModuleIO", "write_grid");

    const int my_rank = GlobalV::MY_RANK;
    const int my_pool = GlobalV::MY_POOL;

    time_t start;
    time_t end;
    std::stringstream ss;

    const int& nx = pgrid.nx;
    const int& ny = pgrid.ny;
    const int& nz = pgrid.nz;
    const int& nxyz = nx * ny * nz;

    start = time(NULL);

    // reduce
    std::vector<double> data_xyz_full(nxyz);    // data to be written
#ifdef __MPI    // reduce to rank 0
    if (my_pool == 0)
    {
        pgrid.reduce(data_xyz_full.data(), data);
    }
    MPI_Barrier(MPI_COMM_WORLD);
#else
    std::memcpy(data_xyz_full.data(), data, nxyz * sizeof(double));
#endif

    // build the info structure
    if (my_rank == 0)
    {
        /// output header for cube file
        ss << "STEP: " << iter << "  Cubefile created from ABACUS. Inner loop is z, followed by y and x" << std::endl;

        ss << nspin << " (nspin) ";
        ss << std::fixed;
        ss << std::setprecision(6);
        if (out_fermi == 1)
        {
            if (PARAM.globalv.two_fermi)
            {
                if (is == 0)
                {
                    ss << ef << " (fermi energy for spin=1, in Ry)" << std::endl;
                }
                else if (is == 1)
                {
                    ss << ef << " (fermi energy for spin=2, in Ry)" << std::endl;
                }
            }
            else
            {
                ss << ef << " (fermi energy, in Ry)" << std::endl;
            }
        }
        else
        {
            ss << std::endl;
        }

        std::vector<std::string> comment(2);
        for (int i = 0;i < 2;++i) { std::getline(ss, comment[i]); }

        double fac = ucell->lat0;
        std::vector<std::vector<double>> axis_vecs =
        {
            {fac * ucell->latvec.e11 / double(nx), fac * ucell->latvec.e12 / double(nx), fac * ucell->latvec.e13 / double(nx)},
            {fac * ucell->latvec.e21 / double(ny), fac * ucell->latvec.e22 / double(ny), fac * ucell->latvec.e23 / double(ny)},
            {fac * ucell->latvec.e31 / double(nz), fac * ucell->latvec.e32 / double(nz), fac * ucell->latvec.e33 / double(nz)}
        };

        std::string element = "";
        std::vector<int> atom_type;
        std::vector<double> atom_charge;
        std::vector<std::vector<double>> atom_pos;
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
                atom_type.push_back(z);
                atom_charge.push_back(ucell->atoms[it].ncpp.zv);
                atom_pos.push_back({ fac * ucell->atoms[it].tau[ia].x, fac * ucell->atoms[it].tau[ia].y, fac * ucell->atoms[it].tau[ia].z });
            }
        }
        write_cube(fn, comment, ucell->nat, { 0.0, 0.0, 0.0 }, { nx, ny, nz }, axis_vecs, atom_type, atom_charge, atom_pos, data_xyz_full, precision);
        end = time(NULL);
        ModuleBase::GlobalFunc::OUT_TIME("write_grid", start, end);
    }

    return;
}

void ModuleIO::write_cube(const std::string& file,
    const std::vector<std::string>& comment,
    const int natom,
    const std::vector<double>& cel_pos,
    const std::vector<int>& nvoxel,
    const std::vector<std::vector<double>>& axis_vecs,
    const std::vector<int>& atom_type,
    const std::vector<double>& atom_charge,
    const std::vector<std::vector<double>>& atom_pos,
    const std::vector<double>& data,
    const int precision,
    const int ndata_line)
{
    std::ofstream ofs(file);

    for (int i = 0;i < 2;++i) { ofs << comment[i] << "\n"; }

    ofs << std::fixed;
    ofs << std::setprecision(1);    // as before
    ofs << natom << " " << cel_pos[0] << " " << cel_pos[1] << " " << cel_pos[2] << " \n";

    ofs << std::setprecision(6);    //as before
    for (int i = 0;i < 3;++i)
    {
        ofs << nvoxel[i] << " " << axis_vecs[i][0] << " " << axis_vecs[i][1] << " " << axis_vecs[i][2] << "\n";
    }

    for (int i = 0;i < natom;++i)
    {
        ofs << " " << atom_type[i] << " " << atom_charge[i] << " " << atom_pos[i][0] << " " << atom_pos[i][1] << " " << atom_pos[i][2] << "\n";
    }

    ofs.unsetf(std::ofstream::fixed);
    ofs << std::setprecision(precision);
    ofs << std::scientific;
    const int nxy = nvoxel[0] * nvoxel[1];
    const int nz = nvoxel[2];
    for (int ixy = 0; ixy < nxy; ++ixy)
    {
        for (int iz = 0;iz < nz;++iz)
        {
            ofs << " " << data[ixy * nz + iz];
            if ((iz + 1) % ndata_line == 0 && iz != nz - 1) { ofs << "\n"; }
        }
        ofs << "\n";
    }
    ofs.close();
}