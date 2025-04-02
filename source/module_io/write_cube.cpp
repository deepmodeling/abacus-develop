#include "module_base/element_name.h"
#include "module_io/cube_io.h"
#include "module_parameter/parameter.h"
#include<vector>
#include "module_hamilt_pw/hamilt_pwdft/parallel_grid.h"

void ModuleIO::write_vdata_palgrid(
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
    ModuleBase::TITLE("ModuleIO", "write_vdata_palgrid");

    const int my_rank = GlobalV::MY_RANK;
    const int my_pool = GlobalV::MY_POOL;

    time_t start;
    time_t end;
    std::stringstream ss;

    const int& nx = pgrid.nx;
    const int& ny = pgrid.ny;
    const int& nz = pgrid.nz;
    const int& nxyz = nx * ny * nz;

    start = time(nullptr);

    // reduce
    std::vector<double> data_xyz_full(nxyz);    // data to be written
#ifdef __MPI    // reduce to rank 0
    if (my_pool == 0 && GlobalV::MY_BNDGROUP == 0)
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
        std::vector<double> dx = { fac * ucell->latvec.e11 / double(nx), fac * ucell->latvec.e12 / double(nx), fac * ucell->latvec.e13 / double(nx) };
        std::vector<double> dy = { fac * ucell->latvec.e21 / double(ny), fac * ucell->latvec.e22 / double(ny), fac * ucell->latvec.e23 / double(ny) };
        std::vector<double> dz = { fac * ucell->latvec.e31 / double(nz), fac * ucell->latvec.e32 / double(nz), fac * ucell->latvec.e33 / double(nz) };

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
        write_cube(fn, comment, ucell->nat, { 0.0, 0.0, 0.0 }, nx, ny, nz, dx, dy, dz, atom_type, atom_charge, atom_pos, data_xyz_full, precision);
        end = time(nullptr);
        ModuleBase::GlobalFunc::OUT_TIME("write_vdata_palgrid", start, end);
    }

    return;
}

void ModuleIO::write_cube(const std::string& file,
    const std::vector<std::string>& comment,
    const int& natom,
    const std::vector<double>& origin,
    const int& nx,
    const int& ny,
    const int& nz,
    const std::vector<double>& dx,
    const std::vector<double>& dy,
    const std::vector<double>& dz,
    const std::vector<int>& atom_type,
    const std::vector<double>& atom_charge,
    const std::vector<std::vector<double>>& atom_pos,
    const std::vector<double>& data,
    const int precision,
    const int ndata_line)
{
    assert(comment.size() >= 2);
    for (int i = 0;i < 2;++i) { assert(comment[i].find("\n") == std::string::npos); }
    assert(origin.size() >= 3);
    assert(dx.size() >= 3);
    assert(dy.size() >= 3);
    assert(dz.size() >= 3);
    assert(atom_type.size() >= natom);
    assert(atom_charge.size() >= natom);
    assert(atom_pos.size() >= natom);
    for (int i = 0;i < natom;++i) { assert(atom_pos[i].size() >= 3); }
    assert(data.size() >= nx * ny * nz);

    std::ofstream ofs(file);

    for (int i = 0;i < 2;++i) { ofs << comment[i] << "\n"; }

    ofs << std::fixed;
    ofs << std::setprecision(1);    // as before

    ofs << natom << " " << origin[0] << " " << origin[1] << " " << origin[2] << " \n";

    ofs << std::setprecision(6);    //as before
    ofs << nx << " " << dx[0] << " " << dx[1] << " " << dx[2] << "\n";
    ofs << ny << " " << dy[0] << " " << dy[1] << " " << dy[2] << "\n";
    ofs << nz << " " << dz[0] << " " << dz[1] << " " << dz[2] << "\n";

    for (int i = 0;i < natom;++i)
    {
        ofs << " " << atom_type[i] << " " << atom_charge[i] << " " << atom_pos[i][0] << " " << atom_pos[i][1] << " " << atom_pos[i][2] << "\n";
    }

    ofs.unsetf(std::ofstream::fixed);
    ofs << std::setprecision(precision);
    ofs << std::scientific;
    const int nxy = nx * ny;
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