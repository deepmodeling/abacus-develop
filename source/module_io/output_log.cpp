#include "output_log.h"

#include "module_parameter/parameter.h"
#include "module_base/constants.h"
#include "module_base/formatter.h"
#include "module_base/global_variable.h"

namespace ModuleIO
{
void output_convergence_after_scf(bool& convergence, double& energy, std::ofstream& ofs_running)
{
    if (convergence)
    {
        ofs_running << "\n charge density convergence is achieved" << std::endl;
        ofs_running << " final etot is " << std::setprecision(11) << energy * ModuleBase::Ry_to_eV << " eV" << std::endl;
    }
    else
    {
        ofs_running << " !! convergence has not been achieved @_@" << std::endl;
        std::cout << " !! CONVERGENCE HAS NOT BEEN ACHIEVED !!" << std::endl;
    }
}

void output_after_relax(bool conv_ion, bool conv_esolver, std::ofstream& ofs_running)
{
    if (conv_ion && !conv_esolver)
    {
        std::cout << "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%" << std::endl;
        std::cout << "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%" << std::endl;
        std::cout << " Relaxation is converged, but the SCF is unconverged! The results are unreliable. " << std::endl;
        std::cout << " It is suggested to increase the maximum SCF step and/or perform the relaxation again."
                  << std::endl;
        std::cout << "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%" << std::endl;
        std::cout << "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%" << std::endl;
        ofs_running << "\n%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%" << std::endl;
        ofs_running << "\n%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%" << std::endl;
        ofs_running << "\n Relaxation is converged, but the SCF is unconverged! The results are unreliable.. "
                    << std::endl;
        ofs_running << "\n It is suggested to increase the maximum SCF step and/or perform the relaxation again. "
                    << std::endl;
        ofs_running << "\n%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%" << std::endl;
        ofs_running << "\n%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%" << std::endl;
    }
}

void output_efermi(bool& convergence, double& efermi, std::ofstream& ofs_running)
{
    if (convergence && PARAM.inp.out_level != "m")
    {
        ofs_running << std::setprecision(16);
        ofs_running << " EFERMI = " << std::setprecision(11) << efermi * ModuleBase::Ry_to_eV << " eV" << std::endl;
    }
}

void output_vacuum_level(const UnitCell* ucell,
                         const double* const* rho,
                         const double* v_elecstat,
                         const int& nx,
                         const int& ny,
                         const int& nz,
                         const int& nxyz,
                         const int& nrxx,
                         const int& nplane,
                         const int& startz_current,
                         std::ofstream& ofs_running)
{
    // determine the vacuum direction
    double vacuum[3] = {0.0};
    for (int dir = 0; dir < 3; dir++)
    {
        std::vector<double> pos;
        for (int it = 0; it < ucell->ntype; ++it)
        {
            for (int ia = 0; ia < ucell->atoms[it].na; ++ia)
            {
                pos.push_back(ucell->atoms[it].taud[ia][dir]);
            }
        }

        std::sort(pos.begin(), pos.end());
        for (int i = 1; i < pos.size(); i++)
        {
            vacuum[dir] = std::max(vacuum[dir], pos[i] - pos[i - 1]);
        }

        // consider the periodic boundary condition
        vacuum[dir] = std::max(vacuum[dir], pos[0] + 1 - pos[pos.size() - 1]);
    }

    // we assume that the cell is a cuboid
    // get the direction with the largest vacuum
    int direction = 2;
    vacuum[0] *= ucell->latvec.e11;
    vacuum[1] *= ucell->latvec.e22;
    vacuum[2] *= ucell->latvec.e33;
    if (vacuum[0] > vacuum[2])
    {
        direction = 0;
    }
    if (vacuum[1] > vacuum[direction])
    {
        direction = 1;
    }

    int length = nz;
    if (direction == 0)
    {
        length = nx;
    }
    else if (direction == 1)
    {
        length = ny;
    }

    // get the average along the direction in real space
    auto average = [](const int& ny,
                      const int& nxyz,
                      const int& nrxx,
                      const int& nplane,
                      const int& startz_current,
                      const int& direction,
                      const int& length,
                      const double* v,
                      double* ave,
                      bool abs) {
        for (int ir = 0; ir < nrxx; ++ir)
        {
            int index = 0;
            if (direction == 0)
            {
                index = ir / (ny * nplane);
            }
            else if (direction == 1)
            {
                int i = ir / (ny * nplane);
                index = ir / nplane - i * ny;
            }
            else if (direction == 2)
            {
                index = ir % nplane + startz_current;
            }

            double value = abs ? std::fabs(v[ir]) : v[ir];

            ave[index] += value;
        }

#ifdef __MPI
        MPI_Allreduce(MPI_IN_PLACE, ave, length, MPI_DOUBLE, MPI_SUM, POOL_WORLD);
#endif

        int surface = nxyz / length;
        for (int i = 0; i < length; ++i)
        {
            ave[i] /= surface;
        }
    };

    // average charge density along direction
    std::vector<double> totchg(nrxx, 0.0);
    for (int ir = 0; ir < nrxx; ++ir)
    {
        totchg[ir] = rho[0][ir];
    }
    if (PARAM.inp.nspin == 2)
    {
        for (int ir = 0; ir < nrxx; ++ir)
        {
            totchg[ir] += rho[1][ir];
        }
    }

    std::vector<double> ave(length, 0.0);
    average(ny, nxyz, nrxx, nplane, startz_current, direction, length, totchg.data(), ave.data(), true);

    // set vacuum to be the point in space where the electronic charge density is the minimum
    // get the index corresponding to the minimum charge density
    int min_index = 0;
    double min_value = 1e9;
    double windows[7] = {0.1, 0.2, 0.3, 0.4, 0.3, 0.2, 0.1};
    for (int i = 0; i < length; i++)
    {
        double sum = 0;
        int temp = i - 3 + length;
        // use a sliding average to smoothen in charge density
        for (int win = 0; win < 7; win++)
        {
            int index = (temp + win) % length;
            sum += ave[index] * windows[win];
        }

        if (sum < min_value)
        {
            min_value = sum;
            min_index = i;
        }
    }

    // average electrostatic potential along direction
    ave.assign(ave.size(), 0.0);
    average(ny, nxyz, nrxx, nplane, startz_current, direction, length, v_elecstat, ave.data(), false);

    // get the vacuum level
    double vacuum_level = ave[min_index] * ModuleBase::Ry_to_eV;
    ofs_running << "The vacuum level is " << vacuum_level << " eV" << std::endl;
}

void print_force(std::ofstream& ofs_running,
                 const UnitCell& cell,
                 const std::string& name,
                 const ModuleBase::matrix& force,
                 bool ry)
{
    const double output_acc = 1.0e-8;
    double fac = 1.0;
    if (!ry)
    {
        fac = ModuleBase::Ry_to_eV / 0.529177;
    }

    std::vector<std::string> atom_label;
    std::vector<double> force_x;
    std::vector<double> force_y;
    std::vector<double> force_z;
    std::string table;
    std::vector<std::string> titles({name, "", "", ""});
    int iat = 0;
    for (int it = 0; it < cell.ntype; it++)
    {
        for (int ia = 0; ia < cell.atoms[it].na; ia++)
        {
            std::string atom_labels = cell.atoms[it].label + std::to_string(ia + 1);
            double fx = std::abs(force(iat, 0)) > output_acc ? force(iat, 0) * fac : 0.0;
            double fy = std::abs(force(iat, 1)) > output_acc ? force(iat, 1) * fac : 0.0;
            double fz = std::abs(force(iat, 2)) > output_acc ? force(iat, 2) * fac : 0.0;
            atom_label.push_back(atom_labels);
            force_x.push_back(fx);
            force_y.push_back(fy);
            force_z.push_back(fz);

            iat++;
        }
    }

    FmtTable fmt(titles, atom_label.size(), {"%10s", "%20.10f", "%20.10f", "%20.10f"});
    fmt << atom_label << force_x << force_y << force_z;
    table = fmt.str();
    ofs_running << table << std::endl;
    if (PARAM.inp.test_force) { std::cout << table << std::endl;
}
}
void print_stress(const std::string& name, const ModuleBase::matrix& scs, const bool screen, const bool ry)
{
    const double output_acc = 1.0e-8;
    double unit_transform = 1;
    std::string title = name;
    std::string unit = "";
    if (ry)
    {
        title += " (a.u.)";
        unit = " a.u.";
    }
    else
    {
        title += " (KBAR)";
        unit = " KBAR";
        unit_transform = ModuleBase::RYDBERG_SI / pow(ModuleBase::BOHR_RADIUS_SI, 3) * 1.0e-8;
    }
    std::vector<double> stress_x;
    std::vector<double> stress_y;
    std::vector<double> stress_z;
    std::string table;
    std::vector<std::string> titles({title, "", ""});
    for (int i = 0; i < 3; i++)
    {
        double sx = scs(i, 0) * unit_transform;
        double sy = scs(i, 1) * unit_transform;
        double sz = scs(i, 2) * unit_transform;
        stress_x.push_back(sx);
        stress_y.push_back(sy);
        stress_z.push_back(sz);
    }

    double pressure = (scs(0, 0) + scs(1, 1) + scs(2, 2)) / 3.0 * unit_transform;

    FmtTable fmt(titles, 3, {"%20.10f", "%20.10f", "%20.10f"});
    fmt << stress_x << stress_y << stress_z;
    table = fmt.str();
    GlobalV::ofs_running << table;
    if (name == "TOTAL-STRESS")
    {
        GlobalV::ofs_running << " TOTAL-PRESSURE: " << std::fixed << std::setprecision(6) << pressure << unit
                             << std::endl
                             << std::endl;
    }
    if (screen)
    {
        std::cout << table;
        if (name == "TOTAL-STRESS")
        {
            std::cout << " TOTAL-PRESSURE: " << std::fixed << std::setprecision(6) << pressure << unit << std::endl
                      << std::endl;
        }
    }
    return;
}

void write_head(std::ofstream& ofs_running, const int& istep, const int& iter, const std::string& basisname)
{
    ofs_running << "\n " << basisname << " ALGORITHM --------------- ION=" << std::setw(4) << istep + 1
                << "  ELEC=" << std::setw(4) << iter << "--------------------------------\n";
}

}// namespace ModuleIO