#include "write_pot.h"

#include "module_base/element_name.h"
#include "module_base/timer.h"
#include "module_elecstate/potentials/H_Hartree_pw.h"
#include "module_elecstate/potentials/efield.h"
#include "module_hamilt_pw/hamilt_pwdft/global.h"
#include "module_io/cube_io.h"

namespace ModuleIO
{

void write_pot(
    const int &out_pot,
    const int &nspin, 
    const std::string &global_out_dir,
#ifdef __MPI
    const int& bz,
    const int& nbz,
    const int& nplane,
    const int& startz_current,
#endif
    const int& nx,
    const int& ny,
    const int& nz,
    const ModuleBase::matrix& v)
{
    ModuleBase::TITLE("ModuleIO", "write_pot");
    if(out_pot == 3)
    {
        for(int is = 0; is < nspin; is++)
        {
            std::stringstream ss;
            ss << global_out_dir << "SPIN" << is+1 << "_POT_INI.cube";
            ModuleIO::write_pot_spin(
                    out_pot,
#ifdef __MPI
					bz,
					nbz,
					nplane,
					startz_current,
#endif
                    is,
                    0, // iter
                    ss.str(),
                    nx,
                    ny,
                    nz,
                    v,
                    11); // precsion
        }
    }

    ModuleBase::TITLE("ModuleIO", "write_pot");
    return;
}



void write_pot_spin(
    const int& out_pot,
#ifdef __MPI
    const int& bz,
    const int& nbz,
    const int& nplane,
    const int& startz_current,
#endif
    const int& is,
    const int& iter,
    const std::string& fn,
    const int& nx,
    const int& ny,
    const int& nz,
    const ModuleBase::matrix& v,
    const int& precision,
    const int& hartree)
{
    ModuleBase::TITLE("ModuleIO", "write_pot_spin");
    if (out_pot != 1 && out_pot != 3)
    {
        return;
    }
    ModuleBase::timer::tick("ModuleIO", "write_pot_spin");

    double* temp_v = nullptr;
    if (is == 0)
    {
        temp_v = v.c;
    }
    else if (is == 1)
    {
        temp_v = &(v.c[nx * ny * nz]);
    }

    double ef_tmp = 0.0;
    int out_fermi = 0;

    ModuleIO::write_cube(
#ifdef __MPI
        bz,
        nbz,
        nplane,
        startz_current,
#endif
        temp_v,
        is,
        GlobalV::NSPIN,
        iter,
        fn,
        nx,
        ny,
        nz,
        ef_tmp,
        &(GlobalC::ucell),
        precision,
        out_fermi);

    ModuleBase::timer::tick("ModuleIO", "write_pot_spin");
    return;
}

void write_elecstat_pot(
#ifdef __MPI
    const int& bz,
    const int& nbz,
#endif
    const std::string& fn,
    ModulePW::PW_Basis* rho_basis,
    const Charge* const chr,
    const UnitCell* ucell,
    const double* v_eff)
{
    ModuleBase::TITLE("ModuleIO", "write_elecstat_pot");
    ModuleBase::timer::tick("ModuleIO", "write_elecstat_pot");

    std::vector<double> v_elecstat(rho_basis->nrxx, 0.0);

    const int nspin = GlobalV::NSPIN;
    const int efield = GlobalV::EFIELD_FLAG;
    const int dip_corr = GlobalV::DIP_COR_FLAG;
    const bool imp_sol = GlobalV::imp_sol;

    //==========================================
    // Hartree potential
    //==========================================
    ModuleBase::matrix vh(nspin, rho_basis->nrxx);
    vh = elecstate::H_Hartree_pw::v_hartree(*ucell, rho_basis, nspin, chr->rho);

    //==========================================
    //! Dipole correction
    //==========================================
    ModuleBase::matrix v_efield;
    if (efield>0 && dip_corr>0)
    {
        v_efield.create(nspin, rho_basis->nrxx);
        v_efield = elecstate::Efield::add_efield(*ucell,
                                                 const_cast<ModulePW::PW_Basis*>(rho_basis),
                                                 nspin,
                                                 chr->rho,
                                                 GlobalC::solvent_model);
    }

    //==========================================
    //! Add hartree potential and local pseudopot
    //==========================================
    for (int ir = 0; ir < rho_basis->nrxx; ir++)
    {
        // the spin index is 0
        v_elecstat[ir] = vh(0, ir) + v_eff[ir];

        if (efield>0 && dip_corr>0)
        {
            v_elecstat[ir] += v_efield(0, ir);
        }
        if(imp_sol == true)
        {
            v_elecstat[ir] += GlobalC::solvent_model.delta_phi[ir];
        }
    }

    //-------------------------------------------
    //! Get the vacuum level of the system
    //-------------------------------------------
    int direction = 2, length = rho_basis->nz;
    if (efield > 0 && dip_corr > 0)
    {
        direction = elecstate::Efield::efield_dir;
    }
    else
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

        // get the direction with the largest vacuum
        direction = 0;
        if (vacuum[1] > vacuum[0])
        {
            direction = 1;
        }
        if (vacuum[2] > vacuum[direction])
        {
            direction = 2;
        }
    }

    if (direction == 0)
    {
        length = rho_basis->nx;
    }
    else if (direction == 1)
    {
        length = rho_basis->ny;
    }

    // get the average along the direction in real space
    auto average = [](ModulePW::PW_Basis* rho_basis, const int direction, double* v, double* ave, bool abs) {
        for (int ir = 0; ir < rho_basis->nrxx; ++ir)
        {
            int index = 0;
            if (direction == 0)
            {
                index = ir / (rho_basis->ny * rho_basis->nplane);
            }
            else if (direction == 1)
            {
                int i = ir / (rho_basis->ny * rho_basis->nplane);
                index = ir / rho_basis->nplane - i * rho_basis->ny;
            }
            else if (direction == 2)
            {
                index = ir % rho_basis->nplane + rho_basis->startz_current;
            }

            double value = abs ? std::fabs(v[ir]) : v[ir];

            ave[index] += value;
        }

        // determine the length along the direction
        int length = rho_basis->nz, surface = rho_basis->nx * rho_basis->ny;
        if (direction == 0)
        {
            length = rho_basis->nx;
            surface = rho_basis->ny * rho_basis->nz;
        }
        else if (direction == 1)
        {
            length = rho_basis->ny;
            surface = rho_basis->nx * rho_basis->nz;
        }

#ifdef __MPI
        MPI_Allreduce(MPI_IN_PLACE, ave, length, MPI_DOUBLE, MPI_SUM, POOL_WORLD);
#endif

        for (int i = 0; i < length; ++i)
        {
            ave[i] /= surface;
        }
    };

    // average charge density along direction
    double* totchg = new double[rho_basis->nrxx];
    for (int ir = 0; ir < rho_basis->nrxx; ++ir)
    {
        totchg[ir] = chr->rho[0][ir];
    }
    if (GlobalV::NSPIN == 2)
    {
        for (int ir = 0; ir < rho_basis->nrxx; ++ir)
        {
            totchg[ir] += chr->rho[1][ir];
        }
    }

    double* ave = new double[length];
    for (int i = 0; i < length; i++)
    {
        ave[i] = 0.0;
    }
    average(rho_basis, direction, totchg, ave, true);

    // set vacuum to be the point in space where the electronic charge density is the minimum
    // get the index corresponding to the minimum charge density
    int min_index = 0;
    double min_value = 1e9;
    double windows[7] = {0.1, 0.2, 0.3, 0.4, 0.3, 0.2, 0.1};
    for (int i = 0; i < length; i++)
    {
        double sum = 0;
        // use a sliding average to smoothen in charge density
        for (int win = 0; win < 7; win++)
        {
            int index = (i - 3 + win + length) % length;
            sum += ave[index] * windows[win];
        }

        if (sum < min_value)
        {
            min_value = sum;
            min_index = i;
        }
    }

    // average electrostatic potential along direction
    for (int i = 0; i < length; i++)
    {
        ave[i] = 0.0;
    }
    average(rho_basis, direction, v_elecstat.data(), ave, false);

    // get the vacuum level
    double vacuum = ave[min_index] * ModuleBase::Ry_to_eV;
    std::cout << " min index = " << min_index << " min value = " << min_value << std::endl;
    GlobalV::ofs_running << "The vacuum level is " << vacuum << " eV" << std::endl;

    delete[] ave;
    delete[] totchg;

    //-------------------------------------------
    //! Write down the electrostatic potential
    //-------------------------------------------
    int precision = 9;
    int is = -1;
    double ef_tmp = 0.0;
    int out_fermi = 0;

    ModuleIO::write_cube(
#ifdef __MPI
        bz,
        nbz,
        rho_basis->nplane,
        rho_basis->startz_current,
#endif
        v_elecstat.data(),
        is,
        nspin,
        0,
        fn,
        rho_basis->nx,
        rho_basis->ny,
        rho_basis->nz,
        ef_tmp,
        &(GlobalC::ucell),
        precision,
        out_fermi);

    ModuleBase::timer::tick("ModuleIO", "write_elecstat_pot");
    return;
}

} // namespace ModuleIO
