#include "module_base/parallel_reduce.h"
#include "module_elecstate/module_charge/charge.h"
#include "module_hamilt_lcao/module_tddft/ELEC_evolve.h"
#include "module_hamilt_pw/hamilt_pwdft/global.h"
#include "module_io/dipole_io.h"

// fuxiang add 2017-03-15
void ModuleIO::write_dipole(const double *rho_save,
                            const int &is,
                            const int &istep,
                            const std::string &fn,
                            const int &precision,
                            const bool for_plot)
{
    ModuleBase::TITLE("ModuleIO", "write_dipole");

    time_t start, end;
    std::ofstream ofs;

    if (GlobalV::MY_RANK == 0)
    {
        start = time(NULL);

        ofs.open(fn.c_str(), ofstream::app);
        if (!ofs)
        {
            ModuleBase::WARNING("ModuleIO", "Can't create Charge File!");
        }
    }

    double bmod[3];
    for (int i = 0; i < 3; i++)
    {
        bmod[i] = prepare(GlobalC::ucell, i);
    }

#ifndef __MPI
    double dipole_elec_x = 0.0, dipole_elec_y = 0.0, dipole_elec_z = 0.0;
    for (int k = 0; k < GlobalC::rhopw->nz; k++)
    {
        for (int j = 0; j < GlobalC::rhopw->ny; j++)
        {
            for (int i = 0; i < GlobalC::rhopw->nx; i++)
            {
                dipole_elec_x += rho_save[i * GlobalC::rhopw->ny * GlobalC::rhopw->nz + j * GlobalC::rhopw->nz + k] * i
                                 * GlobalC::ucell.lat0 * 0.529177 / GlobalC::rhopw->nx;
                dipole_elec_y += rho_save[i * GlobalC::rhopw->ny * GlobalC::rhopw->nz + j * GlobalC::rhopw->nz + k] * j
                                 * GlobalC::ucell.lat0 * 0.529177 / GlobalC::rhopw->ny;
                dipole_elec_z += rho_save[i * GlobalC::rhopw->ny * GlobalC::rhopw->nz + j * GlobalC::rhopw->nz + k] * k
                                 * GlobalC::ucell.lat0 * 0.529177 / GlobalC::rhopw->nz;
            }
        }
    }
    dipole_elec_x *= GlobalC::ucell.omega / static_cast<double>(GlobalC::rhopw->nxyz);
    dipole_elec_y *= GlobalC::ucell.omega / static_cast<double>(GlobalC::rhopw->nxyz);
    dipole_elec_z *= GlobalC::ucell.omega / static_cast<double>(GlobalC::rhopw->nxyz);
    Parallel_Reduce::reduce_double_pool(dipole_elec_x);
    Parallel_Reduce::reduce_double_pool(dipole_elec_y);
    Parallel_Reduce::reduce_double_pool(dipole_elec_z);

    ofs << istep << " " << dipole_elec_x << " " << dipole_elec_y << dipole_elec_z;
#else

            // case 2: > first part rho: processor 0 receive the rho
            // from other processors
            else if (GlobalV::RANK_IN_POOL == 0)
            {
                MPI_Recv(zpiece, nxy, MPI_DOUBLE, which_ip[iz], tag, POOL_WORLD, &ierror);
                // GlobalV::ofs_running << "\n Receieve First number = " << zpiece[0];
            }

            // write data
            if (GlobalV::MY_RANK == 0)
            {
                //	ofs << "\niz=" << iz;
                // mohan update 2011-03-30
                for (int iy = 0; iy < GlobalC::rhopw->ny; iy++)
                {
                    for (int ix = 0; ix < GlobalC::rhopw->nx; ix++)
                    {
                        /*
                                                if(ix<GlobalC::pw.ncx/2)
                                                    {dipole_elec_x +=
                           zpiece[ix*GlobalC::pw.ncy+iy]*ix*GlobalC::ucell.lat0*0.529177/GlobalC::pw.ncx;} else
                                                    {dipole_elec_x +=
                           zpiece[ix*GlobalC::pw.ncy+iy]*(ix-GlobalC::pw.ncx)*GlobalC::ucell.lat0*0.529177/GlobalC::pw.ncx;}
                                                if(iy<GlobalC::pw.ncy/2)
                                                    {dipole_elec_y +=
                           zpiece[ix*GlobalC::pw.ncy+iy]*iy*GlobalC::ucell.lat0*0.529177/GlobalC::pw.ncy;} else
                                                    {dipole_elec_y +=
                           zpiece[ix*GlobalC::pw.ncy+iy]*(iy-GlobalC::pw.ncy)*GlobalC::ucell.lat0*0.529177/GlobalC::pw.ncy;}
                                                if(iz<GlobalC::pw.ncz/2)
                                                    {dipole_elec_z +=
                           zpiece[ix*GlobalC::pw.ncy+iy]*iz*GlobalC::ucell.lat0*0.529177/GlobalC::pw.ncz;} else
                                                    {dipole_elec_z +=
                           zpiece[ix*GlobalC::pw.ncy+iy]*(iz-GlobalC::pw.ncz)*GlobalC::ucell.lat0*0.529177/GlobalC::pw.ncz;}
                        */
                        dipole_elec_x += zpiece[ix * GlobalC::rhopw->ny + iy] * ix * GlobalC::ucell.lat0 * 0.529177
                                         / GlobalC::rhopw->nx;
                        dipole_elec_y += zpiece[ix * GlobalC::rhopw->ny + iy] * iy * GlobalC::ucell.lat0 * 0.529177
                                         / GlobalC::rhopw->ny;
                        dipole_elec_z += zpiece[ix * GlobalC::rhopw->ny + iy] * iz * GlobalC::ucell.lat0 * 0.529177
                                         / GlobalC::rhopw->nz;
                    }
                }
            }
        } // end iz

        delete[] zpiece;

    Parallel_Reduce::reduce_double_pool(dipole_elec[0]);
    Parallel_Reduce::reduce_double_pool(dipole_elec[1]);
    Parallel_Reduce::reduce_double_pool(dipole_elec[2]);
    for (int i = 0; i < 3; ++i)
    {
        dipole_elec[i] *= GlobalC::ucell.lat0 / bmod[i] * GlobalC::ucell.omega / GlobalC::rhopw->nxyz;
    }

        std::cout << std::setprecision(8) << "dipole_elec_x: " << dipole_elec_x << std::endl;
        std::cout << std::setprecision(8) << "dipole_elec_y: " << dipole_elec_y << std::endl;
        std::cout << std::setprecision(8) << "dipole_elec_z: " << dipole_elec_z << std::endl;

        ofs << istep << " " << dipole_elec_x << " " << dipole_elec_y << " " << dipole_elec_z << std::endl;

    double dipole_ion[3] = {0.0};
    double dipole_sum = 0.0;

    for (int i = 0; i < 3; ++i)
    {
        for (int it = 0; it < GlobalC::ucell.ntype; ++it)
        {
            double sum = 0;
            for (int ia = 0; ia < GlobalC::ucell.atoms[it].na; ++ia)
            {

                sum += GlobalC::ucell.atoms[it].taud[ia][i];
            }
            dipole_ion[i] += sum * GlobalC::ucell.atoms[it].ncpp.zv;
        }
        dipole_ion[i] *= GlobalC::ucell.lat0 / bmod[i]; //* ModuleBase::FOUR_PI / GlobalC::ucell.omega;
    }

    std::cout << std::setprecision(8) << "dipole_ion_x: " << dipole_ion[0] << std::endl;
    std::cout << std::setprecision(8) << "dipole_ion_y: " << dipole_ion[1] << std::endl;
    std::cout << std::setprecision(8) << "dipole_ion_z: " << dipole_ion[2] << std::endl;

    double dipole[3] = {0.0};
    for (int i = 0; i < 3; ++i)
    {
        dipole[i] = dipole_ion[i] - dipole_elec[i];
    }
    std::cout << std::setprecision(8) << "dipole_x: " << dipole[0] << std::endl;
    std::cout << std::setprecision(8) << "dipole_y: " << dipole[1] << std::endl;
    std::cout << std::setprecision(8) << "dipole_z: " << dipole[2] << std::endl;
    dipole_sum = sqrt(dipole[0] * dipole[0] + dipole[1] * dipole[1] + dipole[2] * dipole[2]);
    std::cout << std::setprecision(8) << "dipole_sum: " << dipole_sum << std::endl;

#endif

    // calculate ion dipole;
    if (GlobalV::MY_RANK == 0)
    {
        end = time(NULL);
        ModuleBase::GlobalFunc::OUT_TIME("write_dipole", start, end);
        ofs.close();
    }

    return;
}

double ModuleIO::prepare(const UnitCell &cell, int &dir)
{
    double bvec[3] = {0.0};
    double bmod = 0.0;
    if (dir == 0)
    {
        bvec[0] = cell.G.e11;
        bvec[1] = cell.G.e12;
        bvec[2] = cell.G.e13;
    }
    else if (dir == 1)
    {
        bvec[0] = cell.G.e21;
        bvec[1] = cell.G.e22;
        bvec[2] = cell.G.e23;
    }
    else if (dir == 2)
    {
        bvec[0] = cell.G.e31;
        bvec[1] = cell.G.e32;
        bvec[2] = cell.G.e33;
    }
    else
    {
        ModuleBase::WARNING_QUIT("ModuleIO::prepare", "direction is wrong!");
    }
    bmod = sqrt(pow(bvec[0], 2) + pow(bvec[1], 2) + pow(bvec[2], 2));
    return bmod;
}