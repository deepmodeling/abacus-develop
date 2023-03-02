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
    double dipole_elec[3] = {0.0, 0.0, 0.0};

    for (int ir = 0; ir < GlobalC::rhopw->nrxx; ++ir)
    {
        int i = ir / (GlobalC::rhopw->ny * GlobalC::rhopw->nplane);
        int j = ir / GlobalC::rhopw->nplane - i * GlobalC::rhopw->ny;
        int k = ir % GlobalC::rhopw->nplane + GlobalC::rhopw->startz_current;
        double x = (double)i / GlobalC::rhopw->nx;
        double y = (double)j / GlobalC::rhopw->ny;
        double z = (double)k / GlobalC::rhopw->nz;

        dipole_elec[0] += rho_save[ir] * x * GlobalC::ucell.lat0*ModuleBase::BOHR_TO_A;
        dipole_elec[1] += rho_save[ir] * y * GlobalC::ucell.lat0*ModuleBase::BOHR_TO_A;
        dipole_elec[2] += rho_save[ir] * z * GlobalC::ucell.lat0*ModuleBase::BOHR_TO_A;
    }
    dipole_elec[0] *= GlobalC::ucell.omega / static_cast<double>( GlobalC::rhopw->nxyz );
    dipole_elec[1] *= GlobalC::ucell.omega / static_cast<double>( GlobalC::rhopw->nxyz );
	dipole_elec[2] *= GlobalC::ucell.omega / static_cast<double>( GlobalC::rhopw->nxyz );
    Parallel_Reduce::reduce_double_pool(dipole_elec[0]);
    Parallel_Reduce::reduce_double_pool(dipole_elec[1]);
    Parallel_Reduce::reduce_double_pool(dipole_elec[2]);
    double dipole_elec_x_tmp = dipole_elec[0];
	double dipole_elec_y_tmp = dipole_elec[1];
	double dipole_elec_z_tmp = dipole_elec[2];
    dipole_elec[0] = GlobalC::ucell.a1[0]*dipole_elec_x_tmp + GlobalC::ucell.a2[0]*dipole_elec_y_tmp + GlobalC::ucell.a3[0]*dipole_elec_z_tmp;
	dipole_elec[1] = GlobalC::ucell.a1[1]*dipole_elec_x_tmp + GlobalC::ucell.a2[1]*dipole_elec_y_tmp + GlobalC::ucell.a3[1]*dipole_elec_z_tmp;
	dipole_elec[2] = GlobalC::ucell.a1[2]*dipole_elec_x_tmp + GlobalC::ucell.a2[2]*dipole_elec_y_tmp + GlobalC::ucell.a3[2]*dipole_elec_z_tmp;
    
    std::cout << std::setprecision(8) << "dipole_elec_x: " << dipole_elec[0] << std::endl;
    std::cout << std::setprecision(8) << "dipole_elec_y: " << dipole_elec[1] << std::endl;
    std::cout << std::setprecision(8) << "dipole_elec_z: " << dipole_elec[2] << std::endl;

    ofs << std::setprecision(8) << istep << " " << dipole_elec[0] << " " << dipole_elec[1] << " " << dipole_elec[2]
        << std::endl;

    double dipole_ion_x = 0.0, dipole_ion_y = 0.0, dipole_ion_z = 0.0, dipole_sum = 0.0;
    if (GlobalC::ucell.ntype == 1)
    {
        for (int ia = 0; ia < GlobalC::ucell.atoms[0].na; ia++)
        {
            dipole_ion_x
                += GlobalC::ucell.atoms[0].tau[ia].x * GlobalC::ucell.lat0 * 0.529177 * INPUT.td_val_elec_01;
            dipole_ion_y
                += GlobalC::ucell.atoms[0].tau[ia].y * GlobalC::ucell.lat0 * 0.529177 * INPUT.td_val_elec_01;
            dipole_ion_z
                += GlobalC::ucell.atoms[0].tau[ia].z * GlobalC::ucell.lat0 * 0.529177 * INPUT.td_val_elec_01;
        }
    }
    else if (GlobalC::ucell.ntype == 2)
    {
        for (int ia = 0; ia < GlobalC::ucell.atoms[0].na; ia++)
        {
            dipole_ion_x
                += GlobalC::ucell.atoms[0].tau[ia].x * GlobalC::ucell.lat0 * 0.529177 *INPUT.td_val_elec_01;
            dipole_ion_y
                += GlobalC::ucell.atoms[0].tau[ia].y * GlobalC::ucell.lat0 * 0.529177 * INPUT.td_val_elec_01;
            dipole_ion_z
                += GlobalC::ucell.atoms[0].tau[ia].z * GlobalC::ucell.lat0 * 0.529177 * INPUT.td_val_elec_01;
        }
        for (int ia = 0; ia < GlobalC::ucell.atoms[1].na; ia++)
        {
            dipole_ion_x
                += GlobalC::ucell.atoms[1].tau[ia].x * GlobalC::ucell.lat0 * 0.529177 * INPUT.td_val_elec_02;
            dipole_ion_y
                += GlobalC::ucell.atoms[1].tau[ia].y * GlobalC::ucell.lat0 * 0.529177 * INPUT.td_val_elec_02;
            dipole_ion_z
                += GlobalC::ucell.atoms[1].tau[ia].z * GlobalC::ucell.lat0 * 0.529177 * INPUT.td_val_elec_02;
        }
    }
    else if (GlobalC::ucell.ntype == 3)
    {
        for (int ia = 0; ia < GlobalC::ucell.atoms[0].na; ia++)
        {
            dipole_ion_x
                += GlobalC::ucell.atoms[0].tau[ia].x * GlobalC::ucell.lat0 * 0.529177 * INPUT.td_val_elec_01;
            dipole_ion_y
                += GlobalC::ucell.atoms[0].tau[ia].y * GlobalC::ucell.lat0 * 0.529177 * INPUT.td_val_elec_01;
            dipole_ion_z
                += GlobalC::ucell.atoms[0].tau[ia].z * GlobalC::ucell.lat0 * 0.529177 * INPUT.td_val_elec_01;
        }
        for (int ia = 0; ia < GlobalC::ucell.atoms[1].na; ia++)
        {
            dipole_ion_x
                += GlobalC::ucell.atoms[1].tau[ia].x * GlobalC::ucell.lat0 * 0.529177 * INPUT.td_val_elec_02;
            dipole_ion_y
                += GlobalC::ucell.atoms[1].tau[ia].y * GlobalC::ucell.lat0 * 0.529177 * INPUT.td_val_elec_02;
            dipole_ion_z
                += GlobalC::ucell.atoms[1].tau[ia].z * GlobalC::ucell.lat0 * 0.529177 * INPUT.td_val_elec_02;
        }
        for (int ia = 0; ia < GlobalC::ucell.atoms[2].na; ia++)
        {
            dipole_ion_x
                += GlobalC::ucell.atoms[2].tau[ia].x * GlobalC::ucell.lat0 * 0.529177 * INPUT.td_val_elec_03;
            dipole_ion_y
                += GlobalC::ucell.atoms[2].tau[ia].y * GlobalC::ucell.lat0 * 0.529177 * INPUT.td_val_elec_03;
            dipole_ion_z
                += GlobalC::ucell.atoms[2].tau[ia].z * GlobalC::ucell.lat0 * 0.529177 * INPUT.td_val_elec_03;
        }
    }
    else
    {
        std::cout << "atom ntype is too large!" << std::endl;
    }

    std::cout << std::setprecision(8) << "dipole_ion_x: " << dipole_ion_x << std::endl;
    std::cout << std::setprecision(8) << "dipole_ion_y: " << dipole_ion_y << std::endl;
    std::cout << std::setprecision(8) << "dipole_ion_z: " << dipole_ion_z << std::endl;

    double dipole_x = 0.0, dipole_y = 0.0, dipole_z = 0.0;
    dipole_x = dipole_ion_x - dipole_elec[0];
    dipole_y = dipole_ion_y - dipole_elec[1];
    dipole_z = dipole_ion_z - dipole_elec[2];
    std::cout << std::setprecision(8) << "dipole_x: " << dipole_x << std::endl;
    std::cout << std::setprecision(8) << "dipole_y: " << dipole_y << std::endl;
    std::cout << std::setprecision(8) << "dipole_z: " << dipole_z << std::endl;
    dipole_sum = sqrt(dipole_x * dipole_x + dipole_y * dipole_y + dipole_z * dipole_z);
    std::cout << std::setprecision(8) << "dipole_sum: " << dipole_sum << std::endl;

#endif

    if (GlobalV::MY_RANK == 0)
    {
        end = time(NULL);
        ModuleBase::GlobalFunc::OUT_TIME("write_dipole", start, end);
        ofs.close();
    }

    return;
}
