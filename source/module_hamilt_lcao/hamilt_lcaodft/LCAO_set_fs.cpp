#include "module_hamilt_lcao/hamilt_lcaodft/LCAO_domain.h"

namespace LCAO_domain
{

void set_force
(
    const Parallel_Orbitals &pv,
    const int &iw1_all,
    const int &iw2_all,
    const double& vx,
    const double& vy,
    const double& vz,
    const char &dtype,
    double* dsloc_x,
    double* dsloc_y,
    double* dsloc_z,
    double* dhloc_fixed_x,
    double* dhloc_fixed_y,
    double* dhloc_fixed_z)
{
    // use iw1_all and iw2_all to set Hloc
    // becareful! The ir and ic may < 0!!!!!!!!!!!!!!!!
    const int ir = pv.global2local_row(iw1_all);
    const int ic = pv.global2local_col(iw2_all);
    const long index = ir * pv.ncol + ic;
    
    if( index >= pv.nloc)
    {
        std::cout << " iw1_all = " << iw1_all << std::endl;
        std::cout << " iw2_all = " << iw2_all << std::endl;
        std::cout << " ir = " << ir << std::endl;
        std::cout << " ic = " << ic << std::endl;
        std::cout << " index = " << index << std::endl;
        std::cout << " pv.nloc = " << pv.nloc << std::endl;
        ModuleBase::WARNING_QUIT("LCAO_Matrix","set_force");
    }	 

    if (dtype == 'S')
    {
        dsloc_x[index] += vx;
        dsloc_y[index] += vy;
        dsloc_z[index] += vz;
    }
    else if (dtype == 'T')
    {
        // notice, the sign is '-', minus.
        dhloc_fixed_x[index] -= vx;
        dhloc_fixed_y[index] -= vy;
        dhloc_fixed_z[index] -= vz;
    }
    else if (dtype == 'N')
    {
        dhloc_fixed_x[index] += vx;
        dhloc_fixed_y[index] += vy;
        dhloc_fixed_z[index] += vz;
    }

    return;
}

void set_stress
(
    const Parallel_Orbitals &pv,
    const int &iw1_all,
    const int &iw2_all,
    const double& vx,
    const double& vy,
    const double& vz,
    const char &dtype,
    const ModuleBase::Vector3<double> &dtau,
    double* dsloc_11,
    double* dsloc_12,
    double* dsloc_13,
    double* dsloc_22,
    double* dsloc_23,
    double* dsloc_33,
    double* dhloc_fixed_11,
    double* dhloc_fixed_12,
    double* dhloc_fixed_13,
    double* dhloc_fixed_22,
    double* dhloc_fixed_23,
    double* dhloc_fixed_33)
{
    // use iw1_all and iw2_all to set Hloc
    // becareful! The ir and ic may < 0!!!!!!!!!!!!!!!!
    const int ir = pv.global2local_row(iw1_all);
    const int ic = pv.global2local_col(iw2_all);
    const long index = ir * pv.ncol + ic;

    if( index >= pv.nloc)
    {
        std::cout << " iw1_all = " << iw1_all << std::endl;
        std::cout << " iw2_all = " << iw2_all << std::endl;
        std::cout << " ir = " << ir << std::endl;
        std::cout << " ic = " << ic << std::endl;
        std::cout << " index = " << index << std::endl;
        std::cout << " pv.nloc = " << pv.nloc << std::endl;
        ModuleBase::WARNING_QUIT("LCAO_Matrix","set_stress");
    }

    if (dtype == 'S')
    {
        dsloc_11[index] += vx * dtau.x;
        dsloc_12[index] += vx * dtau.y;
        dsloc_13[index] += vx * dtau.z;
        dsloc_22[index] += vy * dtau.y;
        dsloc_23[index] += vy * dtau.z;
        dsloc_33[index] += vz * dtau.z;
    }
    else if (dtype == 'T')
    {
        // notice, the sign is '-', minus.
        dhloc_fixed_11[index] -= vx * dtau.x;
        dhloc_fixed_12[index] -= vx * dtau.y;
        dhloc_fixed_13[index] -= vx * dtau.z;
        dhloc_fixed_22[index] -= vy * dtau.y;
        dhloc_fixed_23[index] -= vy * dtau.z;
        dhloc_fixed_33[index] -= vz * dtau.z;
    }
    else if (dtype == 'N')
    {
        dhloc_fixed_11[index] += vx * dtau.x;
        dhloc_fixed_12[index] += vx * dtau.y;
        dhloc_fixed_13[index] += vx * dtau.z;
        dhloc_fixed_22[index] += vy * dtau.y;
        dhloc_fixed_23[index] += vy * dtau.z;
        dhloc_fixed_33[index] += vz * dtau.z;
    }

    return;
}

}
