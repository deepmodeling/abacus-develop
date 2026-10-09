#include "source_lcao/lcao_domain.h"

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
        ModuleBase::WARNING_QUIT("LCAO_domain","set_force");
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

}
