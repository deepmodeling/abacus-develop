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
    const int ir = pv.para_2d().irloc(iw1_all);
    const int ic = pv.para_2d().icloc(iw2_all);
    const long index = ir * pv.para_2d().ncloc() + ic;
    
    if( index >= pv.para_2d().size())
    {
        std::cout << " iw1_all = " << iw1_all << std::endl;
        std::cout << " iw2_all = " << iw2_all << std::endl;
        std::cout << " ir = " << ir << std::endl;
        std::cout << " ic = " << ic << std::endl;
        std::cout << " index = " << index << std::endl;
        std::cout << " pv.nloc = " << pv.para_2d().size() << std::endl;
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
