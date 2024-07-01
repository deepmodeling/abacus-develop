#include "ORB_gen_tables.h"

#include "ORB_read.h"
#include "module_base/constants.h"
#include "module_base/math_integral.h"
#include "module_base/math_polyint.h"
#include "module_base/timer.h"
#include "module_base/ylm.h"
#include "module_hamilt_lcao/hamilt_lcaodft/center2_orb.h"

ORB_gen_tables::ORB_gen_tables()
{
}
ORB_gen_tables::~ORB_gen_tables()
{
}

/// call in hamilt_linear::init_before_ions.
void ORB_gen_tables::gen_tables(std::ofstream& ofs_in,
                                LCAO_Orbitals& orb,
                                const int& Lmax_exx,
                                const bool& deepks_setorb,
                                const int& nprojmax,
                                const int* nproj,
                                const Numerical_Nonlocal* beta_)
{
    ModuleBase::TITLE("ORB_gen_tables", "gen_tables");
    ModuleBase::timer::tick("ORB_gen_tables", "gen_tables");

    ofs_in << "\n SETUP THE TWO-CENTER INTEGRATION TABLES" << std::endl;


    // PLEASE add explanations for all options of 'orb_num' and 'mode'
    // mohan add 2021-04-03
    // Peize Lin update 2016-01-26
#ifdef __ORBITAL
    int orb_num = 4;
#else
    int orb_num = 2; //
#endif
    int mode = 1; // 1: <phi|phi> and <phi|beta>
    int Lmax_used = 0;
    int Lmax = 0;

    const int ntype = orb.get_ntype();
    int lmax_orb = -1, lmax_beta = -1;
    for (int it = 0; it < ntype; it++)
    {
        lmax_orb = std::max(lmax_orb, orb.Phi[it].getLmax());
        lmax_beta = std::max(lmax_beta, beta_[it].getLmax());
    }
    const double dr = orb.get_dR();
    const double dk = orb.get_dk();
    const int kmesh = orb.get_kmesh() * 4 + 1;
    int Rmesh = static_cast<int>(orb.get_Rmax() / dr) + 4;
    Rmesh += 1 - Rmesh % 2;

    Center2_Orb::init_Table_Spherical_Bessel(orb_num,
                                             mode,
                                             Lmax_used,
                                             Lmax,
                                             Lmax_exx,
                                             lmax_orb,
                                             lmax_beta,
                                             dr,
                                             dk,
                                             kmesh,
                                             Rmesh,
                                             psb_);


    /////////////////////////////
    /// (3) make Gaunt coefficients table
    /////////////////////////////

    const int lmax = (Lmax_used - 1) / 2;
    // MGT.init_Ylm_Gaunt(orb.get_lmax()+1, 0.0,PI,0.0,ModuleBase::TWO_PI);
    MGT.init_Gaunt_CH(lmax);
    // MGT.init_Gaunt(orb.get_lmax()+1);
    MGT.init_Gaunt(lmax);

    ModuleBase::timer::tick("ORB_gen_tables", "gen_tables");
    return;
}

double ORB_gen_tables::get_distance(const ModuleBase::Vector3<double>& R1, const ModuleBase::Vector3<double>& R2) const
{
    assert(this->lat0 > 0.0);
    ModuleBase::Vector3<double> dR = R1 - R2;
    return dR.norm() * this->lat0;
}
