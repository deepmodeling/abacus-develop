//=======================
// AUTHOR : Peize Lin
// DATE :   2023-02-23
//=======================

#include "Matrix_Orbs22.h"

#include "module_base/timer.h"
#include "module_base/tool_title.h"
#include "module_hamilt_pw/hamilt_pwdft/global.h"

void Matrix_Orbs22::init(const int mode, 
                        const UnitCell& ucell,
                        const LCAO_Orbitals& orb, 
                        const double kmesh_times, 
                        const double rmesh_times)
{
    ModuleBase::TITLE("Matrix_Orbs22", "init");
    ModuleBase::timer::tick("Matrix_Orbs22", "init");
    int Lmax_used, Lmax;

    this->lat0   = &ucell.lat0;
    const int ntype = orb.get_ntype();
    int lmax_orb = -1, lmax_beta = -1;
    for (int it = 0; it < ntype; it++)
    {
        lmax_orb = std::max(lmax_orb, orb.Phi[it].getLmax());
        lmax_beta = std::max(lmax_beta, ucell.infoNL.Beta[it].getLmax());
    }
    const double dr = orb.get_dR();
    const double dk = orb.get_dk();
    const int kmesh = orb.get_kmesh() * kmesh_times + 1;
    int Rmesh = static_cast<int>(orb.get_Rmax() * rmesh_times / dr) + 4;
    Rmesh += 1 - Rmesh % 2;

    Center2_Orb::init_Table_Spherical_Bessel(4,
                                             mode,
                                             Lmax_used,
                                             Lmax,
                                             GlobalC::exx_info.info_ri.abfs_Lmax,
                                             lmax_orb,
                                             lmax_beta,
                                             dr,
                                             dk,
                                             kmesh,
                                             Rmesh,
                                             psb_);

    //=========================================
    // (3) make Gaunt coefficients table
    //=========================================
    this->MGT.init_Gaunt_CH(2 * Lmax + 1); // why +1
    this->MGT.init_Gaunt(2 * Lmax + 1);

    ModuleBase::timer::tick("Matrix_Orbs22", "init");
    std::cout << "Matrix_Orbs22::init()::done" << std::endl;
}

void Matrix_Orbs22::init_radial(const std::vector<std::vector<std::vector<Numerical_Orbital_Lm>>>& orb_A1,
                                const std::vector<std::vector<std::vector<Numerical_Orbital_Lm>>>& orb_A2,
                                const std::vector<std::vector<std::vector<Numerical_Orbital_Lm>>>& orb_B1,
                                const std::vector<std::vector<std::vector<Numerical_Orbital_Lm>>>& orb_B2)
{
    ModuleBase::TITLE("Matrix_Orbs22", "init_radial");
    ModuleBase::timer::tick("Matrix_Orbs22", "init_radial");
    assert(orb_A1.size() == orb_A2.size());
    assert(orb_B1.size() == orb_B2.size());
    for (size_t TA = 0; TA != orb_A1.size(); ++TA)
        for (size_t TB = 0; TB != orb_B1.size(); ++TB)
            for (int LA1 = 0; LA1 != orb_A1[TA].size(); ++LA1)
                for (size_t NA1 = 0; NA1 != orb_A1[TA][LA1].size(); ++NA1)
                    for (int LA2 = 0; LA2 != orb_A2[TA].size(); ++LA2)
                        for (size_t NA2 = 0; NA2 != orb_A2[TA][LA2].size(); ++NA2)
                            for (int LB1 = 0; LB1 != orb_B1[TB].size(); ++LB1)
                                for (size_t NB1 = 0; NB1 != orb_B1[TB][LB1].size(); ++NB1)
                                    for (int LB2 = 0; LB2 != orb_B2[TB].size(); ++LB2)
                                        for (size_t NB2 = 0; NB2 != orb_B2[TB][LB2].size(); ++NB2)
                                            center2_orb22_s[TA][TB][LA1][NA1][LA2][NA2][LB1][NB1][LB2].insert(
                                                std::make_pair(NB2,
                                                               Center2_Orb::Orb22(orb_A1[TA][LA1][NA1],
                                                                                  orb_A2[TA][LA2][NA2],
                                                                                  orb_B1[TB][LB1][NB1],
                                                                                  orb_B2[TB][LB2][NB2],
                                                                                  psb_,
                                                                                  this->MGT)));
    ModuleBase::timer::tick("Matrix_Orbs22", "init_radial");
}

void Matrix_Orbs22::init_radial(const LCAO_Orbitals& orb_A1,
                                const LCAO_Orbitals& orb_A2,
                                const LCAO_Orbitals& orb_B1,
                                const LCAO_Orbitals& orb_B2)
{
    ModuleBase::TITLE("Matrix_Orbs22", "init_radial");
    ModuleBase::timer::tick("Matrix_Orbs22", "init_radial");
    assert(orb_A1.get_ntype() == orb_A2.get_ntype());
    assert(orb_B1.get_ntype() == orb_B2.get_ntype());
    for (size_t TA = 0; TA != orb_A1.get_ntype(); ++TA)
        for (size_t TB = 0; TB != orb_B1.get_ntype(); ++TB)
            for (int LA1 = 0; LA1 != orb_A1.Phi[TA].getLmax(); ++LA1)
                for (size_t NA1 = 0; NA1 != orb_A1.Phi[TA].getNchi(LA1); ++NA1)
                    for (int LA2 = 0; LA2 <= orb_A2.Phi[TA].getLmax(); ++LA2)
                        for (size_t NA2 = 0; NA2 != orb_A2.Phi[TA].getNchi(LA2); ++NA2)
                            for (int LB1 = 0; LB1 <= orb_B1.Phi[TB].getLmax(); ++LB1)
                                for (size_t NB1 = 0; NB1 != orb_B1.Phi[TB].getNchi(LB1); ++NB1)
                                    for (int LB2 = 0; LB2 <= orb_B2.Phi[TB].getLmax(); ++LB2)
                                        for (size_t NB2 = 0; NB2 != orb_B2.Phi[TB].getNchi(LB2); ++NB2)
                                            center2_orb22_s[TA][TB][LA1][NA1][LA2][NA2][LB1][NB1][LB2].insert(
                                                std::make_pair(NB2,
                                                               Center2_Orb::Orb22(orb_A1.Phi[TA].PhiLN(LA1, NA1),
                                                                                  orb_A2.Phi[TA].PhiLN(LA2, NA2),
                                                                                  orb_B1.Phi[TB].PhiLN(LB1, NB1),
                                                                                  orb_B2.Phi[TB].PhiLN(LB2, NB2),
                                                                                  psb_,
                                                                                  this->MGT)));
    ModuleBase::timer::tick("Matrix_Orbs22", "init_radial");
}

void Matrix_Orbs22::init_radial_table()
{
    ModuleBase::TITLE("Matrix_Orbs22", "init_radial");
    ModuleBase::timer::tick("Matrix_Orbs22", "init_radial_table");
    for (auto& coA: center2_orb22_s)
        for (auto& coB: coA.second)
            for (auto& coC: coB.second)
                for (auto& coD: coC.second)
                    for (auto& coE: coD.second)
                        for (auto& coF: coE.second)
                            for (auto& coG: coF.second)
                                for (auto& coH: coG.second)
                                    for (auto& coI: coH.second)
                                        for (auto& coJ: coI.second)
                                            coJ.second.init_radial_table();
    ModuleBase::timer::tick("Matrix_Orbs22", "init_radial_table");
}

void Matrix_Orbs22::init_radial_table(const std::map<size_t, std::map<size_t, std::set<double>>>& Rs)
{
    ModuleBase::TITLE("Matrix_Orbs22", "init_radial_table_Rs");
    ModuleBase::timer::tick("Matrix_Orbs22", "init_radial_table");
    const double lat0 = *this->lat0;
    for (const auto& RsA: Rs)
        for (const auto& RsB: RsA.second)
        {
            if (auto* const center2_orb22_sAB = static_cast<std::map<
                    int,
                    std::map<
                        size_t,
                        std::map<
                            int,
                            std::map<
                                size_t,
                                std::map<
                                    int,
                                    std::map<size_t, std::map<int, std::map<size_t, Center2_Orb::Orb22>>>>>>>>* const>(
                    ModuleBase::GlobalFunc::MAP_EXIST(center2_orb22_s, RsA.first, RsB.first)))
            {
                std::set<size_t> radials;
                for (const double& R: RsB.second)
                {
                    const double position = R * lat0 / lcao_dr_;
                    const size_t iq = static_cast<size_t>(position);
                    for (size_t i = 0; i != 4; ++i)
                        radials.insert(iq + i);
                }
                for (auto& coC: *center2_orb22_sAB)
                    for (auto& coD: coC.second)
                        for (auto& coE: coD.second)
                            for (auto& coF: coE.second)
                                for (auto& coG: coF.second)
                                    for (auto& coH: coG.second)
                                        for (auto& coI: coH.second)
                                            for (auto& coJ: coI.second)
                                                coJ.second.init_radial_table();
            }
        }
    ModuleBase::timer::tick("Matrix_Orbs22", "init_radial_table");
}
