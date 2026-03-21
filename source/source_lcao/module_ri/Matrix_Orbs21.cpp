//=======================
// AUTHOR : Peize Lin
// DATE :   2022-08-17
//=======================

#include "Matrix_Orbs21.h"

#include "exx_abfs-construct_orbs.h"
#include "source_base/timer.h"
#include "source_base/tool_title.h"

void Matrix_Orbs21::init(
    const std::vector<std::vector<std::vector<Numerical_Orbital_Lm>>>& orb_A1,
    const std::vector<std::vector<std::vector<Numerical_Orbital_Lm>>>& orb_A2,
    const std::vector<std::vector<std::vector<Numerical_Orbital_Lm>>>& orb_B,
    const UnitCell& ucell,
    const LCAO_Orbitals& orb, 
    const double kmesh_times)
{
    ModuleBase::TITLE("Matrix_Orbs21", "init");
    ModuleBase::timer::tick("Matrix_Orbs21", "init");
    this->lat0 = &ucell.lat0;

    const int Lmax = std::max({
        Exx_Abfs::Construct_Orbs::get_Lmax(orb_A1) + Exx_Abfs::Construct_Orbs::get_Lmax(orb_A2),
        Exx_Abfs::Construct_Orbs::get_Lmax(orb_B) });
    const int Lmax_used = Exx_Abfs::Construct_Orbs::get_Lmax(orb_A1) + Exx_Abfs::Construct_Orbs::get_Lmax(orb_A2) + Exx_Abfs::Construct_Orbs::get_Lmax(orb_B);

    //=========================================
    // (3) make Gaunt coefficients table
    //=========================================
    if(!this->MGT)
        { this->MGT = std::make_shared<ORB_gaunt_table>(); }
    if(this->MGT->get_Lmax_Gaunt_CH() < Lmax)
        { this->MGT->init_Gaunt_CH(Lmax); }
    if(this->MGT->get_Lmax_Gaunt_Coefficients() < Lmax)
        { this->MGT->init_Gaunt(Lmax); }
    
    const double dr = orb.get_dR();
    const double dk = orb.get_dk();
    const int kmesh = orb.get_kmesh() * kmesh_times + 1;
    const double rmax
        = std::min({Exx_Abfs::Construct_Orbs::get_Rmax(orb_A1), Exx_Abfs::Construct_Orbs::get_Rmax(orb_A2)})
        + Exx_Abfs::Construct_Orbs::get_Rmax(orb_B);
   int Rmesh = static_cast<int>(rmax / dr) + 4;                            // extend Rcut, keep dR
    Rmesh += 1 - Rmesh % 2;
    Center2_Orb::init_Table_Spherical_Bessel(Lmax_used,
                                             dr,
                                             dk,
                                             kmesh,
                                             Rmesh,
                                             psb_);

    assert(orb_A1.size() == orb_A2.size());
    for (size_t TA = 0; TA != orb_A1.size(); ++TA) {
        for (size_t TB = 0; TB != orb_B.size(); ++TB) {
            for (int LA1 = 0; LA1 != orb_A1[TA].size(); ++LA1) {
                for (size_t NA1 = 0; NA1 != orb_A1[TA][LA1].size(); ++NA1) {
                    for (int LA2 = 0; LA2 != orb_A2[TA].size(); ++LA2) {
                        for (size_t NA2 = 0; NA2 != orb_A2[TA][LA2].size(); ++NA2) {
                            for (int LB = 0; LB != orb_B[TB].size(); ++LB) {
                                for (size_t NB = 0; NB != orb_B[TB][LB].size(); ++NB) {
                                    center2_orb21_s[TA][TB][LA1][NA1][LA2][NA2][LB].insert(
                                        std::make_pair(
                                            NB,
                                            Center2_Orb::Orb21(
                                                orb_A1[TA][LA1][NA1],
                                                orb_A2[TA][LA2][NA2],
                                                orb_B[TB][LB][NB],
                                                psb_,
                                                *this->MGT)));
    }}}}}}}}
    ModuleBase::timer::tick("Matrix_Orbs21", "init");
}

/*
void Matrix_Orbs21::init_radial(const std::vector<std::vector<std::vector<Numerical_Orbital_Lm>>>& orb_A1,
                                const LCAO_Orbitals& orb_A2,
                                const LCAO_Orbitals& orb_B)
{
    ModuleBase::TITLE("Matrix_Orbs21", "init_radial");
    ModuleBase::timer::tick("Matrix_Orbs21", "init_radial");
    assert(orb_A1.size() == orb_A2.get_ntype());
    for (size_t TA = 0; TA != orb_A1.size(); ++TA) 
    {
        for (size_t TB = 0; TB != orb_B.get_ntype(); ++TB) 
        {
            for (int LA1 = 0; LA1 != orb_A1[TA].size(); ++LA1) 
            {
                for (size_t NA1 = 0; NA1 != orb_A1[TA][LA1].size(); ++NA1) 
                {
                    for (int LA2 = 0; LA2 <= orb_A2.Phi[TA].getLmax(); ++LA2) 
                    {
                        for (size_t NA2 = 0; NA2 != orb_A2.Phi[TA].getNchi(LA2); ++NA2) 
                        {
                            for (int LB = 0; LB <= orb_B.Phi[TB].getLmax(); ++LB) 
                            {
                                for (size_t NB = 0; NB != orb_B.Phi[TB].getNchi(LB); ++NB) 
                                {
                                    center2_orb21_s[TA][TB][LA1][NA1][LA2][NA2][LB].insert(
                                        std::make_pair(NB,
                                                       Center2_Orb::Orb21(orb_A1[TA][LA1][NA1],
                                                                          orb_A2.Phi[TA].PhiLN(LA2, NA2),
                                                                          orb_B.Phi[TB].PhiLN(LB, NB),
                                                                          psb_,
                                                                          *this->MGT)));
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    ModuleBase::timer::tick("Matrix_Orbs21", "init_radial");
}
*/

void Matrix_Orbs21::init_radial_table()
{
    ModuleBase::TITLE("Matrix_Orbs21", "init_radial_table");
    ModuleBase::timer::tick("Matrix_Orbs21", "init_radial_table");
    for (auto& coA: center2_orb21_s) 
    {
        for (auto& coB: coA.second) 
        {
            for (auto& coC: coB.second) 
            {
                for (auto& coD: coC.second) 
                {
                    for (auto& coE: coD.second) 
                    {
                        for (auto& coF: coE.second) 
                        {
                            for (auto& coG: coF.second) 
                            {
                                for (auto& coH: coG.second) 
                                {
                                    coH.second.init_radial_table();
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    ModuleBase::timer::tick("Matrix_Orbs21", "init_radial_table");
}

void Matrix_Orbs21::init_radial_table(const std::map<size_t, std::map<size_t, std::set<double>>>& Rs)
{
    ModuleBase::TITLE("Matrix_Orbs21", "init_radial_table_Rs");
    ModuleBase::timer::tick("Matrix_Orbs21", "init_radial_table");
    const double lat0 = *this->lat0;
    for (const auto& RsA: Rs) {
        for (const auto& RsB: RsA.second)
        {
            if (auto* const center2_orb21_sAB = static_cast<std::map<
                    int,
                    std::map<
                        size_t,
                        std::map<int, std::map<size_t, std::map<int, std::map<size_t, Center2_Orb::Orb21>>>>>>* const>(
                    ModuleBase::GlobalFunc::MAP_EXIST(center2_orb21_s, RsA.first, RsB.first)))
            {
                std::set<size_t> radials;
                for (const double& R: RsB.second)
                {
                    const double position = R * lat0 / lcao_dr_;
                    const size_t iq = static_cast<size_t>(position);
                    for (size_t i = 0; i != 4; ++i) 
                    {
                        radials.insert(iq + i);
                    }
                }
                for (auto& coC: *center2_orb21_sAB) 
                {
                    for (auto& coD: coC.second) 
                    {
                        for (auto& coE: coD.second) 
                        {
                            for (auto& coF: coE.second) 
                            {
                                for (auto& coG: coF.second) 
                                {
                                    for (auto& coH: coG.second) 
                                    {
                                        coH.second.init_radial_table(radials);
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    ModuleBase::timer::tick("Matrix_Orbs21", "init_radial_table");
}
