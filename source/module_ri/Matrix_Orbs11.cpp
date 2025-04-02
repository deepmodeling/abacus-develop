//=======================
// AUTHOR : Peize Lin
// DATE :   2022-08-17
//=======================

#include "Matrix_Orbs11.h"

#include "module_base/timer.h"
#include "module_base/tool_title.h"
#include "module_hamilt_pw/hamilt_pwdft/global.h"

void Matrix_Orbs11::init(const int mode, 
                         const UnitCell& ucell,
                         const LCAO_Orbitals& orb, 
                         const double kmesh_times, 
                         const double rmesh_times)
{
    ModuleBase::TITLE("Matrix_Orbs11", "init");
    ModuleBase::timer::tick("Matrix_Orbs11", "init");

    int Lmax_used, Lmax;
    this->lat0 = &ucell.lat0;
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

    Center2_Orb::init_Table_Spherical_Bessel(2,
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
    this->MGT.init_Gaunt_CH(Lmax);
    this->MGT.init_Gaunt(Lmax);

    ModuleBase::timer::tick("Matrix_Orbs11", "init");
}

void Matrix_Orbs11::init_radial(const std::vector<std::vector<std::vector<Numerical_Orbital_Lm>>>& orb_A,
                                const std::vector<std::vector<std::vector<Numerical_Orbital_Lm>>>& orb_B)
{
    ModuleBase::TITLE("Matrix_Orbs11", "init_radial");
    ModuleBase::timer::tick("Matrix_Orbs11", "init_radial");
    for (size_t TA = 0; TA != orb_A.size(); ++TA) {
        for (size_t TB = 0; TB != orb_B.size(); ++TB) {
            for (int LA = 0; LA != orb_A[TA].size(); ++LA) {
                for (size_t NA = 0; NA != orb_A[TA][LA].size(); ++NA) {
                    for (int LB = 0; LB != orb_B[TB].size(); ++LB) {
                        for (size_t NB = 0; NB != orb_B[TB][LB].size(); ++NB) {
                            center2_orb11_s[TA][TB][LA][NA][LB].insert(std::make_pair(
                                NB,
                                Center2_Orb::Orb11(orb_A[TA][LA][NA], orb_B[TB][LB][NB], psb_, this->MGT)));
                        }
                    }
                }
            }
        }
    }
    ModuleBase::timer::tick("Matrix_Orbs11", "init_radial");
}

void Matrix_Orbs11::init_radial(const LCAO_Orbitals& orb_A, const LCAO_Orbitals& orb_B)
{
    ModuleBase::TITLE("Matrix_Orbs11", "init_radial");
    ModuleBase::timer::tick("Matrix_Orbs11", "init_radial");
    for (size_t TA = 0; TA != orb_A.get_ntype(); ++TA) {
        for (size_t TB = 0; TB != orb_B.get_ntype(); ++TB) {
            for (int LA = 0; LA <= orb_A.Phi[TA].getLmax(); ++LA) {
                for (size_t NA = 0; NA != orb_A.Phi[TA].getNchi(LA); ++NA) {
                    for (int LB = 0; LB <= orb_B.Phi[TB].getLmax(); ++LB) {
                        for (size_t NB = 0; NB != orb_B.Phi[TB].getNchi(LB); ++NB) {
                            center2_orb11_s[TA][TB][LA][NA][LB].insert(
                                std::make_pair(NB,
                                               Center2_Orb::Orb11(orb_A.Phi[TA].PhiLN(LA, NA),
                                                                  orb_B.Phi[TB].PhiLN(LB, NB),
                                                                  psb_,
                                                                  this->MGT)));
                        }
                    }
                }
            }
        }
    }
    ModuleBase::timer::tick("Matrix_Orbs11", "init_radial");
}

void Matrix_Orbs11::init_radial_table()
{
    ModuleBase::TITLE("Matrix_Orbs11", "init_radial_table");
    ModuleBase::timer::tick("Matrix_Orbs11", "init_radial_table");
    for (auto& coA: center2_orb11_s) {
        for (auto& coB: coA.second) {
            for (auto& coC: coB.second) {
                for (auto& coD: coC.second) {
                    for (auto& coE: coD.second) {
                        for (auto& coF: coE.second) {
                            coF.second.init_radial_table();
                        }
                    }
                }
            }
        }
    }
    ModuleBase::timer::tick("Matrix_Orbs11", "init_radial_table");
}

void Matrix_Orbs11::init_radial_table(const std::map<size_t, std::map<size_t, std::set<double>>>& Rs)
{
    ModuleBase::TITLE("Matrix_Orbs11", "init_radial_table_Rs");
    ModuleBase::timer::tick("Matrix_Orbs11", "init_radial_table");
    const double lat0 = *this->lat0;
    for (const auto& RsA: Rs) {
        for (const auto& RsB: RsA.second)
        {
            if (auto* const center2_orb11_sAB = static_cast<
                    std::map<int, std::map<size_t, std::map<int, std::map<size_t, Center2_Orb::Orb11>>>>* const>(
                    ModuleBase::GlobalFunc::MAP_EXIST(center2_orb11_s, RsA.first, RsB.first)))
            {
                std::set<size_t> radials;
                for (const double& R: RsB.second)
                {
                    const double position = R * lat0 / lcao_dr_;
                    const size_t iq = static_cast<size_t>(position);
                    for (size_t i = 0; i != 4; ++i) {
                        radials.insert(iq + i);
                    }
                }
                for (auto& coC: *center2_orb11_sAB) {
                    for (auto& coD: coC.second) {
                        for (auto& coE: coD.second) {
                            for (auto& coF: coE.second) {
                                coF.second.init_radial_table(radials);
                            }
                        }
                    }
                }
            }
        }
}
    ModuleBase::timer::tick("Matrix_Orbs11", "init_radial_table");
}
