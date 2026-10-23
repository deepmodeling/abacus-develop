#include "pos_op_basis.h"

#include "source_base/mathzone_add1.h"
#include "source_base/timer.h"
#include "source_base/tool_quit.h"
#include "source_cell/nonlocal_info_base.h"

PosOpBasis::PosOpBasis()
{
}

PosOpBasis::~PosOpBasis()
{
}

void PosOpBasis::setup_tables(const LCAO_Orbitals& orb)
{
    const int ntype = orb.get_ntype();
    int lmax_orb = -1;
    for (int it = 0; it < ntype; it++)
    {
        lmax_orb = std::max(lmax_orb, orb.Phi[it].getLmax());
    }
    const double dr = orb.get_dR();
    const double dk = orb.get_dk();
    const int kmesh = orb.get_kmesh() * 4 + 1;
    int Rmesh = static_cast<int>(orb.get_Rmax() / dr) + 4;
    Rmesh += 1 - Rmesh % 2;

    const int Lmax = lmax_orb + 1;
    const int Lmax_used = 2 * lmax_orb + 1;
    Center2_Orb::init_Table_Spherical_Bessel(Lmax_used, dr, dk, kmesh, Rmesh, psb_);
    ModuleBase::Ylm::set_coefficients();
    MGT.init_Gaunt_CH(Lmax);
    MGT.init_Gaunt(Lmax);
}

void PosOpBasis::build_orbs(const UnitCell& ucell, const LCAO_Orbitals& orb, bool cal_force)
{
    int orb_r_ntype = 0;
    int mat_Nr = orb.Phi[0].PhiLN(0, 0).getNr();
    int count_Nr = 0;

    orbs.resize(orb.get_ntype());
    for (int T = 0; T < orb.get_ntype(); ++T)
    {
        count_Nr = orb.Phi[T].PhiLN(0, 0).getNr();
        if (count_Nr > mat_Nr)
        {
            mat_Nr = count_Nr;
            orb_r_ntype = T;
        }

        orbs[T].resize(orb.Phi[T].getLmax() + 1);
        for (int L = 0; L <= orb.Phi[T].getLmax(); ++L)
        {
            orbs[T][L].resize(orb.Phi[T].getNchi(L));
            for (int N = 0; N < orb.Phi[T].getNchi(L); ++N)
            {
                const auto& orb_origin = orb.Phi[T].PhiLN(L, N);
                orbs[T][L][N].set_orbital_info(orb_origin.getLabel(),
                                               orb_origin.getType(),
                                               orb_origin.getL(),
                                               orb_origin.getChi(),
                                               orb_origin.getNr(),
                                               orb_origin.getRab(),
                                               orb_origin.getRadial(),
                                               Numerical_Orbital_Lm::Psi_Type::Psi,
                                               orb_origin.getPsi(),
                                               static_cast<int>(orb_origin.getNk() * 4) | 1,
                                               orb_origin.getDk(),
                                               orb_origin.getDruniform(),
                                               false,
                                               true,
                                               cal_force);
            }
        }
    }

    orb_r.set_orbital_info(orbs[orb_r_ntype][0][0].getLabel(),  // atom label
                           orb_r_ntype,                         // atom type
                           1,                                   // angular momentum L
                           1,                                   // number of orbitals of this L , just N
                           orbs[orb_r_ntype][0][0].getNr(),     // number of radial mesh
                           orbs[orb_r_ntype][0][0].getRab(),    // the mesh interval in radial mesh
                           orbs[orb_r_ntype][0][0].getRadial(), // radial mesh value(a.u.)
                           Numerical_Orbital_Lm::Psi_Type::Psi,
                           orbs[orb_r_ntype][0][0].getRadial(), // radial wave function
                           orbs[orb_r_ntype][0][0].getNk(),
                           orbs[orb_r_ntype][0][0].getDk(),
                           orbs[orb_r_ntype][0][0].getDruniform(),
                           false,
                           true,
                           cal_force);
}

void PosOpBasis::build_nonlocal_orbs(const UnitCell& ucell, const LCAO_Orbitals& orb, bool cal_force)
{
    const NonlocalInfoBase& infoNL_ = *ucell.infoNL;

    orbs_nonlocal.resize(orb.get_ntype());
    for (int T = 0; T < orb.get_ntype(); ++T)
    {
        const int nproj = infoNL_.get_nproj(T);
        orbs_nonlocal[T].resize(nproj);
        for (int ip = 0; ip < nproj; ip++)
        {
            int nr = infoNL_.get_proj_Nr(T, ip);
            double dr_uniform = 0.01;
            int nr_uniform
                = static_cast<int>((infoNL_.get_proj_radial(T, ip)[nr - 1] - infoNL_.get_proj_radial(T, ip)[0]) / dr_uniform) + 1;
            double* rad = new double[nr_uniform];
            double* rab = new double[nr_uniform];
            for (int ir = 0; ir < nr_uniform; ir++)
            {
                rad[ir] = ir * dr_uniform;
                rab[ir] = dr_uniform;
            }
            double* y2 = new double[nr];
            double* Beta_r_uniform = new double[nr_uniform];
            double* dbeta_uniform = new double[nr_uniform];
            ModuleBase::Mathzone_Add1::SplineD2(infoNL_.get_proj_radial(T, ip),
                                                infoNL_.get_proj_beta_r(T, ip),
                                                nr,
                                                0.0,
                                                0.0,
                                                y2);
            ModuleBase::Mathzone_Add1::Cubic_Spline_Interpolation(infoNL_.get_proj_radial(T, ip),
                                                                  infoNL_.get_proj_beta_r(T, ip),
                                                                  y2,
                                                                  nr,
                                                                  rad,
                                                                  nr_uniform,
                                                                  Beta_r_uniform,
                                                                  dbeta_uniform);

            // linear extrapolation at the zero point
            if (infoNL_.get_proj_radial(T, ip)[0] > 1e-10)
            {
                double slope = (infoNL_.get_proj_beta_r(T, ip)[1] - infoNL_.get_proj_beta_r(T, ip)[0])
                               / (infoNL_.get_proj_radial(T, ip)[1] - infoNL_.get_proj_radial(T, ip)[0]);
                Beta_r_uniform[0] = infoNL_.get_proj_beta_r(T, ip)[0] - slope * infoNL_.get_proj_radial(T, ip)[0];
            }

            // Here, the operation beta_r / r is performed. To avoid divergence at r=0, beta_r(0) is set to beta_r(1).
            // However, this may introduce issues, so caution is needed.
            for (int ir = 1; ir < nr_uniform; ir++)
            {
                Beta_r_uniform[ir] = Beta_r_uniform[ir] / rad[ir];
            }
            Beta_r_uniform[0] = Beta_r_uniform[1];

            orbs_nonlocal[T][ip].set_orbital_info(infoNL_.get_label(T),
                                                  infoNL_.get_type(T),
                                                  infoNL_.get_proj_L(T, ip),
                                                  1,
                                                  nr_uniform,
                                                  rab,
                                                  rad,
                                                  Numerical_Orbital_Lm::Psi_Type::Psi,
                                                  Beta_r_uniform,
                                                  static_cast<int>(infoNL_.get_proj_Nk(T, ip) * 4) | 1,
                                                  infoNL_.get_proj_dk(T, ip),
                                                  infoNL_.get_proj_dr_uniform(T, ip),
                                                  false,
                                                  true,
                                                  cal_force);

            delete[] rad;
            delete[] rab;
            delete[] y2;
            delete[] Beta_r_uniform;
            delete[] dbeta_uniform;
        }
    }
}

void PosOpBasis::build_iw_map(const UnitCell& ucell, int nlocal)
{
    int map_size = nlocal;
    int required_orbitals = 0;
    for (int it = 0; it < ucell.ntype; ++it)
    {
        required_orbitals += ucell.atoms[it].nw * ucell.atoms[it].na;
    }
    map_size = std::max(map_size, required_orbitals);

    iw2it.resize(map_size);
    iw2ia.resize(map_size);
    iw2iL.resize(map_size);
    iw2iN.resize(map_size);
    iw2im.resize(map_size);

    int iw = 0;
    for (int it = 0; it < ucell.ntype; it++)
    {
        for (int ia = 0; ia < ucell.atoms[it].na; ia++)
        {
            for (int iL = 0; iL < ucell.atoms[it].nwl + 1; iL++)
            {
                for (int iN = 0; iN < ucell.atoms[it].l_nchi[iL]; iN++)
                {
                    for (int im = 0; im < (2 * iL + 1); im++)
                    {
                        iw2it[iw] = it;
                        iw2ia[iw] = ia;
                        iw2iL[iw] = iL;
                        iw2iN[iw] = iN;
                        iw2im[iw] = im;
                        iw++;
                    }
                }
            }
        }
    }
}

void PosOpBasis::build(const UnitCell& ucell, const LCAO_Orbitals& orb, bool cal_force, int nlocal)
{
    ModuleBase::TITLE("PosOpBasis", "build");
    ModuleBase::timer::start("PosOpBasis", "build");

    setup_tables(orb);
    build_orbs(ucell, orb, cal_force);

    // build center2_orb11 and center2_orb21_r tables
    for (int TA = 0; TA < orb.get_ntype(); ++TA)
    {
        for (int TB = 0; TB < orb.get_ntype(); ++TB)
        {
            for (int LA = 0; LA <= orb.Phi[TA].getLmax(); ++LA)
            {
                for (int NA = 0; NA < orb.Phi[TA].getNchi(LA); ++NA)
                {
                    for (int LB = 0; LB <= orb.Phi[TB].getLmax(); ++LB)
                    {
                        for (int NB = 0; NB < orb.Phi[TB].getNchi(LB); ++NB)
                        {
                            center2_orb11[TA][TB][LA][NA][LB].insert(
                                std::make_pair(NB, Center2_Orb::Orb11(orbs[TA][LA][NA], orbs[TB][LB][NB], psb_, MGT)));
                        }
                    }
                }
            }
        }
    }

    for (int TA = 0; TA < orb.get_ntype(); ++TA)
    {
        for (int TB = 0; TB < orb.get_ntype(); ++TB)
        {
            for (int LA = 0; LA <= orb.Phi[TA].getLmax(); ++LA)
            {
                for (int NA = 0; NA < orb.Phi[TA].getNchi(LA); ++NA)
                {
                    for (int LB = 0; LB <= orb.Phi[TB].getLmax(); ++LB)
                    {
                        for (int NB = 0; NB < orb.Phi[TB].getNchi(LB); ++NB)
                        {
                            center2_orb21_r[TA][TB][LA][NA][LB].insert(
                                std::make_pair(NB, Center2_Orb::Orb21(orbs[TA][LA][NA], orb_r, orbs[TB][LB][NB], psb_, MGT)));
                        }
                    }
                }
            }
        }
    }

    for (auto& co1: center2_orb11)
    {
        for (auto& co2: co1.second)
        {
            for (auto& co3: co2.second)
            {
                for (auto& co4: co3.second)
                {
                    for (auto& co5: co4.second)
                    {
                        for (auto& co6: co5.second)
                        {
                            co6.second.init_radial_table();
                        }
                    }
                }
            }
        }
    }

    for (auto& co1: center2_orb21_r)
    {
        for (auto& co2: co1.second)
        {
            for (auto& co3: co2.second)
            {
                for (auto& co4: co3.second)
                {
                    for (auto& co5: co4.second)
                    {
                        for (auto& co6: co5.second)
                        {
                            co6.second.init_radial_table();
                        }
                    }
                }
            }
        }
    }

    build_iw_map(ucell, nlocal);

    ModuleBase::timer::end("PosOpBasis", "build");
}

void PosOpBasis::build_nonlocal(const UnitCell& ucell, const LCAO_Orbitals& orb, bool cal_force, int nlocal)
{
    ModuleBase::TITLE("PosOpBasis", "build_nonlocal");
    ModuleBase::timer::start("PosOpBasis", "build_nonlocal");

    setup_tables(orb);
    build_orbs(ucell, orb, cal_force);
    build_nonlocal_orbs(ucell, orb, cal_force);

    const NonlocalInfoBase& infoNL_ = *ucell.infoNL;

    // build center2_orb11_nonlocal and center2_orb21_r_nonlocal tables
    for (int TA = 0; TA < orb.get_ntype(); ++TA)
    {
        for (int TB = 0; TB < orb.get_ntype(); ++TB)
        {
            for (int LA = 0; LA <= orb.Phi[TA].getLmax(); ++LA)
            {
                for (int NA = 0; NA < orb.Phi[TA].getNchi(LA); ++NA)
                {
                    for (int ip = 0; ip < infoNL_.get_nproj(TB); ip++)
                    {
                        center2_orb11_nonlocal[TA][TB][LA][NA].insert(
                            std::make_pair(ip, Center2_Orb::Orb11(orbs[TA][LA][NA], orbs_nonlocal[TB][ip], psb_, MGT)));
                    }
                }
            }
        }
    }

    for (int TA = 0; TA < orb.get_ntype(); ++TA)
    {
        for (int TB = 0; TB < orb.get_ntype(); ++TB)
        {
            for (int LA = 0; LA <= orb.Phi[TA].getLmax(); ++LA)
            {
                for (int NA = 0; NA < orb.Phi[TA].getNchi(LA); ++NA)
                {
                    for (int ip = 0; ip < infoNL_.get_nproj(TB); ip++)
                    {
                        center2_orb21_r_nonlocal[TA][TB][LA][NA].insert(
                            std::make_pair(ip, Center2_Orb::Orb21(orbs[TA][LA][NA], orb_r, orbs_nonlocal[TB][ip], psb_, MGT)));
                    }
                }
            }
        }
    }

    for (auto& co1: center2_orb11_nonlocal)
    {
        for (auto& co2: co1.second)
        {
            for (auto& co3: co2.second)
            {
                for (auto& co4: co3.second)
                {
                    for (auto& co5: co4.second)
                    {
                        co5.second.init_radial_table();
                    }
                }
            }
        }
    }

    for (auto& co1: center2_orb21_r_nonlocal)
    {
        for (auto& co2: co1.second)
        {
            for (auto& co3: co2.second)
            {
                for (auto& co4: co3.second)
                {
                    for (auto& co5: co4.second)
                    {
                        co5.second.init_radial_table();
                    }
                }
            }
        }
    }

    build_iw_map(ucell, nlocal);

    ModuleBase::timer::end("PosOpBasis", "build_nonlocal");
}
