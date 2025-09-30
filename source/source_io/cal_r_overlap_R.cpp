#include "cal_r_overlap_R.h"

#include "source_io/module_parameter/parameter.h"
#include "source_base/parallel_reduce.h"
#include "source_base/timer.h"
#include "source_cell/module_neighbor/sltk_grid_driver.h"
#include "source_pw/module_pwdft/global.h"
#include "source_base/mathzone_add1.h"

cal_r_overlap_R::cal_r_overlap_R()
{
}

cal_r_overlap_R::~cal_r_overlap_R()
{
}

void cal_r_overlap_R::initialize_orb_table(const UnitCell& ucell,
                                           const LCAO_Orbitals& orb)
{
    int Lmax_used = 0;
    int Lmax = 0;
    int exx_lmax = 0;
#ifdef __EXX
    exx_lmax = GlobalC::exx_info.info_ri.abfs_Lmax;
#endif

    const int ntype = orb.get_ntype();
    int lmax_orb = -1, lmax_beta = -1;
    for (int it = 0; it < ntype; it++)
    {
        lmax_orb = std::max(lmax_orb, orb.Phi[it].getLmax());
        lmax_beta = std::max(lmax_beta, ucell.infoNL.Beta[it].getLmax());
    }
    const double dr = orb.get_dR();
    const double dk = orb.get_dk();
    const int kmesh = orb.get_kmesh() * 4 + 1;
    int Rmesh = static_cast<int>(orb.get_Rmax() / dr) + 4;
    Rmesh += 1 - Rmesh % 2;

    Center2_Orb::init_Table_Spherical_Bessel(2,
                                             3,
                                             Lmax_used,
                                             Lmax,
                                             exx_lmax,
                                             lmax_orb,
                                             lmax_beta,
                                             dr,
                                             dk,
                                             kmesh,
                                             Rmesh,
                                             psb_);
    ModuleBase::Ylm::set_coefficients();
    MGT.init_Gaunt_CH(Lmax);
    MGT.init_Gaunt(Lmax);
}

void cal_r_overlap_R::construct_orbs_and_orb_r(const UnitCell& ucell,
                                               const LCAO_Orbitals& orb)
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
                                               static_cast<int>(orb_origin.getNk() * kmesh_times) | 1,
                                               orb_origin.getDk(),
                                               orb_origin.getDruniform(),
                                               false,
                                               true,
                                               PARAM.inp.cal_force);
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
                           PARAM.inp.cal_force);

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
                            center2_orb21_r[TA][TB][LA][NA][LB].insert(std::make_pair(
                                NB,
                                Center2_Orb::Orb21(orbs[TA][LA][NA], orb_r, orbs[TB][LB][NB], psb_, MGT)));
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

    iw2it.resize(PARAM.globalv.nlocal);
    iw2ia.resize(PARAM.globalv.nlocal);
    iw2iL.resize(PARAM.globalv.nlocal);
    iw2iN.resize(PARAM.globalv.nlocal);
    iw2im.resize(PARAM.globalv.nlocal);

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

void cal_r_overlap_R::construct_orbs_and_nonlocal_and_orb_r(const UnitCell& ucell,const LCAO_Orbitals& orb)
{
    const InfoNonlocal& infoNL_ = ucell.infoNL;

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
                                               static_cast<int>(orb_origin.getNk() * kmesh_times) | 1,
                                               orb_origin.getDk(),
                                               orb_origin.getDruniform(),
                                               false,
                                               true,
                                               PARAM.inp.cal_force);
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
                           PARAM.inp.cal_force);

    orbs_nonlocal.resize(orb.get_ntype());
    for (int T = 0; T < orb.get_ntype(); ++T)
    {
        const int nproj = infoNL_.nproj[T];
        orbs_nonlocal[T].resize(nproj);
        for (int ip = 0; ip < nproj; ip++)
        {
            int nr = infoNL_.Beta[T].Proj[ip].getNr();
            double dr_uniform = 0.01;
	        int nr_uniform = static_cast<int>((infoNL_.Beta[T].Proj[ip].getRadial(nr-1) - infoNL_.Beta[T].Proj[ip].getRadial(0))/dr_uniform) + 1;
            double* rad = new double[nr_uniform];
            double* rab = new double[nr_uniform];
            for (int ir = 0; ir < nr_uniform; ir++)
            {
                rad[ir] = ir*dr_uniform;
                rab[ir] = dr_uniform;
            }
            double* y2 = new double[nr];
            double* Beta_r_uniform = new double[nr_uniform];
            double* dbeta_uniform = new double[nr_uniform];
            ModuleBase::Mathzone_Add1::SplineD2(infoNL_.Beta[T].Proj[ip].getRadial(), infoNL_.Beta[T].Proj[ip].getBeta_r(), nr, 0.0, 0.0, y2);
            ModuleBase::Mathzone_Add1::Cubic_Spline_Interpolation(
                infoNL_.Beta[T].Proj[ip].getRadial(), 
                infoNL_.Beta[T].Proj[ip].getBeta_r(), 
                y2, 
                nr, 
                rad, 
                nr_uniform, 
                Beta_r_uniform, 
                dbeta_uniform
            );

            // linear extrapolation at the zero point
            if (infoNL_.Beta[T].Proj[ip].getRadial(0) > 1e-10)
            {
                double slope = (infoNL_.Beta[T].Proj[ip].getBeta_r(1) - infoNL_.Beta[T].Proj[ip].getBeta_r(0)) / (infoNL_.Beta[T].Proj[ip].getRadial(1) - infoNL_.Beta[T].Proj[ip].getRadial(0));
                Beta_r_uniform[0] = infoNL_.Beta[T].Proj[ip].getBeta_r(0) - slope * infoNL_.Beta[T].Proj[ip].getRadial(0);
            }

            // Here, the operation beta_r / r is performed. To avoid divergence at r=0, beta_r(0) is set to beta_r(1). 
            // However, this may introduce issues, so caution is needed.
            for (int ir = 1; ir < nr_uniform; ir++)
            {
                Beta_r_uniform[ir] = Beta_r_uniform[ir] / rad[ir];
            }
            Beta_r_uniform[0] = Beta_r_uniform[1];

            orbs_nonlocal[T][ip].set_orbital_info(infoNL_.Beta[T].getLabel(),
                                                  infoNL_.Beta[T].getType(),
                                                  infoNL_.Beta[T].Proj[ip].getL(),
                                                  1,
                                                  nr_uniform,
                                                  rab,
                                                  rad,
                                                  Numerical_Orbital_Lm::Psi_Type::Psi,
                                                  Beta_r_uniform,
                                                  static_cast<int>(infoNL_.Beta[T].Proj[ip].getNk() * kmesh_times) | 1,
                                                  infoNL_.Beta[T].Proj[ip].getDk(),
                                                  infoNL_.Beta[T].Proj[ip].getDruniform(),
                                                  false,
                                                  true,
                                                  PARAM.inp.cal_force);

            delete [] rad;
            delete [] rab;
            delete [] y2;
            delete [] Beta_r_uniform;
            delete [] dbeta_uniform;
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
                    for (int ip = 0; ip < infoNL_.nproj[TB]; ip++)
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
                    for (int ip = 0; ip < infoNL_.nproj[TB]; ip++)
                    {
                        center2_orb21_r_nonlocal[TA][TB][LA][NA].insert(std::make_pair(
                            ip,
                            Center2_Orb::Orb21(orbs[TA][LA][NA], orb_r, orbs_nonlocal[TB][ip], psb_, MGT)));
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

    iw2it.resize(PARAM.globalv.nlocal);
    iw2ia.resize(PARAM.globalv.nlocal);
    iw2iL.resize(PARAM.globalv.nlocal);
    iw2iN.resize(PARAM.globalv.nlocal);
    iw2im.resize(PARAM.globalv.nlocal);

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

void cal_r_overlap_R::init(const UnitCell& ucell,const Parallel_Orbitals& pv, const LCAO_Orbitals& orb)
{
    ModuleBase::TITLE("cal_r_overlap_R", "init");
    ModuleBase::timer::tick("cal_r_overlap_R", "init");
    this->ParaV = &pv;

    initialize_orb_table(ucell,orb);
    construct_orbs_and_orb_r(ucell,orb);

    ModuleBase::timer::tick("cal_r_overlap_R", "init");
    return;
}

void cal_r_overlap_R::init_nonlocal(const UnitCell& ucell,const Parallel_Orbitals& pv, const LCAO_Orbitals& orb)
{
    ModuleBase::TITLE("cal_r_overlap_R", "init_nonlocal");
    ModuleBase::timer::tick("cal_r_overlap_R", "init_nonlocal");
    this->ParaV = &pv;

    initialize_orb_table(ucell,orb);
    construct_orbs_and_nonlocal_and_orb_r(ucell,orb);

    ModuleBase::timer::tick("cal_r_overlap_R", "init_nonlocal");
    return;
}

ModuleBase::Vector3<double> cal_r_overlap_R::get_psi_r_psi(const ModuleBase::Vector3<double>& R1,
                                                           const int& T1,
                                                           const int& L1,
                                                           const int& m1,
                                                           const int& N1,
                                                           const ModuleBase::Vector3<double>& R2,
                                                           const int& T2,
                                                           const int& L2,
                                                           const int& m2,
                                                           const int& N2)
{
    ModuleBase::Vector3<double> origin_point(0.0, 0.0, 0.0);
    double factor = sqrt(ModuleBase::FOUR_PI / 3.0);
    const ModuleBase::Vector3<double>& distance = R2 - R1;

    double overlap_o = center2_orb11[T1][T2][L1][N1][L2].at(N2).cal_overlap(origin_point, distance, m1, m2);

    double overlap_x = -1 * factor
                       * center2_orb21_r[T1][T2][L1][N1][L2].at(N2).cal_overlap(origin_point,
                                                                                distance,
                                                                                m1,
                                                                                1,
                                                                                m2); // m =  1

    double overlap_y = -1 * factor
                       * center2_orb21_r[T1][T2][L1][N1][L2].at(N2).cal_overlap(origin_point,
                                                                                distance,
                                                                                m1,
                                                                                2,
                                                                                m2); // m = -1

    double overlap_z = factor
                       * center2_orb21_r[T1][T2][L1][N1][L2].at(N2).cal_overlap(origin_point,
                                                                                distance,
                                                                                m1,
                                                                                0,
                                                                                m2); // m =  0

    ModuleBase::Vector3<double> temp_prp
        = ModuleBase::Vector3<double>(overlap_x, overlap_y, overlap_z) + R1 * overlap_o;

    return temp_prp;
}

void cal_r_overlap_R::get_psi_r_beta(const UnitCell& ucell,
                                     std::vector<std::vector<double>>& nlm,
                                     const ModuleBase::Vector3<double>& R1,
                                     const int& T1,
                                     const int& L1,
                                     const int& m1,
                                     const int& N1,
                                     const ModuleBase::Vector3<double>& R2,
                                     const int& T2)
{
    ModuleBase::Vector3<double> origin_point(0.0, 0.0, 0.0);
    double factor = sqrt(ModuleBase::FOUR_PI / 3.0);
    const ModuleBase::Vector3<double>& distance = R2 - R1;
    const InfoNonlocal& infoNL_ = ucell.infoNL;
    const int nproj = infoNL_.nproj[T2];
    nlm.resize(4);
    if (nproj == 0)
    {
        for(int i = 0;i < 4;i++)
        {
            nlm[i].resize(1);
        }
        return;
    }

    int natomwfc = 0;
    for (int ip = 0; ip < nproj; ip++)
    {
        const int L2 = infoNL_.Beta[T2].Proj[ip].getL(); // mohan add 2021-05-07
        natomwfc += 2 * L2 + 1;
    }
    for(int i = 0;i < 4;i++)
    {
        nlm[i].resize(natomwfc);
    }
    int index = 0;
    for (int ip = 0; ip < nproj; ip++)
    {
        const int L2 = infoNL_.Beta[T2].Proj[ip].getL();
        for (int m2 = 0; m2 < 2 * L2 + 1; m2++)
        {
            double overlap_o
                = center2_orb11_nonlocal[T1][T2][L1][N1].at(ip).cal_overlap(origin_point, distance, m1, m2);

            double overlap_x = -1 * factor
                               * center2_orb21_r_nonlocal[T1][T2][L1][N1].at(ip).cal_overlap(origin_point,
                                                                                             distance,
                                                                                             m1,
                                                                                             1,
                                                                                             m2); // m =  1

            double overlap_y = -1 * factor
                               * center2_orb21_r_nonlocal[T1][T2][L1][N1].at(ip).cal_overlap(origin_point,
                                                                                             distance,
                                                                                             m1,
                                                                                             2,
                                                                                             m2); // m = -1

            double overlap_z = factor
                               * center2_orb21_r_nonlocal[T1][T2][L1][N1].at(ip).cal_overlap(origin_point,
                                                                                             distance,
                                                                                             m1,
                                                                                             0,
                                                                                             m2); // m =  0

            //nlm[index] = ModuleBase::Vector3<double>(overlap_x, overlap_y, overlap_z) + R1 * overlap_o;

            //nlm[index] = ModuleBase::Vector3<double>(overlap_o, overlap_y, overlap_z);// + R1 * overlap_o;
            nlm[0][index] = overlap_o;
            nlm[1][index] = overlap_x + (R1 * overlap_o).x;
            nlm[2][index] = overlap_y + (R1 * overlap_o).y;
            nlm[3][index] = overlap_z + (R1 * overlap_o).z;
            index++;
        }
    }
}


void cal_r_overlap_R::out_rR(const UnitCell& ucell, const Grid_Driver& gd, const int& istep)
{
    ModuleBase::TITLE("cal_r_overlap_R", "out_rR");
    ModuleBase::timer::tick("cal_r_overlap_R", "out_rR");

    int step = istep;
    // set R coor range
    int R_minX = int(-gd.getGlayerX_minus());
    int R_minY = int(-gd.getGlayerY_minus());
    int R_minZ = int(-gd.getGlayerZ_minus());

    int R_x = gd.getGlayerX() + gd.getGlayerX_minus();
    int R_y = gd.getGlayerY() + gd.getGlayerY_minus();
    int R_z = gd.getGlayerZ() + gd.getGlayerZ_minus();

    std::set<Abfs::Vector3_Order<int>> all_R_coor;
    for (int ix = 0; ix < R_x; ix++)
    {
        for (int iy = 0; iy < R_y; iy++)
        {
            for (int iz = 0; iz < R_z; iz++)
            {
                Abfs::Vector3_Order<int> temp_R(ix + R_minX, iy + R_minY, iz + R_minZ);
                all_R_coor.insert(temp_R);
            }
        }
    }

    // calculate rR matrix
    ModuleBase::Vector3<double> tau1, tau2, dtau;
    ModuleBase::Vector3<double> origin_point(0.0, 0.0, 0.0);
    double factor = sqrt(ModuleBase::FOUR_PI / 3.0);
    int output_R_number = 0;

    std::stringstream tem1;
    tem1 << PARAM.globalv.global_out_dir << "tmp-rr.csr";
    std::ofstream ofs_tem1;
    std::ifstream ifs_tem1;

    if (GlobalV::DRANK == 0)
    {
        if (binary)
        {
            ofs_tem1.open(tem1.str().c_str(), std::ios::binary);
        }
        else
        {
            ofs_tem1.open(tem1.str().c_str());
        }
    }

    for (auto& R_coor: all_R_coor)
    {
        std::map<size_t, std::map<size_t, double>> psi_r_psi_sparse[3];

        int dRx = R_coor.x;
        int dRy = R_coor.y;
        int dRz = R_coor.z;

        ModuleBase::Vector3<double> R_car = ModuleBase::Vector3<double>(dRx, dRy, dRz) * ucell.latvec;

        int ir, ic;
        for (int iw1 = 0; iw1 < PARAM.globalv.nlocal; iw1++)
        {
            ir = this->ParaV->global2local_row(iw1);
            if (ir >= 0)
            {
                for (int iw2 = 0; iw2 < PARAM.globalv.nlocal; iw2++)
                {
                    ic = this->ParaV->global2local_col(iw2);
                    if (ic >= 0)
                    {
                        int orb_index_row = iw1 / PARAM.globalv.npol;
                        int orb_index_col = iw2 / PARAM.globalv.npol;

                        // The off-diagonal term in SOC calculaiton is zero, and the two diagonal terms are the same
                        int new_index = iw1 - PARAM.globalv.npol * orb_index_row
                                        + (iw2 - PARAM.globalv.npol * orb_index_col) * PARAM.globalv.npol;

                        if (new_index == 0 || new_index == 3)
                        {
                            int it1 = iw2it[orb_index_row];
                            int ia1 = iw2ia[orb_index_row];
                            int iN1 = iw2iN[orb_index_row];
                            int iL1 = iw2iL[orb_index_row];
                            int im1 = iw2im[orb_index_row];

                            int it2 = iw2it[orb_index_col];
                            int ia2 = iw2ia[orb_index_col];
                            int iN2 = iw2iN[orb_index_col];
                            int iL2 = iw2iL[orb_index_col];
                            int im2 = iw2im[orb_index_col];

                            ModuleBase::Vector3<double> r_distance
                                = (ucell.atoms[it2].tau[ia2] - ucell.atoms[it1].tau[ia1] + R_car)
                                  * ucell.lat0;

                            double overlap_o = center2_orb11[it1][it2][iL1][iN1][iL2].at(iN2).cal_overlap(origin_point,
                                                                                                          r_distance,
                                                                                                          im1,
                                                                                                          im2);

                            double overlap_x
                                = -1 * factor
                                  * center2_orb21_r[it1][it2][iL1][iN1][iL2].at(iN2).cal_overlap(origin_point,
                                                                                                 r_distance,
                                                                                                 im1,
                                                                                                 1,
                                                                                                 im2); // m =  1

                            double overlap_y
                                = -1 * factor
                                  * center2_orb21_r[it1][it2][iL1][iN1][iL2].at(iN2).cal_overlap(origin_point,
                                                                                                 r_distance,
                                                                                                 im1,
                                                                                                 2,
                                                                                                 im2); // m = -1

                            double overlap_z
                                = factor
                                  * center2_orb21_r[it1][it2][iL1][iN1][iL2].at(iN2).cal_overlap(origin_point,
                                                                                                 r_distance,
                                                                                                 im1,
                                                                                                 0,
                                                                                                 im2); // m =  0

                            ModuleBase::Vector3<double> temp_prp
                                = ModuleBase::Vector3<double>(overlap_x, overlap_y, overlap_z)
                                  + ucell.atoms[it1].tau[ia1] * ucell.lat0 * overlap_o;

                            if (std::abs(temp_prp.x) > sparse_threshold)
                            {
                                psi_r_psi_sparse[0][iw1][iw2] = temp_prp.x;
                            }

                            if (std::abs(temp_prp.y) > sparse_threshold)
                            {
                                psi_r_psi_sparse[1][iw1][iw2] = temp_prp.y;
                            }

                            if (std::abs(temp_prp.z) > sparse_threshold)
                            {
                                psi_r_psi_sparse[2][iw1][iw2] = temp_prp.z;
                            }
                        }
                    }
                }
            }
        }

        int rR_nonzero_num[3] = {0, 0, 0};
        for (int direction = 0; direction < 3; ++direction)
        {
            for (auto& row_loop: psi_r_psi_sparse[direction])
            {
                rR_nonzero_num[direction] += row_loop.second.size();
            }
        }

        Parallel_Reduce::reduce_all(rR_nonzero_num, 3);

        if (rR_nonzero_num[0] || rR_nonzero_num[1] || rR_nonzero_num[2])
        {
            output_R_number++;

            if (binary)
            {
                ofs_tem1.write(reinterpret_cast<char*>(&dRx), sizeof(int));
                ofs_tem1.write(reinterpret_cast<char*>(&dRy), sizeof(int));
                ofs_tem1.write(reinterpret_cast<char*>(&dRz), sizeof(int));
            }
            else
            {
                ofs_tem1 << dRx << " " << dRy << " " << dRz << std::endl;
            }

            for (int direction = 0; direction < 3; ++direction)
            {
                if (GlobalV::DRANK == 0)
                {
                    if (binary)
                    {
                        ofs_tem1.write(reinterpret_cast<char*>(&rR_nonzero_num[direction]), sizeof(int));
                    }
                    else
                    {
                        ofs_tem1 << rR_nonzero_num[direction] << std::endl;
                    }
                }

                if (rR_nonzero_num[direction])
                {
                    ModuleIO::output_single_R(ofs_tem1,
                                              psi_r_psi_sparse[direction],
                                              sparse_threshold,
                                              binary,
                                              *(this->ParaV));
                }
                else
                {
                    // do nothing
                }
            }
        }
    }

    if (GlobalV::DRANK == 0)
    {
        std::ofstream out_r;
        std::stringstream ssr;
        if (PARAM.inp.calculation == "md" && !PARAM.inp.out_app_flag)
        {
            ssr << PARAM.globalv.global_matrix_dir
                << "rrg" << step << ".csr";
        }
        else
        {
            ssr << PARAM.globalv.global_out_dir << "rr.csr";
        }

        if (binary) // .dat
        {
            ofs_tem1.close();
            int nlocal = PARAM.globalv.nlocal;
            if (PARAM.inp.calculation == "md" && PARAM.inp.out_app_flag && step)
            {
                out_r.open(ssr.str().c_str(), std::ios::binary | std::ios::app);
            }
            else
            {
                out_r.open(ssr.str().c_str(), std::ios::binary);
            }
            out_r.write(reinterpret_cast<char*>(&step), sizeof(int));
            out_r.write(reinterpret_cast<char*>(&nlocal), sizeof(int));
            out_r.write(reinterpret_cast<char*>(&output_R_number), sizeof(int));

            ifs_tem1.open(tem1.str().c_str(), std::ios::binary);
            out_r << ifs_tem1.rdbuf();
            ifs_tem1.close();
            out_r.close();
        }
        else // .txt
        {
            ofs_tem1.close();
            if (PARAM.inp.calculation == "md" && PARAM.inp.out_app_flag && step)
            {
                out_r.open(ssr.str().c_str(), std::ios::app);
            }
            else
            {
                out_r.open(ssr.str().c_str());
            }
            out_r << "STEP: " << step << std::endl;
            out_r << "Matrix Dimension of r(R): " << PARAM.globalv.nlocal << std::endl;
            out_r << "Matrix number of r(R): " << output_R_number << std::endl;

            ifs_tem1.open(tem1.str().c_str());
            out_r << ifs_tem1.rdbuf();
            ifs_tem1.close();
            out_r.close();
        }

        std::remove(tem1.str().c_str());
    }

    ModuleBase::timer::tick("cal_r_overlap_R", "out_rR");
    return;
}

void cal_r_overlap_R::out_rR_other(const UnitCell& ucell, const int& istep, const std::set<Abfs::Vector3_Order<int>>& output_R_coor)
{
    ModuleBase::TITLE("cal_r_overlap_R", "out_rR_other");
    ModuleBase::timer::tick("cal_r_overlap_R", "out_rR_other");

    // calculate rR matrix
    ModuleBase::Vector3<double> tau1, tau2, dtau;
    ModuleBase::Vector3<double> origin_point(0.0, 0.0, 0.0);
    double factor = sqrt(ModuleBase::FOUR_PI / 3.0);
    int output_R_number = output_R_coor.size();
    int step = istep;

    std::ofstream out_r;
    std::stringstream ssr;
    if (PARAM.inp.calculation == "md" && !PARAM.inp.out_app_flag)
    {
        ssr << PARAM.globalv.global_matrix_dir
            << "rrg" << step << ".csr";
    }
    else
    {
        ssr << PARAM.globalv.global_out_dir << "rr.csr";
    }

    if (GlobalV::DRANK == 0)
    {
        if (binary)
        {
            int nlocal = PARAM.globalv.nlocal;
            if (PARAM.inp.calculation == "md" && PARAM.inp.out_app_flag && step)
            {
                out_r.open(ssr.str().c_str(), std::ios::binary | std::ios::app);
            }
            else
            {
                out_r.open(ssr.str().c_str(), std::ios::binary);
            }
            out_r.write(reinterpret_cast<char*>(&step), sizeof(int));
            out_r.write(reinterpret_cast<char*>(&nlocal), sizeof(int));
            out_r.write(reinterpret_cast<char*>(&output_R_number), sizeof(int));
        }
        else
        {
            if (PARAM.inp.calculation == "md" && PARAM.inp.out_app_flag && step)
            {
                out_r.open(ssr.str().c_str(), std::ios::app);
            }
            else
            {
                out_r.open(ssr.str().c_str());
            }
            out_r << "STEP: " << step << std::endl;
            out_r << "Matrix Dimension of r(R): " << PARAM.globalv.nlocal << std::endl;
            out_r << "Matrix number of r(R): " << output_R_number << std::endl;
        }
    }

    for (auto& R_coor: output_R_coor)
    {
        std::map<size_t, std::map<size_t, double>> psi_r_psi_sparse[3];

        int dRx = R_coor.x;
        int dRy = R_coor.y;
        int dRz = R_coor.z;

        ModuleBase::Vector3<double> R_car = ModuleBase::Vector3<double>(dRx, dRy, dRz) * ucell.latvec;

        int ir = 0;
        int ic = 0;
        for (int iw1 = 0; iw1 < PARAM.globalv.nlocal; iw1++)
        {
            ir = this->ParaV->global2local_row(iw1);
            if (ir >= 0)
            {
                for (int iw2 = 0; iw2 < PARAM.globalv.nlocal; iw2++)
                {
                    ic = this->ParaV->global2local_col(iw2);
                    if (ic >= 0)
                    {
                        int orb_index_row = iw1 / PARAM.globalv.npol;
                        int orb_index_col = iw2 / PARAM.globalv.npol;

                        // The off-diagonal term in SOC calculaiton is zero, and the two diagonal terms are the same
                        int new_index = iw1 - PARAM.globalv.npol * orb_index_row
                                        + (iw2 - PARAM.globalv.npol * orb_index_col) * PARAM.globalv.npol;

                        if (new_index == 0 || new_index == 3)
                        {
                            int it1 = iw2it[orb_index_row];
                            int ia1 = iw2ia[orb_index_row];
                            int iN1 = iw2iN[orb_index_row];
                            int iL1 = iw2iL[orb_index_row];
                            int im1 = iw2im[orb_index_row];

                            int it2 = iw2it[orb_index_col];
                            int ia2 = iw2ia[orb_index_col];
                            int iN2 = iw2iN[orb_index_col];
                            int iL2 = iw2iL[orb_index_col];
                            int im2 = iw2im[orb_index_col];

                            ModuleBase::Vector3<double> r_distance
                                = (ucell.atoms[it2].tau[ia2] - ucell.atoms[it1].tau[ia1] + R_car)
                                  * ucell.lat0;

                            double overlap_o = center2_orb11[it1][it2][iL1][iN1][iL2].at(iN2).cal_overlap(origin_point,
                                                                                                          r_distance,
                                                                                                          im1,
                                                                                                          im2);

                            double overlap_x
                                = -1 * factor
                                  * center2_orb21_r[it1][it2][iL1][iN1][iL2].at(iN2).cal_overlap(origin_point,
                                                                                                 r_distance,
                                                                                                 im1,
                                                                                                 1,
                                                                                                 im2); // m =  1

                            double overlap_y
                                = -1 * factor
                                  * center2_orb21_r[it1][it2][iL1][iN1][iL2].at(iN2).cal_overlap(origin_point,
                                                                                                 r_distance,
                                                                                                 im1,
                                                                                                 2,
                                                                                                 im2); // m = -1

                            double overlap_z
                                = factor
                                  * center2_orb21_r[it1][it2][iL1][iN1][iL2].at(iN2).cal_overlap(origin_point,
                                                                                                 r_distance,
                                                                                                 im1,
                                                                                                 0,
                                                                                                 im2); // m =  0

                            ModuleBase::Vector3<double> temp_prp
                                = ModuleBase::Vector3<double>(overlap_x, overlap_y, overlap_z)
                                  + ucell.atoms[it1].tau[ia1] * ucell.lat0 * overlap_o;

                            if (std::abs(temp_prp.x) > sparse_threshold)
                            {
                                psi_r_psi_sparse[0][iw1][iw2] = temp_prp.x;
                            }

                            if (std::abs(temp_prp.y) > sparse_threshold)
                            {
                                psi_r_psi_sparse[1][iw1][iw2] = temp_prp.y;
                            }

                            if (std::abs(temp_prp.z) > sparse_threshold)
                            {
                                psi_r_psi_sparse[2][iw1][iw2] = temp_prp.z;
                            }
                        }
                    }
                }
            }
        }

        int rR_nonzero_num[3] = {0, 0, 0};
        for (int direction = 0; direction < 3; ++direction)
        {
            for (auto& row_loop: psi_r_psi_sparse[direction])
            {
                rR_nonzero_num[direction] += row_loop.second.size();
            }
        }

        Parallel_Reduce::reduce_all(rR_nonzero_num, 3);

        if (binary) // .dat
        {
            out_r.write(reinterpret_cast<char*>(&dRx), sizeof(int));
            out_r.write(reinterpret_cast<char*>(&dRy), sizeof(int));
            out_r.write(reinterpret_cast<char*>(&dRz), sizeof(int));
        }
        else // .txt
        {
            out_r << dRx << " " << dRy << " " << dRz << std::endl;
        }

        for (int direction = 0; direction < 3; ++direction)
        {
            if (GlobalV::DRANK == 0)
            {
                if (binary)
                {
                    out_r.write(reinterpret_cast<char*>(&rR_nonzero_num[direction]), sizeof(int));
                }
                else
                {
                    out_r << rR_nonzero_num[direction] << std::endl;
                }
            }

            if (rR_nonzero_num[direction])
            {
                ModuleIO::output_single_R(out_r, psi_r_psi_sparse[direction], sparse_threshold, binary, *(this->ParaV));
            }
            else
            {
                // do nothing
            }
        }
    }

    if (GlobalV::DRANK == 0)
    {
        out_r.close();
    }

    ModuleBase::timer::tick("cal_r_overlap_R", "out_rR_other");
    return;
}
