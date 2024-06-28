#include "ORB_gen_tables.h"

#include "ORB_read.h"
#include "module_base/constants.h"
#include "module_base/math_integral.h"
#include "module_base/math_polyint.h"
#include "module_base/timer.h"
#include "module_base/ylm.h"

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

    //////////////////////////////
    /// (1) MOT: make overlap table.
    //////////////////////////////
    MOT.allocate(orb.get_ntype(), orb.get_lmax(), orb.get_kmesh(), orb.get_Rmax(), orb.get_dR(), orb.get_dk());

    // OV: overlap
    MOT.init_OV_Tpair(orb);
    MOT.init_OV_Opair(orb);

    //////////////////////////////
    /// (2) init Ylm Coef
    //////////////////////////////
    // liaochen add 2010/4/29
    ModuleBase::Ylm::set_coefficients();

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

    MOT.init_Table_Spherical_Bessel(orb_num, mode, Lmax_used, Lmax, Lmax_exx, orb, beta_, MOT.pSB);

    // calculate S(R) for interpolation
    MOT.init_Table(orb);

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

void ORB_gen_tables::snap_psipsi(const LCAO_Orbitals& orb,
                                 double olm[],
                                 const int& job,    // 0, 1
                                 const char& dtype, // derivative type: S or T
                                 const ModuleBase::Vector3<double>& R1,
                                 const int& T1,
                                 const int& L1,
                                 const int& m1,
                                 const int& N1,
                                 const ModuleBase::Vector3<double>& R2,
                                 const int& T2,
                                 const int& L2,
                                 const int& m2,
                                 const int& N2,
                                 bool cal_syns,
                                 double dmax) const
{
    // ModuleBase::TITLE("ORB_gen_tables","snap_psipsi");
    // ModuleBase::timer::tick ("ORB_gen_tables", "snap_psipsi");
    if (job != 0 && job != 1)
    {
        ModuleBase::WARNING_QUIT("ORB_gen_tables::snap_psipsi", "job must be equal to 0 or 1!");
    }

    Numerical_Orbital_AtomRelation noar;
    noar.set_position(R1, R2);
    assert(this->lat0 > 0.0);

    /// (1) get distance between R1 and R2 (a.u.)
    /// judge if there exist overlap
    double distance = noar.get_distance() * this->lat0;

    const double Rcut1 = orb.Phi[T1].getRcut();
    const double Rcut2 = orb.Phi[T2].getRcut(); // caoyu modified 2021-05-08

    if (job == 0)
    {
        ModuleBase::GlobalFunc::ZEROS(olm, 1);
    }
    else if (job == 1)
    {
        ModuleBase::GlobalFunc::ZEROS(olm, 3);
    }

    if (distance > (Rcut1 + Rcut2))
        return;

    /// if distance == 0,
    /// \f[ \int psi(r) psi(r-R)\f] dr independent of R if R == 0.
    /// distance += tiny1 avoid overflow during calculation.
    const double tiny1 = 1e-12;
    double tiny2 = 1e-10;
    if (distance < tiny1)
        distance += tiny1;

    if (cal_syns)
        tiny2 = dmax;

    /// (2) if there exist overlap, calculate the mesh number
    /// between two atoms
    const int rmesh = ORB_table_phi::get_rmesh(Rcut1, Rcut2, MOT.dr); // caoyu modified 2021-05-08

    /// (3) Find three dimension of 'Table_S' or 'Table_T'.
    /// -dim1 : type pairs,
    /// -dim2 : radial orbital pairs,
    /// -dim3 : find lmax between T1 and T2, and get lmax*2+1
    int dim2 = 0;

    int dim1 = this->MOT.OV_Tpair(T1, T2);
    int dim3 = this->MOT.OV_L2plus1(T1, T2); // 2*lmax+1
    if (T1 <= T2)
    {
        dim2 = this->MOT.OV_Opair(dim1, L1, L2, N1, N2);
    }
    else
    {
        dim2 = this->MOT.OV_Opair(dim1, L2, L1, N2, N1);
    }

    // Gaunt Index
    const int gindex1 = L1 * L1 + m1;
    const int gindex2 = L2 * L2 + m2;

    // Peize Lin change rly, grly 2016-08-26
    std::vector<double> rly;
    std::vector<std::vector<double>> grly;

    //	double *ylm = new double[nlm];
    //	dR = R1 - R2;
    double arr_dR[3];
    arr_dR[0] = noar.getX() * this->lat0;
    arr_dR[1] = noar.getY() * this->lat0;
    arr_dR[2] = noar.getZ() * this->lat0;

    // double xdr = arr_dR[0] / distance;
    // double ydr = arr_dR[1] / distance;
    // double zdr = arr_dR[2] / distance;

    //=======================
    // *r**l*Ylm_real
    // include its derivations
    //=======================
    if (job == 0)
    {
        //		Ylm::rlylm(dim3, arr_dR[0], arr_dR[1], arr_dR[2], rly);
        //		Ylm::sph_harm (dim3-1, xdr, ydr, zdr, rly);
        ModuleBase::Ylm::rl_sph_harm(dim3 - 1, arr_dR[0], arr_dR[1], arr_dR[2], rly);
    }
    else
    {
        //		Ylm::rlylm(dim3, arr_dR[0], arr_dR[1], arr_dR[2], rly, grly);
        ModuleBase::Ylm::grad_rl_sph_harm(dim3 - 1, arr_dR[0], arr_dR[1], arr_dR[2], rly, grly);
    }

    switch (dtype)
    {
    case 'S':
        for (int L = 0; L < dim3; L++) // maxL = dim3-1
        {
            //===========================================================
            // triangle rule for L and sum of L, L1, L2 should be even
            //===========================================================
            int AL = L1 + L2;
            int SL = std::abs(L1 - L2);

            if ((L > AL) || (L < SL) || ((L - SL) % 2 == 1))
                continue;

            double Interp_Slm = 0.0;
            double Interp_dSlm = 0.0;
            double tmpOlm0 = 0.0;
            double tmpOlm1 = 0.0;

            // prefactor
            double i_exp = pow(-1.0, (L1 - L2 - L) / 2);
            double rl = pow(distance, L);

            if (distance > tiny2)
            {
                Interp_Slm = i_exp
                             * ModuleBase::PolyInt::Polynomial_Interpolation(MOT.Table_SR[0][dim1][dim2][L],
                                                                             rmesh,
                                                                             MOT.dr,
                                                                             distance);
                Interp_Slm /= rl;
            }
            else // distance = 0.0;
            {
                Interp_Slm = i_exp * MOT.Table_SR[0][dim1][dim2][L][0];
            }

            if (job == 1) // calculate the derivative.
            {
                if (distance > tiny2)
                {
                    Interp_dSlm = i_exp
                                  * ModuleBase::PolyInt::Polynomial_Interpolation(MOT.Table_SR[1][dim1][dim2][L],
                                                                                  rmesh,
                                                                                  MOT.dr,
                                                                                  distance);
                    Interp_dSlm = Interp_dSlm / pow(distance, L) - Interp_Slm * L / distance;
                }
                else
                {
                    Interp_dSlm = 0.0;
                }
            }

            for (int m = 0; m < 2 * L + 1; m++)
            {
                int gindex = L * L + m;
                //			double tmpGaunt1 = MGT.Get_Gaunt_SH(L1, m1, L2, m2, L, m);
                double tmpGaunt = MGT.Gaunt_Coefficients(gindex1, gindex2, gindex);

                tmpOlm0 = Interp_Slm * tmpGaunt;

                if (job == 1)
                {
                    tmpOlm1 = Interp_dSlm * tmpGaunt;
                }

                switch (job)
                {
                case 0: // calculate overlap.
                {
                    olm[0] += tmpOlm0 * rly[MGT.get_lm_index(L, m)];

                    /*
                        if( abs ( tmpOlm0 * rly[ MGT.get_lm_index(L, m) ] ) > 1.0e-3 )
                        {
                        std::cout << " L=" << L << " m=" << m << " tmpOlm0=" << tmpOlm0
                        << " rly=" << rly[ MGT.get_lm_index(L, m) ]
                        << " r=" << olm[0]
                        << std::endl;
                        }
                        */
                    break;
                }
                case 1: // calculate gradient.
                {
                    for (int ir = 0; ir < 3; ir++)
                    {
                        olm[ir] += tmpOlm0 * grly[MGT.get_lm_index(L, m)][ir]
                                   + tmpOlm1 * rly[MGT.get_lm_index(L, m)] * arr_dR[ir] / distance;
                    }
                    break;
                }
                default:
                    break;
                }
            } // m
        }
        break;

    case 'T':
        for (int L = 0; L < dim3; L++)
        {
            int AL = L1 + L2;
            int SL = std::abs(L1 - L2);

            if ((L > AL) || (L < SL) || ((L - SL) % 2 == 1))
                continue;

            double Interp_Tlm, Interp_dTlm, tmpKem0, tmpKem1;
            Interp_Tlm = Interp_dTlm = tmpKem0 = tmpKem1 = 0.0;

            // pre-fac
            double i_exp = pow(-1.0, (L1 - L2 - L) / 2);

            double rl = pow(distance, L);
            if (distance > tiny2)
            {
                Interp_Tlm = i_exp
                             * ModuleBase::PolyInt::Polynomial_Interpolation(MOT.Table_TR[0][dim1][dim2][L],
                                                                             rmesh,
                                                                             MOT.dr,
                                                                             distance);
                Interp_Tlm /= rl;
            }
            else
                Interp_Tlm = i_exp * MOT.Table_TR[0][dim1][dim2][L][0];

            if (job == 1)
            {
                if (distance > tiny2)
                {
                    Interp_dTlm = i_exp
                                  * ModuleBase::PolyInt::Polynomial_Interpolation(MOT.Table_TR[1][dim1][dim2][L],
                                                                                  rmesh,
                                                                                  MOT.dr,
                                                                                  distance);
                    Interp_dTlm = Interp_dTlm / rl - Interp_Tlm * L / distance;
                }
                else
                    Interp_dTlm = 0.0;
            }

            for (int m = 0; m < 2 * L + 1; m++)
            {
                int gindex = L * L + m;
                //	double tmpGaunt = MGT.Get_Gaunt_SH(L1, m1, L2, m2, L, m);
                double tmpGaunt = MGT.Gaunt_Coefficients(gindex1, gindex2, gindex);

                tmpKem0 = Interp_Tlm * tmpGaunt;
                if (job == 1)
                {
                    tmpKem1 = Interp_dTlm * tmpGaunt;
                }

                switch (job)
                {
                case 0: {
                    olm[0] += tmpKem0 * rly[MGT.get_lm_index(L, m)];
                    break;
                }
                case 1: {
                    for (int ir = 0; ir < 3; ir++)
                    {
                        olm[ir] += tmpKem0 * grly[MGT.get_lm_index(L, m)][ir]
                                   + tmpKem1 * rly[MGT.get_lm_index(L, m)] * arr_dR[ir] / distance;
                    }
                    break;
                }
                default:
                    break;
                }
            } // end T: m
        }     // end T: :
        break;
    }
    //	ModuleBase::timer::tick ("ORB_gen_tables", "snap_psipsi");
    return;
}

double ORB_gen_tables::get_distance(const ModuleBase::Vector3<double>& R1, const ModuleBase::Vector3<double>& R2) const
{
    assert(this->lat0 > 0.0);
    ModuleBase::Vector3<double> dR = R1 - R2;
    return dR.norm() * this->lat0;
}
