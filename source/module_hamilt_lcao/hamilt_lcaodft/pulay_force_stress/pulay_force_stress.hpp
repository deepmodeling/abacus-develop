#pragma once
#include "pulay_force_stress.h"
#include "module_hamilt_lcao/hamilt_lcaodft/stress_tools.h"
namespace PulayForceStress
{
    // inline std::vector<double> single_derivative(
    //     const TwoCenterBundle& c2,
    //     const char& dtype,
    //     const int& T1, const int& L1, const int& N1, const int& m1,
    //     const int& T2, const int& L2, const int& N2, const int& m2,
    //     const ModuleBase::Vector3<double>& dtau_lat0)
    // {
    //     std::vector<double> dhs_xyz(3);
    //     // convert m (0,1,...2l) to M (-l, -l+1, ..., l-1, l)
    //     const int M1 = (m1 % 2 == 0) ? -m1 / 2 : (m1 + 1) / 2;
    //     const int M2 = (m2 % 2 == 0) ? -m2 / 2 : (m2 + 1) / 2;
    //     switch (dtype)
    //     {
    //     case 'S':
    //         two_center_bundle.overlap_orb->calculate(T1, L1, N1, M1, T2, L2, N2, M2, dtau_lat0, nullptr, dhs_xyz.data());
    //         break;
    //     case 'T':
    //         two_center_bundle.kinetic_orb->calculate(T1, L1, N1, M1, T2, L2, N2, M2, dtau_lat0, nullptr, dhs_xyz.data());
    //         break;
    //         // how about nonlocal?
    //     default: // not supposed to happen
    //         ModuleBase::WARNING_QUIT("LCAO_domain::build_ST_new", "dtype must be S or T");
    //     }
    // }

    // inline std::tuple<std::vector<std::vector<double>>, std::vector<std::vector<std::vector<double>>>> get_dHS
    // (
    //     const TwoCenterIntegrator& c2,
    //     const char& dtype,
    //     const bool& isstress
    // )   // replace build_ST_new
    // {
    //     std::vector<std::vector<double>> dHSx;
    //     std::vector<std::vector<std::vector<double>>> dHSxy;
    //     return std::make_tuple(dHSx, dHSxy);
    // }
    // template<typename TK, typename TR>
    // void cal_pulay_fs(
    //     ModuleBase::matrix& f,  ///< [out] force
    //     ModuleBase::matrix& s,  ///< [out] stress
    //     const elecstate::DensityMatrix<TK, TR>& dm,  ///< [in] density matrix
    //     const UnitCell& ucell,  ///< [in] unit cell
    //     const TwoCenterBundle& c2,
    //     const bool& isstress)
    // {
    //     const int nspin = GlobalV::NSPIN;
    //     const int nbasis = GlobalV::NLOCAL;
    //     // any better container?
    //     const auto& dHSs = get_dHSx(c2);
    //     std::vector<std::vector<double>> dHSx = std::get<0>(dHSs);
    //     std::vector<std::vector<std::vector<double>>> dHSxy = std::get<1>(dHSs);
    //     for (int i = 0; i < nbasis; i++)
    //     {
    //         const int iat = ucell.iwt2iat[i];
    //         for (int j = 0; j < nbasis; j++)
    //         {
    //             const int mu = pv.global2local_row(j);
    //             const int nu = pv.global2local_col(i);

    //             if (mu >= 0 && nu >= 0)
    //             {
    //                 const int index = mu * pv.ncol + nu;
    //                 double sum = 0.0;
    //                 for (int is = 0; is < nspin; ++is)
    //                 {
    //                     sum += dm.get_DMK(is + 1, 0, nu, mu);
    //                 }
    //                 for (int i = 0;i < 3;++i) { f(iat, i) += 2 * sum * dHSx[i][index]; }
    //                 if (isstress)
    //                 {
    //                     for (int i = 0;i < 3;++i) { for (int j = i;j < 3;++j) { s(i, j) += 2 * sum * dHSxy[i][j][index]; } }
    //                 }
    //             }
    //         }
    //     }

    //     if (isstress)
    //     {
    //         StressTools::stress_fill(ucell.lat0, ucell.omega, s);
    //     }
    // }
    template<typename TK, typename TR>
    void cal_pulay_fs(
        ModuleBase::matrix& f,  ///< [out] force
        ModuleBase::matrix& s,  ///< [out] stress
        const elecstate::DensityMatrix<TK, TR>& dm,  ///< [in] density matrix
        const UnitCell& ucell,  ///< [in] unit cell
        const elecstate::Potential* pot, ///< [in] potential on grid
        typename TGint<TK>::type& gint,
        const bool& isforce,
        const bool& isstress,
        const bool& set_dmr_gint)
    {
        if (set_dmr_gint) { gint.transfer_DM2DtoGrid(dm.get_DMR_vector()); }    // 2d block to grid

        const int nspin = GlobalV::NSPIN;
        for (int is = 0; is < nspin; ++is)
        {
            const double* vr_eff1 = pot->get_effective_v(is);
            const double* vofk_eff1 = nullptr;

            if (XC_Functional::get_func_type() == 3 || XC_Functional::get_func_type() == 5)
            {
                vofk_eff1 = pot->get_effective_vofk(is);

                Gint_inout inout(
                    is,
                    vr_eff1,
                    vofk_eff1,
                    isforce,
                    isstress,
                    &f,
                    &s,
                    Gint_Tools::job_type::force_meta);

                gint.cal_gint(&inout);
            }
            else
            {
                Gint_inout inout(
                    is,
                    vr_eff1,
                    isforce,
                    isstress,
                    &f,
                    &s,
                    Gint_Tools::job_type::force);

                gint.cal_gint(&inout);
            }
        }

        if (isstress) { StressTools::stress_fill(-1.0, ucell.omega, s); }
    }
}