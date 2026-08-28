#include "dftu_lcao.h"
#include "dftu_occup.h"

#include "source_base/matrix.h"  // occ_mat uses ModuleBase::matrix::operator()
#include "source_base/tool_quit.h"
#include "source_base/tool_title.h"
#include "source_base/timer.h"

#include <complex>
#include <vector>

Plus_U::Plus_U()
{}

Plus_U::~Plus_U()
{}

void Plus_U::init(UnitCell& cell,
                const Parallel_Orbitals* pv,
                const int npol,
                const int nspin,
                const std::vector<int>& orbital_corr,
                const bool yukawa_potential,
                const double yukawa_lambda,
                const std::string& global_readin_dir,
                const std::string& global_out_dir,
                const std::string& init_chg,
                const int nlocal,
                const bool gamma_only_local,
                const std::string& ks_solver,
                const bool cal_force,
                const bool cal_stress,
                const std::string& device,
                const int kpar,
                const std::vector<double>& hubbard_u,
                const double uramping,
                const int occ_mat_ctrl,
                const int mixing_dftu
#ifdef __LCAO
                , const LCAO_Orbitals* orb
#endif
                )
{
    ModuleBase::TITLE("Plus_U", "init");

    this->paraV = pv;

    this->yukawa_lambda = yukawa_lambda;
    this->npol = npol;
    this->nlocal = nlocal;
    this->gamma_only_local = gamma_only_local;
    this->ks_solver = ks_solver;
    this->cal_force = cal_force;
    this->cal_stress = cal_stress;

#ifdef __LCAO
    ptr_orb_ = orb;
    if(ptr_orb_ != nullptr)
    {
        orb_cutoff_ = orb->cutoffs();
    }
    ucell = &cell;
#endif

    if (pv != nullptr)
    {
        const int global_rows = pv->get_global_row_size();
        const int global_cols = pv->get_global_col_size();
        if (global_rows != global_cols)
        {
            ModuleBase::WARNING_QUIT("Plus_U::init", "Global row and column dimensions do not match");
        }
        if (nlocal != global_rows)
        {
            ModuleBase::WARNING_QUIT("Plus_U::init", "nlocal does not match global matrix dimension");
        }
    }

    this->init_base(cell,
                    npol,
                    nspin,
                    orbital_corr,
                    yukawa_potential,
                    global_readin_dir,
                    global_out_dir,
                    init_chg,
                    device,
                    kpar,
                    hubbard_u,
                    uramping,
                    occ_mat_ctrl,
                    mixing_dftu);
    return;
}

#ifdef __LCAO

void Plus_U::cal_energy_correction(const UnitCell& ucell,
                                 const int istep)
{
    ModuleBase::TITLE("Plus_U", "cal_energy_correction");
    ModuleBase::timer::start("Plus_U", "cal_energy_correction");
    if (!is_occ_mat_initialized())
    {
        ModuleBase::timer::end("Plus_U", "cal_energy_correction");
        return;
    }

    // mohan update 20251106
    this->energy_u = 0.0;

    double energy_dc = 0.0;

    for (int T = 0; T < ucell.ntype; T++)
    {
        const int NL = ucell.atoms[T].nwl + 1;
        const int LC = get_orbital_corr(T);
        for (int I = 0; I < ucell.atoms[T].na; I++)
        {
            if (LC == -1)
            {
                continue;
            }

            const int iat = ucell.itia2iat(T, I);
            const int L = get_orbital_corr(T);

            for (int l = 0; l < NL; l++)
            {
                if (l != get_orbital_corr(T))
                {
                    continue;
                }

                const int N = ucell.atoms[T].l_nchi[l];

                const int m_tot = 2 * l + 1;

                // part 1: calculate the DFT+U energy correction
                for (int n = 0; n < N; n++)
                {
                    if (n != 0)
                    {
                        continue;
                    }

                    if (this->nspin == 1 || this->nspin == 2)
                    {
                        for (int spin = 0; spin < 2; spin++)
                        {
                            double nm_trace = 0.0;
                            double nm2_trace = 0.0;

                            for (int m0 = 0; m0 < 2 * l + 1; m0++)
                            {
                                nm_trace += this->occ_mat[iat][l][n][spin](m0, m0);
                                for (int m1 = 0; m1 < 2 * l + 1; m1++)
                                {
                                    nm2_trace += this->occ_mat[iat][l][n][spin](m0, m1)
                                                 * this->occ_mat[iat][l][n][spin](m1, m0);
                                }
                            }
                            if (use_yukawa_)
                            {
                                this->energy_u += 0.5 * (this->U_Yukawa[T][l][n] - this->J_Yukawa[T][l][n])
                                            * (nm_trace - nm2_trace);
                            }
                            else
                            {
                                this->energy_u += 0.5 * this->u_current[T] * (nm_trace - nm2_trace);
                            }
                        }
                    }
                    else if (this->nspin == 4)
                    {
                        double nm_trace = 0.0;
                        double nm2_trace = 0.0;

                        for (int m0 = 0; m0 < 2 * l + 1; m0++)
                        {
                            for (int ipol0 = 0; ipol0 < this->npol; ipol0++)
                            {
                                const int m0_all = m0 + (2 * l + 1) * ipol0;
                                nm_trace += this->occ_mat[iat][l][n][0](m0_all, m0_all);

                                for (int m1 = 0; m1 < 2 * l + 1; m1++)
                                {
                                    for (int ipol1 = 0; ipol1 < this->npol; ipol1++)
                                    {
                                        int m1_all = m1 + (2 * l + 1) * ipol1;

                                        nm2_trace += this->occ_mat[iat][l][n][0](m0_all, m1_all)
                                                     * this->occ_mat[iat][l][n][0](m1_all, m0_all);
                                    }
                                }
                            }
                        }
                        if (use_yukawa_)
                        {
                            this->energy_u += 0.5 * (this->U_Yukawa[T][l][n] - this->J_Yukawa[T][l][n]) 
                              * (nm_trace - nm2_trace);
                        }
                        else
                        {
                            this->energy_u += 0.5 * this->u_current[T] * (nm_trace - nm2_trace);
                        }
                    }

                    for (int m1 = 0; m1 < 2 * l + 1; m1++)
                    {
                        for (int ipol1 = 0; ipol1 < this->npol; ipol1++)
                        {
                            const int m1_all = m1 + ipol1 * (2 * l + 1);
                            for (int m2 = 0; m2 < 2 * l + 1; m2++)
                            {
                                for (int ipol2 = 0; ipol2 < this->npol; ipol2++)
                                {
                                    const int m2_all = m2 + ipol2 * (2 * l + 1);

                                    if (this->nspin == 1 || this->nspin == 2)
                                    {
                                        for (int is = 0; is < 2; is++)
                                        {
                                            double pot_onsite = 0.0;
                                            pot_onsite = get_onebody_eff_pot(T, iat, l, n, is, m1_all, m2_all, false);
                                            energy_dc += pot_onsite * this->occ_mat[iat][l][n][is](m1_all, m2_all);
                                        }
                                    }
                                    else if (this->nspin == 4)
                                    {
                                        double pot_onsite = 0.0;
                                        pot_onsite = get_onebody_eff_pot(T, iat, l, n, 0, m1_all, m2_all, false);
                                        energy_dc += pot_onsite * this->occ_mat[iat][l][n][0](m1_all, m2_all);
                                    }
                                }
                            }
                        }
                    }
                } // end n
            }     // end L
        }         // end I
    }             // end T

    // substract the double counting energy_dc included in band energy eband
    this->energy_u -= energy_dc;

    ModuleBase::timer::end("Plus_U", "cal_energy_correction");
    return;
}

#endif

// uramping_update() and u_converged() are now implemented in
// dftu_base.cpp as Plus_U_Base methods (inherited by Plus_U).

#ifdef __LCAO

void Plus_U::set_dmr(const elecstate::DensityMatrix<std::complex<double>, double>* dmr)
{
    this->dm_in_dftu_cd = dmr;
    return;
}

void Plus_U::set_dmr(const elecstate::DensityMatrix<double, double>* dmr)
{
    this->dm_in_dftu_d = dmr;
    return;
}

const hamilt::HContainer<double>* Plus_U::get_dmr(int ispin) const
{
    if (this->dm_in_dftu_d != nullptr)
    {
        return this->dm_in_dftu_d->get_DMR_pointer(ispin + 1);
    }
    else if (this->dm_in_dftu_cd != nullptr)
    {
        return this->dm_in_dftu_cd->get_DMR_pointer(ispin + 1);
    }
    else
    {
        return nullptr;
    }
}

namespace DFTU_LCAO {

//! dftu occupation matrix for gamma only using dm(double)
template <>
void cal_occ_mat(const int iter,
                 const UnitCell& ucell,
                 const std::vector<std::vector<double>>& dm,
                 const K_Vectors& kv,
                 const double& mixing_beta,
                 hamilt::Hamilt<double>* p_ham,
                 Plus_U& dftu)
{
    dftu.cal_occ_mat_gamma(iter, ucell, dm, mixing_beta, p_ham);
}

//! dftu occupation matrix for multiple k-points using dm(complex)
template <>
void cal_occ_mat(const int iter,
                 const UnitCell& ucell,
                 const std::vector<std::vector<std::complex<double>>>& dm,
                 const K_Vectors& kv,
                 const double& mixing_beta,
                 hamilt::Hamilt<std::complex<double>>* p_ham,
                 Plus_U& dftu)
{
    dftu.cal_occ_mat_k(iter, ucell, dm, kv, mixing_beta, p_ham);
}

} // namespace DFTU_LCAO

#endif
