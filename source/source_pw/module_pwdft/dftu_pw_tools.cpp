#include "source_pw/module_pwdft/dftu_pw_tools.h"

#include "source_pw/module_pwdft/dftu_base.h"
#include "source_cell/unitcell.h"
#include "source_base/parallel_reduce.h"
#include "source_base/global_variable.h"

namespace pw {

void pauli_to_spin_basis(std::complex<double>* pot_onsite, int m_size)
{
    const int size = m_size * m_size;
    for (int m1 = 0; m1 < m_size; m1++)
    {
        for (int m2 = 0; m2 < m_size; m2++)
        {
            int index[4];
            index[0] = m1 * m_size + m2;
            index[1] = m1 * m_size + m2 + size;
            index[2] = m1 * m_size + m2 + size * 2;
            index[3] = m1 * m_size + m2 + size * 3;
            std::complex<double> pot_onsite_tmp[4];
            for (int i = 0; i < 4; i++)
            {
                pot_onsite_tmp[i] = pot_onsite[index[i]];
            }
            pot_onsite[index[0]] = 0.5 * (pot_onsite_tmp[0] + pot_onsite_tmp[3]);
            pot_onsite[index[3]] = 0.5 * (pot_onsite_tmp[0] - pot_onsite_tmp[3]);
            pot_onsite[index[1]] = 0.5 * (pot_onsite_tmp[1] + std::complex<double>(0.0, 1.0) * pot_onsite_tmp[2]);
            pot_onsite[index[2]] = 0.5 * (pot_onsite_tmp[1] - std::complex<double>(0.0, 1.0) * pot_onsite_tmp[2]);
        }
    }
}

double compute_pot_onsite_spinor(
    std::complex<double>* pot_onsite,
    const double* occ,
    double u_value,
    double diag_coeff,
    double weight_eu,
    int m_size)
{
    double energy_u = 0.0;
    const int m_size2 = m_size * m_size;
    for (int is = 0; is < 4; ++is)
    {
        int start = is * m_size2;
        double diag = (is == 0) ? diag_coeff : 0.0;
        for (int m1 = 0; m1 < m_size; m1++)
        {
            for (int m2 = 0; m2 < m_size; m2++)
            {
                pot_onsite[start + m1 * m_size + m2] = u_value *
                    (diag * (m1 == m2) - occ[start + m2 * m_size + m1]);
                energy_u += u_value * weight_eu
                    * occ[start + m2 * m_size + m1]
                    * occ[start + m1 * m_size + m2];
            }
        }
    }
    pauli_to_spin_basis(pot_onsite, m_size);
    return energy_u;
}

double compute_pot_onsite_scalar(
    std::complex<double>* pot_onsite,
    const double* occ,
    double u_value,
    double diag_coeff,
    double weight_eu,
    int m_size)
{
    double energy_u = 0.0;
    for (int m1 = 0; m1 < m_size; m1++)
    {
        for (int m2 = 0; m2 < m_size; m2++)
        {
            pot_onsite[m1 * m_size + m2] = u_value *
                (diag_coeff * (m1 == m2) - occ[m2 * m_size + m1]);
            energy_u += u_value * weight_eu
                * occ[m2 * m_size + m1]
                * occ[m1 * m_size + m2];
        }
    }
    return energy_u;
}

void accumulate_occ_spinor(
    double* occ_mat_out,
    const std::complex<double>* becp,
    int nbands,
    int npol,
    int nkb,
    int begin_ih,
    int m_begin,
    int tlp1,
    const ModuleBase::matrix& wg,
    int ik)
{
    const int tlp1_2 = tlp1 * tlp1;
    for (int ib = 0; ib < nbands; ib++)
    {
        const double weight = wg(ik, ib);
        int ind_m1m2 = 0;
        for (int m1 = 0; m1 < tlp1; m1++)
        {
            const int index_m1 = ib * npol * nkb + begin_ih + m_begin + m1;
            for (int m2 = 0; m2 < tlp1; m2++)
            {
                const int index_m2 = ib * npol * nkb + begin_ih + m_begin + m2;
                std::complex<double> occ[4];
                occ[0] = weight * std::conj(becp[index_m1]) * becp[index_m2];
                occ[1] = weight * std::conj(becp[index_m1]) * becp[index_m2 + nkb];
                occ[2] = weight * std::conj(becp[index_m1 + nkb]) * becp[index_m2];
                occ[3] = weight * std::conj(becp[index_m1 + nkb]) * becp[index_m2 + nkb];
                occ_mat_out[ind_m1m2] += (occ[0] + occ[3]).real();
                occ_mat_out[ind_m1m2 + tlp1_2] += (occ[1] + occ[2]).real();
                occ_mat_out[ind_m1m2 + 2 * tlp1_2] += (occ[1] - occ[2]).imag();
                occ_mat_out[ind_m1m2 + 3 * tlp1_2] += (occ[0] - occ[3]).real();
                ind_m1m2++;
            }
        }
    }
}

void accumulate_occ_scalar(
    double* occ_mat_out,
    const std::complex<double>* becp,
    int nbands,
    int nkb,
    int begin_ih,
    int m_begin,
    int tlp1,
    const ModuleBase::matrix& wg,
    int ik)
{
    for (int ib = 0; ib < nbands; ib++)
    {
        const double weight = wg(ik, ib);
        int ind_m1m2 = 0;
        for (int m1 = 0; m1 < tlp1; m1++)
        {
            const int index_m1 = ib * nkb + begin_ih + m_begin + m1;
            for (int m2 = 0; m2 < tlp1; m2++)
            {
                const int index_m2 = ib * nkb + begin_ih + m_begin + m2;
                occ_mat_out[ind_m1m2] += weight * (std::conj(becp[index_m1]) * becp[index_m2]).real();
                ind_m1m2++;
            }
        }
    }
}

void reduce_occ_mat(const UnitCell& cell,
                    const int nspin,
                    const int kpar,
                    const std::vector<int>& l_channel,
                    OccupationMatrix& occmat)
{
    for(int iat = 0; iat < cell.nat; iat++)
    {
        const int it = cell.iat2it[iat];
        const int target_l = l_channel[it];
        if(target_l == -1)
        {
            continue;
        }
        const int size = (2 * target_l + 1) * (2 * target_l + 1);

        if(nspin != 4)
        {
            Parallel_Reduce::reduce_double_allpool(kpar,
                    GlobalV::NPROC_IN_POOL,
                    occmat.mat(iat, target_l, 0, 0).c,
                    size);
            if(nspin == 2)
            {
                Parallel_Reduce::reduce_double_allpool(kpar,
                        GlobalV::NPROC_IN_POOL,
                        occmat.mat(iat, target_l, 0, 1).c,
                        size);
            }
        }
        else
        {
            Parallel_Reduce::reduce_double_allpool(kpar,
                    GlobalV::NPROC_IN_POOL,
                    occmat.mat(iat, target_l, 0, 0).c,
                    size * 4);
        }
    }
}

void compute_pot_uterm_and_energy(const UnitCell& cell,
                                  const int nspin,
                                  const std::vector<double>& u_current,
                                  const std::vector<int>& l_channel,
                                  const std::vector<int>& uterm_mat_index,
                                  const OccupationMatrix& occmat,
                                  std::vector<std::complex<double>>& uterm_mat,
                                  double& energy_u)
{
    energy_u = 0.0;
    const double weight_eu = (nspin == 1) ? 1.0 : (nspin == 2) ? 0.5 : 0.25;
    const double diag_coeff = (nspin == 4) ? 1.0 : 0.5;
    // calculate pot_onsite and energy (occ_mat already reduced above)
    for(int iat = 0; iat < cell.nat; iat++)
    {
        const int it = cell.iat2it[iat];
        const int target_l = l_channel[it];
        if(target_l == -1)
        {
            continue;
        }
        const int size = (2 * target_l + 1) * (2 * target_l + 1);

        //update effective potential
        const double u_value = u_current[it];
        std::complex<double>* pot_onsite_iat = &(uterm_mat[uterm_mat_index[iat]]);
        const int m_size = 2 * target_l + 1;

        if(nspin == 4)
        {
            // pot_onsite is stored as 4 contiguous Pauli blocks per atom:
            //   is=0: charge channel (identity), Hubbard U contributes the
            //         diagonal term diag_coeff*delta(m1,m2)
            //   is=1,2,3: spin channels (sigma_x/y/z), no U diagonal term
            // The occupation matrix occ_mat[...][0][0].c packs all 4 blocks
            // contiguously, each of size m_size*m_size.
            energy_u += compute_pot_onsite_spinor(
                pot_onsite_iat,
                occmat.mat(iat, target_l, 0, 0).c,
                u_value, diag_coeff, weight_eu, m_size);
        }
        else // nspin=1 or nspin=2
        {
            // spin-up channel
            energy_u += compute_pot_onsite_scalar(
                pot_onsite_iat,
                occmat.mat(iat, target_l, 0, 0).c,
                u_value, diag_coeff, weight_eu, m_size);
            // spin-down channel for nspin=2
            if(nspin == 2)
            {
                std::complex<double>* pot_onsite_iat1 = &(uterm_mat[uterm_mat.size()/2 + uterm_mat_index[iat]]);
                energy_u += compute_pot_onsite_scalar(
                    pot_onsite_iat1,
                    occmat.mat(iat, target_l, 0, 1).c,
                    u_value, diag_coeff, weight_eu, m_size);
            }
        }
    }
}

} // namespace pw
