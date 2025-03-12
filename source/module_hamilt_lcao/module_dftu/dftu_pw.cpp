#include "dftu.h"
#include "module_hamilt_pw/hamilt_pwdft/onsite_projector.h"
#include "module_base/parallel_reduce.h"
#include "module_parameter/parameter.h"
#include "module_base/timer.h"


namespace ModuleDFTU
{
DFTU* DFTU::get_instance()
{
    return &GlobalC::dftu;
}
/// calculate occupation matrix for DFT+U
void DFTU::cal_occ_pw(const int iter, const void* psi_in, const ModuleBase::matrix& wg_in, const UnitCell& cell, const double& mixing_beta)
{
    ModuleBase::timer::tick("DFTU", "cal_occ_pw");
    this->copy_locale(cell);
    this->zero_locale(cell);

    if(PARAM.inp.device == "cpu")
    {
        auto* onsite_p = projectors::OnsiteProjector<double, base_device::DEVICE_CPU>::get_instance();
        const psi::Psi<std::complex<double>>* psi_p = (const psi::Psi<std::complex<double>>*)psi_in;
        // loop over k-points to calculate Mi of \sum_{k,i,l,m}<Psi_{k,i}|alpha_{l,m}><alpha_{l,m}|Psi_{k,i}>
        const int nbands = psi_p->get_nbands();
        for(int ik = 0; ik < psi_p->get_nk(); ik++)
        {
            psi_p->fix_k(ik);
            onsite_p->tabulate_atomic(ik);

            onsite_p->overlap_proj_psi(nbands*psi_p->get_npol(), psi_p->get_pointer());
            const std::complex<double>* becp = onsite_p->get_h_becp();
            // becp(nbands*npol , nkb)
            // mag = wg * \sum_{nh}becp * becp
            int nkb = onsite_p->get_size_becp() / nbands / psi_p->get_npol();
            int begin_ih = 0;
            for(int iat = 0; iat < cell.nat; iat++)
            {
                const int it = cell.iat2it[iat];
                const int nh = onsite_p->get_nh(iat);
                const int target_l = this->orbital_corr[it];
                if(target_l == -1)
                {
                    begin_ih += nh;
                    continue;
                }
                // m = l^2, l^2+1, ..., (l+1)^2-1
                const int m_begin = target_l * target_l;
                const int tlp1 = 2 * target_l + 1;
                const int tlp1_2 = tlp1 * tlp1;
                for(int ib = 0;ib<nbands;ib++)
                {
                    const double weight = wg_in(ik, ib);
                    int ind_m1m2 = 0;
                    for(int m1 = 0; m1 < tlp1; m1++)
                    {
                        const int index_m1 = ib*2*nkb + begin_ih + m_begin + m1;
                        for(int m2 = 0; m2 < tlp1; m2++)
                        {
                            const int index_m2 = ib*2*nkb + begin_ih + m_begin + m2;
                            std::complex<double> occ[4];
                            occ[0] = weight * conj(becp[index_m1]) * becp[index_m2];
                            occ[1] = weight * conj(becp[index_m1]) * becp[index_m2 + nkb];
                            occ[2] = weight * conj(becp[index_m1 + nkb]) * becp[index_m2];
                            occ[3] = weight * conj(becp[index_m1 + nkb]) * becp[index_m2 + nkb];
                            this->locale[iat][target_l][0][0].c[ind_m1m2] += (occ[0] + occ[3]).real();
                            this->locale[iat][target_l][0][0].c[ind_m1m2 + tlp1_2] += (occ[1] + occ[2]).real();
                            this->locale[iat][target_l][0][0].c[ind_m1m2 + 2 * tlp1_2] += (occ[1] - occ[2]).imag();
                            this->locale[iat][target_l][0][0].c[ind_m1m2 + 3 * tlp1_2] += (occ[0] - occ[3]).real();
                            ind_m1m2++;
                        }
                    }
                }// ib
                begin_ih += nh;
            }// iat
        }// ik
    }
#if defined(__CUDA) || defined(__ROCM)
    else
    {
        auto* onsite_p = projectors::OnsiteProjector<double, base_device::DEVICE_GPU>::get_instance();
        const psi::Psi<std::complex<double>, base_device::DEVICE_GPU>* psi_p = (const psi::Psi<std::complex<double>, base_device::DEVICE_GPU>*)psi_in;
        // loop over k-points to calculate Mi of \sum_{k,i,l,m}<Psi_{k,i}|alpha_{l,m}><alpha_{l,m}|Psi_{k,i}>
        const int nbands = psi_p->get_nbands();
        for(int ik = 0; ik < psi_p->get_nk(); ik++)
        {
            psi_p->fix_k(ik);
            onsite_p->tabulate_atomic(ik);

            onsite_p->overlap_proj_psi(nbands*psi_p->get_npol(), psi_p->get_pointer());
            const std::complex<double>* becp = onsite_p->get_h_becp();
            // becp(nbands*npol , nkb)
            // mag = wg * \sum_{nh}becp * becp
            int nkb = onsite_p->get_size_becp() / nbands / psi_p->get_npol();
            int begin_ih = 0;
            for(int iat = 0; iat < cell.nat; iat++)
            {
                const int it = cell.iat2it[iat];
                const int nh = onsite_p->get_nh(iat);
                const int target_l = this->orbital_corr[it];
                if(target_l == -1)
                {
                    begin_ih += nh;
                    continue;
                }
                // m = l^2, l^2+1, ..., (l+1)^2-1
                const int m_begin = target_l * target_l;
                const int tlp1 = 2 * target_l + 1;
                const int tlp1_2 = tlp1 * tlp1;
                for(int ib = 0;ib<nbands;ib++)
                {
                    const double weight = wg_in(ik, ib);
                    int ind_m1m2 = 0;
                    for(int m1 = 0; m1 < tlp1; m1++)
                    {
                        const int index_m1 = ib*2*nkb + begin_ih + m_begin + m1;
                        for(int m2 = 0; m2 < tlp1; m2++)
                        {
                            const int index_m2 = ib*2*nkb + begin_ih + m_begin + m2;
                            std::complex<double> occ[4];
                            occ[0] = weight * conj(becp[index_m1]) * becp[index_m2];
                            occ[1] = weight * conj(becp[index_m1]) * becp[index_m2 + nkb];
                            occ[2] = weight * conj(becp[index_m1 + nkb]) * becp[index_m2];
                            occ[3] = weight * conj(becp[index_m1 + nkb]) * becp[index_m2 + nkb];
                            this->locale[iat][target_l][0][0].c[ind_m1m2] += (occ[0] + occ[3]).real();
                            this->locale[iat][target_l][0][0].c[ind_m1m2 + tlp1_2] += (occ[1] + occ[2]).real();
                            this->locale[iat][target_l][0][0].c[ind_m1m2 + 2 * tlp1_2] += (occ[1] - occ[2]).imag();
                            this->locale[iat][target_l][0][0].c[ind_m1m2 + 3 * tlp1_2] += (occ[0] - occ[3]).real();
                            ind_m1m2++;
                        }
                    }
                }// ib
                begin_ih += nh;
            }// iat
        }// ik
    }
#endif

    this->EU = 0.0;
    // reduce mag from all k-pools
    for(int iat = 0; iat < cell.nat; iat++)
    {
        const int it = cell.iat2it[iat];
        const int target_l = this->orbital_corr[it];
        if(target_l == -1)
        {
            continue;
        }
        const int size = (2 * target_l + 1) * (2 * target_l + 1);
        Parallel_Reduce::reduce_double_allpool(PARAM.inp.kpar, PARAM.globalv.nproc_in_pool, this->locale[iat][target_l][0][0].c, size * PARAM.inp.nspin);
        //update effective potential
        const double u_value = this->U[it];
        std::complex<double>* vu_iat = &(this->eff_pot_pw[this->eff_pot_pw_index[iat]]);
        const int m_size = 2 * target_l + 1;
        for (int m1 = 0; m1 < m_size; m1++)
        {
            for (int m2 = 0; m2 < m_size; m2++)
            {
                vu_iat[m1 * m_size + m2] = u_value * (1.0 * (m1 == m2) - this->locale[iat][target_l][0][0].c[m2 * m_size + m1]);
                this->EU += u_value * 0.25 * this->locale[iat][target_l][0][0].c[m2 * m_size + m1] * this->locale[iat][target_l][0][0].c[m1 * m_size + m2];
            }
        }
        for (int is = 1; is < 4; ++is)
        {
            int start = is * m_size * m_size;
            for (int m1 = 0; m1 < m_size; m1++)
            {
                for (int m2 = 0; m2 < m_size; m2++)
                {
                    vu_iat[start + m1 * m_size + m2] = u_value * (0 - this->locale[iat][target_l][0][0].c[start + m2 * m_size + m1]);
                    this->EU += u_value * 0.25 * this->locale[iat][target_l][0][0].c[start + m2 * m_size + m1] * this->locale[iat][target_l][0][0].c[start + m1 * m_size + m2];
                }
            }
        }
        // transfer from Pauli matrix representation to spin representation 
        for (int m1 = 0; m1 < m_size; m1++)
        {
            for (int m2 = 0; m2 < m_size; m2++)
            {
                int index[4];
                index[0] = m1 * m_size + m2;
                index[1] = m1 * m_size + m2 + size;
                index[2] = m1 * m_size + m2 + size * 2;
                index[3] = m1 * m_size + m2 + size * 3;
                std::complex<double> vu_tmp[4];
                for (int i = 0; i < 4; i++)
                {
                    vu_tmp[i] = vu_iat[index[i]];
                }
                vu_iat[index[0]] = 0.5 * (vu_tmp[0] + vu_tmp[3]);
                vu_iat[index[3]] = 0.5 * (vu_tmp[0] - vu_tmp[3]);
                vu_iat[index[1]] = 0.5 * (vu_tmp[1] + std::complex<double>(0.0, 1.0) * vu_tmp[2]);
                vu_iat[index[2]] = 0.5 * (vu_tmp[1] - std::complex<double>(0.0, 1.0) * vu_tmp[2]);
            }
        }
    }

    if(mixing_dftu && initialed_locale)
    {
        this->mix_locale(cell, mixing_beta);
    }
    // update effective potential
    ModuleBase::timer::tick("DFTU", "cal_occ_pw");
}
/// calculate the local DFT+U effective potential matrix for PW base.
void DFTU::cal_VU_pot_pw(const int spin)
{

}

} // namespace ModuleDFTU