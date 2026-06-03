#include "dftu.h"
#include "source_pw/module_pwdft/onsite_proj.h"
#include "source_base/parallel_reduce.h"
#include "source_io/module_parameter/parameter.h"
#include "source_base/timer.h"
#include "source_base/parallel_global.h"

/// calculate occupation matrix for DFT+U (PW basis)
///
/// nspin=1 (npol=1): single spin channel; locale[iat][l][n][0] only;
///   eff_pot_pw has one block of tlp1^2 per atom.
///
/// nspin=2 (npol=1): two spin channels stored separately:
///   locale[iat][l][n][0] = spin-up, locale[iat][l][n][1] = spin-down;
///   becp indices: ib*nkb + begin_ih + m (same formula for both spins);
///   spin channel selected by `isk[ik]` (not ik >= nk/2, which fails for kpar>1);
///
/// nspin=4 (npol=2): spinor calculation;
///   locale has a single matrix of size (2*tlp1) x (2*tlp1) per atom
///   storing all 4 Pauli blocks contiguously.
void Plus_U::cal_occ_pw(const int iter, 
		const void* psi_in, 
		const ModuleBase::matrix& wg_in, 
		const UnitCell& cell, 
		Charge_Mixing* p_chgmix,
		const int* isk)
{
    ModuleBase::timer::start("Plus_U", "cal_occ_pw");
    this->copy_locale(cell);
    this->zero_locale(cell);

    if(PARAM.inp.device == "cpu")
    {
        auto* onsite_p = projectors::OnsiteProjector<double, base_device::DEVICE_CPU>::get_instance();
        const psi::Psi<std::complex<double>>* psi_p = (const psi::Psi<std::complex<double>>*)psi_in;
        const int nbands = psi_p->get_nbands();
        const int npol = psi_p->get_npol();
        for(int ik = 0; ik < psi_p->get_nk(); ik++)
        {
            int is = (PARAM.inp.nspin == 2) ? isk[ik] : 0;
            psi_p->fix_k(ik);
            onsite_p->tabulate_atomic(ik);

            onsite_p->overlap_proj_psi(nbands*npol, psi_p->get_pointer());
            const std::complex<double>* becp = onsite_p->get_h_becp();
            int nkb = onsite_p->get_size_becp() / nbands / npol;

            int begin_ih = 0;
            for(int iat = 0; iat < cell.nat; iat++)
            {
                const int it = cell.iat2it[iat];
                const int nh = onsite_p->get_nh(iat);
                const int target_l = get_orbital_corr(it);
                if(!has_correlated_orbital(it))
                {
                    begin_ih += nh;
                    continue;
                }
                const int m_begin = target_l * target_l;
                const int tlp1 = 2 * target_l + 1;
                const int tlp1_2 = tlp1 * tlp1;
                if(PARAM.inp.nspin == 4)
                {
                    for(int ib = 0;ib<nbands;ib++)
                    {
                        const double weight = wg_in(ik, ib);
                        int ind_m1m2 = 0;
                        for(int m1 = 0; m1 < tlp1; m1++)
                        {
                            const int index_m1 = ib*npol*nkb + begin_ih + m_begin + m1;
                            for(int m2 = 0; m2 < tlp1; m2++)
                            {
                                const int index_m2 = ib*npol*nkb + begin_ih + m_begin + m2;
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
                }
                else // nspin=1 or nspin=2
                {
                    for(int ib = 0;ib<nbands;ib++)
                    {
                        const double weight = wg_in(ik, ib);
                        int ind_m1m2 = 0;
                        for(int m1 = 0; m1 < tlp1; m1++)
                        {
                            const int index_m1 = ib*nkb + begin_ih + m_begin + m1;
                            for(int m2 = 0; m2 < tlp1; m2++)
                            {
                                const int index_m2 = ib*nkb + begin_ih + m_begin + m2;
                                this->locale[iat][target_l][0][is].c[ind_m1m2] += weight * (conj(becp[index_m1]) * becp[index_m2]).real();
                                ind_m1m2++;
                            }
                        }
                    }// ib
                }
                begin_ih += nh;
            }// iat

        }// ik
    }
#if defined(__CUDA) || defined(__ROCM)
    else
    {
        auto* onsite_p = projectors::OnsiteProjector<double, base_device::DEVICE_GPU>::get_instance();
        const psi::Psi<std::complex<double>, base_device::DEVICE_GPU>* psi_p = (const psi::Psi<std::complex<double>, base_device::DEVICE_GPU>*)psi_in;
        const int nbands = psi_p->get_nbands();
        const int npol = psi_p->get_npol();
        for(int ik = 0; ik < psi_p->get_nk(); ik++)
        {
            int is = (PARAM.inp.nspin == 2) ? isk[ik] : 0;
            psi_p->fix_k(ik);
            onsite_p->tabulate_atomic(ik);

            onsite_p->overlap_proj_psi(nbands*npol, psi_p->get_pointer());
            const std::complex<double>* becp = onsite_p->get_h_becp();
            int nkb = onsite_p->get_size_becp() / nbands / npol;
            int begin_ih = 0;
            for(int iat = 0; iat < cell.nat; iat++)
            {
                const int it = cell.iat2it[iat];
                const int nh = onsite_p->get_nh(iat);
                const int target_l = get_orbital_corr(it);
                if(!has_correlated_orbital(it))
                {
                    begin_ih += nh;
                    continue;
                }
                const int m_begin = target_l * target_l;
                const int tlp1 = 2 * target_l + 1;
                const int tlp1_2 = tlp1 * tlp1;
                if(PARAM.inp.nspin == 4)
                {
                    for(int ib = 0;ib<nbands;ib++)
                    {
                        const double weight = wg_in(ik, ib);
                        int ind_m1m2 = 0;
                        for(int m1 = 0; m1 < tlp1; m1++)
                        {
                            const int index_m1 = ib*npol*nkb + begin_ih + m_begin + m1;
                            for(int m2 = 0; m2 < tlp1; m2++)
                            {
                                const int index_m2 = ib*npol*nkb + begin_ih + m_begin + m2;
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
                }
                else // nspin=1 or nspin=2
                {
                    for(int ib = 0;ib<nbands;ib++)
                    {
                        const double weight = wg_in(ik, ib);
                        int ind_m1m2 = 0;
                        for(int m1 = 0; m1 < tlp1; m1++)
                        {
                            const int index_m1 = ib*nkb + begin_ih + m_begin + m1;
                            for(int m2 = 0; m2 < tlp1; m2++)
                            {
                                const int index_m2 = ib*nkb + begin_ih + m_begin + m2;
                                this->locale[iat][target_l][0][is].c[ind_m1m2] += weight * (conj(becp[index_m1]) * becp[index_m2]).real();
                                ind_m1m2++;
                            }
                        }
                    }// ib
                }
                begin_ih += nh;
            }// iat
        }// ik
    }
#endif

    // reduce locale from all k-pools
    for(int iat = 0; iat < cell.nat; iat++)
    {
        const int it = cell.iat2it[iat];
        const int target_l = get_orbital_corr(it);
        if(!has_correlated_orbital(it))
        {
            continue;
        }
        const int size = (2 * target_l + 1) * (2 * target_l + 1);

        if(PARAM.inp.nspin != 4)
        {
            Parallel_Reduce::reduce_double_allpool(PARAM.inp.kpar, 
                    GlobalV::NPROC_IN_POOL, 
                    this->locale[iat][target_l][0][0].c, 
                    size);
            if(PARAM.inp.nspin == 2)
            {
                Parallel_Reduce::reduce_double_allpool(PARAM.inp.kpar, 
                        GlobalV::NPROC_IN_POOL, 
                        this->locale[iat][target_l][0][1].c, 
                        size);
            }
        }
        else
        {
            Parallel_Reduce::reduce_double_allpool(PARAM.inp.kpar, 
                    GlobalV::NPROC_IN_POOL, 
                    this->locale[iat][target_l][0][0].c, 
                    size * 4);
        }

        // save locale matrix for this iat to uom_array
        if(this->uom_array.size() != 0)
        {
            for(int mm=0;mm<size;mm++)
            {
                this->uom_array[eff_pot_pw_index[iat]+mm] = this->locale[iat][target_l][0][0].c[mm];
            }
            if(PARAM.inp.nspin == 2)
            {
                const int half_size = this->uom_array.size() / 2;
                for(int mm=0;mm<size;mm++)
                {
                    this->uom_array[half_size + eff_pot_pw_index[iat]+mm] = this->locale[iat][target_l][0][1].c[mm];
                }
            }
        }
    }

    // mixing
    if(is_mixing_enabled() && p_chgmix != nullptr)
    {
        p_chgmix->mix_uom(this->uom_array, this->uom_save);
        this->set_locale(cell);
    }

    Plus_U::energy_u = 0.0;
    const double weight_eu = (PARAM.inp.nspin == 1) ? 0.25 : (PARAM.inp.nspin == 2) ? 0.5 : 0.25;
    const double diag_coeff = (PARAM.inp.nspin == 4) ? 1.0 : 0.5;
    // calculate VU and energy (locale already reduced above)
    for(int iat = 0; iat < cell.nat; iat++)
    {
        const int it = cell.iat2it[iat];
        const int target_l = get_orbital_corr(it);
        if(!has_correlated_orbital(it))
        {
            continue;
        }
        const int size = (2 * target_l + 1) * (2 * target_l + 1);

        //update effective potential
        const double u_value = this->U[it];
        std::complex<double>* vu_iat = &(this->eff_pot_pw[this->eff_pot_pw_index[iat]]);
        const int m_size = 2 * target_l + 1;

        if(PARAM.inp.nspin == 4)
        {
            for (int m1 = 0; m1 < m_size; m1++)
            {
                for (int m2 = 0; m2 < m_size; m2++)
                {
                    vu_iat[m1 * m_size + m2] = u_value * 
                      (diag_coeff * (m1 == m2) - this->locale[iat][target_l][0][0].c[m2 * m_size + m1]);
                    Plus_U::energy_u += u_value * weight_eu * this->locale[iat][target_l][0][0].c[m2 * m_size + m1] 
                             * this->locale[iat][target_l][0][0].c[m1 * m_size + m2];
                }
            }
            for (int is = 1; is < 4; ++is)
            {
                int start = is * m_size * m_size;
                for (int m1 = 0; m1 < m_size; m1++)
                {
                    for (int m2 = 0; m2 < m_size; m2++)
                    {
                        vu_iat[start + m1 * m_size + m2] = u_value * 
                          (0 - this->locale[iat][target_l][0][0].c[start + m2 * m_size + m1]);
                        Plus_U::energy_u += u_value * weight_eu 
                                 * this->locale[iat][target_l][0][0].c[start + m2 * m_size + m1] 
                                 * this->locale[iat][target_l][0][0].c[start + m1 * m_size + m2];
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
        else // nspin=1 or nspin=2
        {
            // spin-up channel
            for (int m1 = 0; m1 < m_size; m1++)
            {
                for (int m2 = 0; m2 < m_size; m2++)
                {
                    vu_iat[m1 * m_size + m2] = u_value * 
                      (diag_coeff * (m1 == m2) - this->locale[iat][target_l][0][0].c[m2 * m_size + m1]);
                    Plus_U::energy_u += u_value * weight_eu * this->locale[iat][target_l][0][0].c[m2 * m_size + m1] 
                             * this->locale[iat][target_l][0][0].c[m1 * m_size + m2];
                }
            }
            // spin-down channel for nspin=2
            if(PARAM.inp.nspin == 2)
            {
                std::complex<double>* vu_iat1 = &(this->eff_pot_pw[this->eff_pot_pw.size()/2 + this->eff_pot_pw_index[iat]]);
                for (int m1 = 0; m1 < m_size; m1++)
                {
                    for (int m2 = 0; m2 < m_size; m2++)
                    {
                        vu_iat1[m1 * m_size + m2] = u_value * 
                          (diag_coeff * (m1 == m2) - this->locale[iat][target_l][0][1].c[m2 * m_size + m1]);
                        Plus_U::energy_u += u_value * weight_eu * this->locale[iat][target_l][0][1].c[m2 * m_size + m1] 
                                 * this->locale[iat][target_l][0][1].c[m1 * m_size + m2];
                    }
                }
            }
        }
    }

    ModuleBase::timer::end("Plus_U", "cal_occ_pw");
}

