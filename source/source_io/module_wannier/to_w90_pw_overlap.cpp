#include "to_w90_pw.h"

#include "source_io/module_parameter/parameter.h"
#include "source_base/parallel_reduce.h"

void toW90_PW::unkdotkb(
    const psi::Psi<std::complex<double>>& psi_pw,
    const ModulePW::PW_Basis_K* wfcpw,
    const int& cal_ik,
    const int& cal_ikb,
    const ModuleBase::Vector3<double> G,
    ModuleBase::ComplexMatrix &Mmn
)
{
    Mmn.create(num_bands, num_bands);

    for (int m = 0; m < num_bands; m++)
    {
        int im = cal_band_index[m];
        std::vector<std::complex<double>> phase(wfcpw->nmaxgr);

        // get the phase value in realspace
        for (int ig = 0; ig < wfcpw->npwk[cal_ik]; ig++)
        {
            // if wfcpw->getgdirect(cal_ik, ig) == phase_G
            ModuleBase::Vector3<double> temp_G = wfcpw->getgdirect(cal_ik, ig) - G;
            if (temp_G.norm() < 1e-9)
            {
                phase[ig] = std::complex<double>(1.0, 0.0);
                break;
            }
        }

        wfcpw->recip2real(phase.data(), phase.data(), cal_ik);

        if (nspin_ == 4)
        {
            // (1) set value
            std::vector<std::complex<double>> psir_up(wfcpw->nmaxgr);
            std::vector<std::complex<double>> psir_dn(wfcpw->nmaxgr);

            // (2) fft and get value
            // int npw_ik = wfcpw->npwk[cal_ik];
            int npwx = wfcpw->npwk_max;
            wfcpw->recip2real(&psi_pw(cal_ik, im, 0), psir_up.data(), cal_ik);
            // wfcpw->recip2real(&psi_pw(cal_ik, im, npw_ik), psir_dn, cal_ik);
            wfcpw->recip2real(&psi_pw(cal_ik, im, npwx), psir_dn.data(), cal_ik);
            for (int ir = 0; ir < wfcpw->nrxx; ir++)
            {
                psir_up[ir] *= phase[ir];
                psir_dn[ir] *= phase[ir];
            }

            wfcpw->real2recip(psir_up.data(), psir_up.data(), cal_ikb);
            wfcpw->real2recip(psir_dn.data(), psir_dn.data(), cal_ikb);

            for (int n = 0; n < num_bands; n++)
            {
                int in = cal_band_index[n];

                if (!gamma_only_wannier)
                {
                    std::complex<double> result_tem(0.0, 0.0);

                    // int npw_ikb = wfcpw->npwk[cal_ikb];
                    int pwNumberMax = wfcpw->npwk_max;

                    // Can be accelerated using lapack
                    for (int ig = 0; ig < pwNumberMax; ig++)
                    {
                        result_tem = result_tem + conj(psir_up[ig]) * psi_pw(cal_ikb, in, ig)
                                                   + conj(psir_dn[ig]) * psi_pw(cal_ikb, in, ig+pwNumberMax);
                    }
                    Parallel_Reduce::reduce_all(result_tem);
                    Mmn(m, n) = result_tem;
                }
                else
                {
                    // GlobalV::ofs_running << "gamma only test" << std::endl;
                }

            }

        }
        else
        {
            // (1) set value
            std::vector<std::complex<double>> psir(wfcpw->nmaxgr);

            // (2) fft and get value
            wfcpw->recip2real(&psi_pw(cal_ik, im, 0), psir.data(), cal_ik);
            for (int ir = 0; ir < wfcpw->nrxx; ir++)
            {
                psir[ir] *= phase[ir];
            }

            wfcpw->real2recip(psir.data(), psir.data(), cal_ikb);

            for (int n = 0; n < num_bands; n++)
            {
                int in = cal_band_index[n];

                if (!gamma_only_wannier)
                {
                    std::complex<double> result_tem(0.0, 0.0);

                    // Can be accelerated using lapack
                    for (int ig = 0; ig < wfcpw->npwk[cal_ikb]; ig++)
                    {
                        result_tem = result_tem + conj(psir[ig]) * psi_pw(cal_ikb, in, ig);
                    }

                    Parallel_Reduce::reduce_all(result_tem);
                    Mmn(m, n) = result_tem;
                }
                else
                {
                    // GlobalV::ofs_running << "gamma only test" << std::endl;
                }

            }

        }

    }

}

void toW90_PW::unkdotW_A(
    const psi::Psi<std::complex<double>>& psi_pw,
    const ModulePW::PW_Basis_K* wfcpw,
    const int& ik,
    const std::vector<ModuleBase::matrix> &radial_in_q,
    ModuleBase::ComplexMatrix &Amn
)
{
    Amn.create(num_bands, num_wannier);

    int npw = wfcpw->npwk[ik];
    int npwx = wfcpw->npwk_max;
    ModuleBase::ComplexMatrix trial_orbitals;
    produce_trial_in_pw(psi_pw, ik, wfcpw, radial_in_q, trial_orbitals);

    for (int iw = 0; iw < num_wannier; iw++)
    {
        for (int ib_w = 0; ib_w < num_bands; ib_w++)
        {
            int ib = cal_band_index[ib_w];

            if (nspin_ != 4)
            {
                for (int ig = 0; ig < npw; ig++)
                {
                    Amn(ib_w, iw) += conj(psi_pw(ik, ib, ig)) * trial_orbitals(iw, ig);
                }
            }
            else
            {
                for (int ig = 0; ig < npw; ig++)
                {
                    Amn(ib_w, iw) += up_con[iw] * conj(psi_pw(ik, ib, ig)) * trial_orbitals(iw, ig)
                                   + dn_con[iw] * conj(psi_pw(ik, ib, ig+npwx)) * trial_orbitals(iw, ig);
                }
            }
        }
    }

    Parallel_Reduce::reduce_all(Amn.c, Amn.size);

}
