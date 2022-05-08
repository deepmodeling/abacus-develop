#include "elecstate_pw.h"

#include "module_base/constants.h"
#include "src_parallel/parallel_reduce.h"
#include "src_pw/global.h"

namespace elecstate
{

void ElecStatePW::psiToRho(const psi::Psi<std::complex<double>>& psi)
{
    this->calculate_weights();
    for (int ik = 0; ik < psi.get_nk(); ++ik)
    {
        psi.fix_k(ik);
        this->updateRhoK(psi);
    }
    this->parallelK();
}

void ElecStatePW::updateRhoK(const psi::Psi<std::complex<double>>& psi)
{
    const int current_k = psi.get_current_k();
    if (current_k == 0)
    {
        this->eband = 0.0;
    }
    this->eband += this->eBandK(current_k);

    this->rhoBandK(psi);

    return;
}

/*void ElecStatePW::getNewRho()
{
    return;
}*/

void ElecStatePW::parallelK()
{
#ifdef __MPI
    charge->rho_mpi();
    if (GlobalV::CALCULATION != "scf-sto" && GlobalV::CALCULATION != "relax-sto"
        && GlobalV::CALCULATION != "md-sto") // qinarui add it temporarily.
    {
        //==================================
        // Reduce all the Energy in each cpu
        //==================================
        this->eband /= GlobalV::NPROC_IN_POOL;
        Parallel_Reduce::reduce_double_all(this->eband);
    }
#endif
    return;
}

void ElecStatePW::rhoBandK(const psi::Psi<std::complex<double>>& psi)
{
    ModuleBase::TITLE("ElecStatePW", "rhoBandK");

    // used for plane wavefunction FFT3D to real space
    static std::vector<std::complex<double>> wfcr;
    wfcr.resize(this->charge->nrxx);
    // used for plane wavefunction FFT3D to real space, non-collinear spin case
    static std::vector<std::complex<double>> wfcr_another_spin;
    if (GlobalV::NSPIN == 4)
        wfcr_another_spin.resize(this->charge->nrxx);

    int ik = psi.get_current_k();
    int npw = psi.get_current_nbas();
    int current_spin = 0;
    if (GlobalV::NSPIN == 2)
    {
        current_spin = this->basis->Klist->isk[ik];
    }
    int nbands = psi.get_nbands();
    int* igk = &(GlobalC::wf.igk(ik, 0));

    //  here we compute the band energy: the sum of the eigenvalues
    if (GlobalV::NSPIN == 4)
    {
        int npwx = npw / 2;
        npw = this->basis->Klist->ngk[ik];
        for (int ibnd = 0; ibnd < nbands; ibnd++)
        {
            ///
            /// only occupied band should be calculated.
            ///
            if (this->wg(ik, ibnd) < ModuleBase::threshold_wg)
                continue;
            ModuleBase::GlobalFunc::ZEROS(wfcr.data(), this->charge->nrxx);
            for (int ig = 0; ig < npw; ig++)
            {
                wfcr[this->basis->ig2fftw[igk[ig]]] = psi(ibnd, ig);
            }
            const_cast<PW_Basis*>(this->basis)->FFT_wfc.FFT3D(wfcr.data(), 1);

            ModuleBase::GlobalFunc::ZEROS(wfcr_another_spin.data(), this->charge->nrxx);
            for (int ig = 0; ig < npw; ig++)
            {
                wfcr_another_spin[this->basis->ig2fftw[igk[ig]]] = psi(ibnd, ig + npwx);
            }
            const_cast<PW_Basis*>(this->basis)->FFT_wfc.FFT3D(wfcr_another_spin.data(), 1);

            const double w1 = this->wg(ik, ibnd) / GlobalC::ucell.omega;

            // Increment the charge density in chr.rho for real space
            for (int ir = 0; ir < this->charge->nrxx; ir++)
            {
                this->charge->rho[0][ir] += w1 * (norm(wfcr[ir]) + norm(wfcr_another_spin[ir]));
            }
            // In this case, calculate the three components of the magnetization
            if (GlobalV::DOMAG)
            {
                for (int ir = 0; ir < this->charge->nrxx; ir++)
                {
                    this->charge->rho[1][ir] += w1 * 2.0
                                                * (wfcr[ir].real() * wfcr_another_spin[ir].real()
                                                   + wfcr[ir].imag() * wfcr_another_spin[ir].imag());
                    this->charge->rho[2][ir] += w1 * 2.0
                                                * (wfcr[ir].real() * wfcr_another_spin[ir].imag()
                                                   - wfcr_another_spin[ir].real() * wfcr[ir].imag());
                    this->charge->rho[3][ir] += w1 * (norm(wfcr[ir]) - norm(wfcr_another_spin[ir]));
                }
            }
            else if (GlobalV::DOMAG_Z)
            {
                for (int ir = 0; ir < this->charge->nrxx; ir++)
                {
                    this->charge->rho[1][ir] = 0;
                    this->charge->rho[2][ir] = 0;
                    this->charge->rho[3][ir] += w1 * (norm(wfcr[ir]) - norm(wfcr_another_spin[ir]));
                }
            }
            else
                for (int is = 1; is < 4; is++)
                {
                    for (int ir = 0; ir < this->charge->nrxx; ir++)
                        this->charge->rho[is][ir] = 0;
                }
        }
    }
    else
        for (int ibnd = 0; ibnd < nbands; ibnd++)
        {
            ///
            /// only occupied band should be calculated.
            ///
            if (this->wg(ik, ibnd) < ModuleBase::threshold_wg)
                continue;

            ModuleBase::GlobalFunc::ZEROS(wfcr.data(), this->charge->nrxx);
            for (int ig = 0; ig < npw; ig++)
            {
                wfcr[this->basis->ig2fftw[igk[ig]]] = psi(ibnd, ig);
            }
            const_cast<PW_Basis*>(this->basis)->FFT_wfc.FFT3D(wfcr.data(), 1);

            const double w1 = this->wg(ik, ibnd) / GlobalC::ucell.omega;

            if (w1 != 0.0)
            {
                for (int ir = 0; ir < this->charge->nrxx; ir++)
                {
                    this->charge->rho[current_spin][ir] += w1 * norm(wfcr[ir]);
                }
            }

            // kinetic energy density
            if (XC_Functional::get_func_type() == 3)
            {
                for (int j = 0; j < 3; j++)
                {
                    ModuleBase::GlobalFunc::ZEROS(wfcr.data(), this->charge->nrxx);
                    for (int ig = 0; ig < npw; ig++)
                    {
                        double fact
                            = this->basis->get_GPlusK_cartesian_projection(ik, igk[ig], j) * GlobalC::ucell.tpiba;
                        wfcr[this->basis->ig2fftw[igk[ig]]] = psi(ibnd, ig) * complex<double>(0.0, fact);
                    }
                    const_cast<PW_Basis*>(this->basis)->FFT_wfc.FFT3D(wfcr.data(), 1);
                    for (int ir = 0; ir < this->charge->nrxx; ir++)
                    {
                        this->charge->kin_r[current_spin][ir] += w1 * norm(wfcr[ir]);
                    }
                }
            }
        }

    return;
}

} // namespace elecstate