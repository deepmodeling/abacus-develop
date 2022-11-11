#include "elecstate_pw.h"

#include "module_base/constants.h"
#include "src_parallel/parallel_reduce.h"
#include "src_pw/global.h"
#include "module_base/timer.h"

namespace elecstate
{

void ElecStatePW::psiToRho(const psi::Psi<std::complex<double>>& psi)
{
    ModuleBase::TITLE("ElecStatePW", "psiToRho");
    ModuleBase::timer::tick("ElecStatePW", "psiToRho");
    this->calculate_weights();

    this->calEBand();

    for(int is=0; is<GlobalV::NSPIN; is++)
	{
		ModuleBase::GlobalFunc::ZEROS(this->charge->rho[is], this->charge->nrxx);
		if (XC_Functional::get_func_type() == 3)
		{
            ModuleBase::GlobalFunc::ZEROS(this->charge->kin_r[is], this->charge->nrxx);
        }
	}

    for (int ik = 0; ik < psi.get_nk(); ++ik)
    {
        psi.fix_k(ik);
        this->updateRhoK(psi);
    }
    this->parallelK();
    ModuleBase::timer::tick("ElecStatePW", "psiToRho");
}

void ElecStatePW::updateRhoK(const psi::Psi<std::complex<double>>& psi)
{
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
    if(GlobalV::ESOLVER_TYPE == "sdft") //qinarui add it 2021-7-21
	{
		this->eband /= GlobalV::NPROC_IN_POOL;
		MPI_Allreduce(MPI_IN_PLACE, &this->eband, 1, MPI_DOUBLE, MPI_SUM , STO_WORLD);
	}
#endif
    return;
}

void ElecStatePW::rhoBandK(const psi::Psi<std::complex<double>>& psi)
{
    ModuleBase::TITLE("ElecStatePW", "rhoBandK");

    // used for plane wavefunction FFT3D to real space
    static std::vector<std::complex<double>> wfcr;
    wfcr.resize(this->basis->nmaxgr);
    // used for plane wavefunction FFT3D to real space, non-collinear spin case
    static std::vector<std::complex<double>> wfcr_another_spin;
    if (GlobalV::NSPIN == 4)
        wfcr_another_spin.resize(this->charge->nrxx);

    int ik = psi.get_current_k();
    int npw = psi.get_current_nbas();
    int current_spin = 0;
    if (GlobalV::NSPIN == 2)
    {
        current_spin = this->klist->isk[ik];
    }
    int nbands = psi.get_nbands();
    //  here we compute the band energy: the sum of the eigenvalues
    if (GlobalV::NSPIN == 4)
    {
        int npwx = npw / 2;
        npw = this->klist->ngk[ik];
        for (int ibnd = 0; ibnd < nbands; ibnd++)
        {
            ///
            /// only occupied band should be calculated.
            ///
            if (this->wg(ik, ibnd) < ModuleBase::threshold_wg)
                continue;


            this->basis->recip2real(&psi(ibnd,0), wfcr.data(), ik);

            this->basis->recip2real(&psi(ibnd,npwx), wfcr_another_spin.data(), ik);

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
    {
        for (int ibnd = 0; ibnd < nbands; ibnd++)
        {
            ///
            /// only occupied band should be calculated.
            ///
            if (this->wg(ik, ibnd) < ModuleBase::threshold_wg)
                continue;

            this->basis->recip2real(&psi(ibnd,0), wfcr.data(), ik);

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
                            = this->basis->getgpluskcar(ik, ig)[j] * GlobalC::ucell.tpiba;
                        wfcr[ig] = psi(ibnd, ig) * complex<double>(0.0, fact);
                    }

                    this->basis->recip2real(wfcr.data(), wfcr.data(), ik);
                    
                    for (int ir = 0; ir < this->charge->nrxx; ir++)
                    {
                        this->charge->kin_r[current_spin][ir] += w1 * norm(wfcr[ir]);
                    }
                }
            }
        }
    }

    return;
}

} // namespace elecstate